"""Run one eval track: a Claude agent with the genomics MCP server's tools.

    uv run --extra eval python evals/run.py reasoning --limit 5          # pilot
    uv run --extra eval python evals/run.py tool_use --variant v1-prompt  # a named variant

Writes .claude/hillclimb/genomics-<track>/<variant>/{results.jsonl,errors.jsonl,traces/}
in the layout the claude-api skill's report builder reads. Resumable: rerunning
skips (case, rep) pairs already in results.jsonl.

The server runs as one stdio subprocess per track (shared by all concurrent
cases, so upstream rate limits are paced in one place) against the frozen
cache that build_dataset.py wrote. The reasoning track starts it with
GENOMICS_MCP_HIDE_CLINVAR=1.
"""

import argparse
import asyncio
import json
import math
import os
import random
import sys
import time
from pathlib import Path
from typing import Any

import anthropic
from mcp import Client, StdioServerParameters

sys.path.insert(0, str(Path(__file__).parent))
from grade import GRADERS, METRICS, TIERS  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "evals" / "data"
MAX_TURNS = 30
CASE_TIMEOUT_S = 900
MAX_TOKENS = 16000

STRENGTHS = ["stand-alone", "very strong", "strong", "moderate", "supporting"]

SUBMIT_TOOLS = {
    "reasoning": {
        "name": "submit_answer",
        "description": "Submit your final classification. Call exactly once, after gathering evidence.",
        "input_schema": {
            "type": "object",
            "properties": {
                "classification": {"type": "string", "enum": TIERS},
                "criteria": {
                    "type": "array",
                    "description": "ACMG/AMP criteria applied, e.g. PM2 at supporting strength.",
                    "items": {
                        "type": "object",
                        "properties": {
                            "code": {"type": "string"},
                            "strength": {"type": "string", "enum": STRENGTHS},
                            "evidence": {"type": "string"},
                        },
                        "required": ["code", "strength", "evidence"],
                        "additionalProperties": False,
                    },
                },
                "not_evaluable": {"type": "string", "description": "Criteria that need data you do not have."},
                "summary": {"type": "string"},
            },
            "required": ["classification", "criteria", "not_evaluable", "summary"],
            "additionalProperties": False,
        },
    },
    "tool_use": {
        "name": "submit_answer",
        "description": "Submit your final answer. Call exactly once.",
        "input_schema": {
            "type": "object",
            "properties": {
                "answer": {"type": "string", "description": "The answer value only, no explanation. Use an empty string if input_invalid."},
                "input_invalid": {"type": "boolean", "description": "True if the question's variant cannot be resolved as given (e.g. reference allele mismatch)."},
            },
            "required": ["answer", "input_invalid"],
            "additionalProperties": False,
        },
    },
    "retrieval": {
        "name": "submit_answer",
        "description": "Submit what ClinVar reports. Call exactly once.",
        "input_schema": {
            "type": "object",
            "properties": {
                "clinvar_accession": {"type": "string", "description": "VCV accession, e.g. VCV000055627."},
                "classification": {"type": "string", "description": "ClinVar's aggregate germline classification, verbatim."},
                "review_stars": {"type": "integer", "description": "ClinVar review status as 0-4 stars."},
                "has_conflicting_submissions": {"type": "boolean", "description": "True if ClinVar's review status reports conflicting classifications."},
            },
            "required": ["clinvar_accession", "classification", "review_stars", "has_conflicting_submissions"],
            "additionalProperties": False,
        },
    },
}

SYSTEM = {
    "reasoning": (
        "You are a variant scientist classifying a germline variant with the ACMG/AMP framework, using only the "
        "genomics tools provided. ClinVar is not available: reason from primary evidence (consequence, in-silico "
        "predictions, population frequency, gene constraint). Apply criteria only at the strength the evidence "
        "supports, and do not assume case-level, segregation or functional data you were not given. When you are "
        "done, call submit_answer once."
    ),
    "tool_use": (
        "Answer the question using the genomics tools provided. Report values exactly as the tools return them, "
        "without rounding. If the variant in the question cannot be resolved as written, set input_invalid. "
        "Call submit_answer once with the answer."
    ),
    "retrieval": (
        "Look up what ClinVar currently reports for the variant using the genomics tools provided, and report it "
        "faithfully. Call submit_answer once."
    ),
}


def user_prompt(track: str, case: dict) -> str:
    inp = case["input"]
    if track == "reasoning":
        return f"Classify {inp['variant']} ({inp['gene']}) for {inp['condition']}."
    if track == "tool_use":
        return inp["question"]
    return f"What does ClinVar report for {inp['variant']}?"


def anthropic_tools(mcp_tools, track: str) -> list[dict]:
    tools = [{"name": t.name, "description": t.description or "", "input_schema": t.input_schema} for t in mcp_tools]
    tools.append({**SUBMIT_TOOLS[track], "strict": True})
    return tools


def _result_text(res) -> str:
    parts = [c.text for c in res.content if getattr(c, "type", "") == "text"]
    return "\n".join(parts) if parts else json.dumps(res.structured_content)


class Runner:
    def __init__(self, args: argparse.Namespace):
        self.args = args
        self.track = args.track
        self.out = ROOT / ".claude" / "hillclimb" / f"genomics-{self.track}" / args.variant
        (self.out / "traces").mkdir(parents=True, exist_ok=True)
        self.client = anthropic.AsyncAnthropic(max_retries=4)
        self.sem = asyncio.Semaphore(args.concurrency)
        self.write_lock = asyncio.Lock()

    def done_keys(self) -> set[tuple[str, int]]:
        path = self.out / "results.jsonl"
        if not path.exists():
            return set()
        return {(r["prompt_id"], r["rep"]) for r in map(json.loads, path.read_text().splitlines())}

    async def append(self, name: str, row: dict) -> None:
        async with self.write_lock:
            with open(self.out / name, "a") as f:
                f.write(json.dumps(row) + "\n")

    async def run_case(self, mcp: Client, tools: list[dict], system: str, case: dict, rep: int) -> None:
        async with self.sem:
            t0 = time.monotonic()
            trace: list[dict] = [{"role": "system", "content": system}]
            usage = {"input_tokens": 0, "output_tokens": 0, "cache_read_input_tokens": 0, "cache_creation_input_tokens": 0}
            state = {"tool_calls": 0, "model": None, "stop_reason": None}
            try:
                answer = await asyncio.wait_for(self._loop(mcp, tools, system, case, trace, usage, state), CASE_TIMEOUT_S)
            except Exception as e:  # noqa: BLE001 - every failure gets a class, never a silent zero
                failure = "timeout" if isinstance(e, TimeoutError) else (
                    "refusal" if isinstance(e, Refusal) else "harness_error")
                await self.append("errors.jsonl", {
                    "prompt_id": case["id"], "rep": rep, "failure_class": failure, "error": f"{type(e).__name__}: {e}"[:500],
                    "model": state["model"], "usage": usage, "latency_s": round(time.monotonic() - t0, 1),
                })
                print(f"  {case['id']} rep{rep}: {failure} ({e})", flush=True)
                return
            trace_path = self.out / "traces" / f"{case['id']}_rep{rep}.json"
            trace_path.write_text(json.dumps(trace, indent=1, default=str))
            grade = GRADERS[self.track](case["expected"], answer) if answer is not None else {
                m["id"]: 0.0 for m in METRICS[self.track]}
            row = {
                "prompt_id": case["id"], "rep": rep, "prompt": user_prompt(self.track, case), "tags": case["tags"],
                "status": "truncated" if state["stop_reason"] == "max_tokens" else "ok",
                "stop_reason": state["stop_reason"], "grade": grade, "answer": answer,
                "expected": case["expected"], "model": state["model"], "usage": usage,
                "tool_calls": state["tool_calls"], "latency_s": round(time.monotonic() - t0, 1),
                "meta": {**case.get("meta", {}), "no_answer": answer is None},
            }
            await self.append("results.jsonl", row)
            headline = METRICS[self.track][0]["id"]
            print(f"  {case['id']} rep{rep}: {headline}={grade[headline]:.0f} tools={state['tool_calls']} "
                  f"{row['latency_s']}s", flush=True)

    async def _loop(self, mcp: Client, tools, system, case, trace, usage, state) -> dict | None:
        prompt = user_prompt(self.track, case)
        trace.append({"role": "user", "content": prompt})
        messages: list[dict] = [{"role": "user", "content": prompt}]
        nudged = False
        for _ in range(MAX_TURNS):
            resp = await self.client.messages.create(
                model=self.args.model, max_tokens=MAX_TOKENS, system=system, tools=tools, messages=messages,
                output_config={"effort": self.args.effort},
                cache_control={"type": "ephemeral"},
            )
            state["model"], state["stop_reason"] = resp.model, resp.stop_reason
            if not resp.model.startswith(self.args.model):
                raise RuntimeError(f"served model {resp.model} != requested {self.args.model}")
            for k in usage:
                usage[k] += getattr(resp.usage, k, 0) or 0
            if resp.stop_reason == "refusal":
                raise Refusal(getattr(resp.stop_details, "category", None))
            messages.append({"role": "assistant", "content": resp.content})  # keep thinking blocks intact

            results, submitted = [], None
            for block in resp.content:
                if block.type == "text" and block.text.strip():
                    trace.append({"role": "assistant", "content": block.text})
                elif block.type == "tool_use":
                    trace.append({"role": "tool_call", "name": block.name, "content": json.dumps(block.input, indent=2)})
                    if block.name == "submit_answer":
                        submitted = block.input
                        results.append({"type": "tool_result", "tool_use_id": block.id, "content": "Recorded."})
                        continue
                    state["tool_calls"] += 1
                    res = await mcp.call_tool(block.name, block.input)
                    text = _result_text(res)
                    trace.append({"role": "tool_result", "content": text})
                    results.append({"type": "tool_result", "tool_use_id": block.id, "content": text, "is_error": bool(res.is_error)})
            if submitted is not None:
                return submitted
            if resp.stop_reason == "max_tokens":
                return None
            if results:
                messages.append({"role": "user", "content": results})
                continue
            if nudged:  # ended twice without submitting
                return None
            nudged = True
            messages.append({"role": "user", "content": "Please call submit_answer with your final answer."})
        return None


class Refusal(Exception):
    pass


def write_state(track: str) -> None:
    flow = ROOT / ".claude" / "hillclimb" / f"genomics-{track}"
    path = flow / "_state.json"
    if not path.exists():
        path.write_text(json.dumps({
            "flow": f"genomics-{track}",
            "metrics": METRICS[track],
            "perf_fields": [{"id": "latency_s", "label": "latency", "unit": "s"}, {"id": "tool_calls", "label": "tool calls"}],
            "prices": {"claude-opus-5-5": {"in": 4.0, "out": 20.0}, "claude-sonnet-5-5": {"in": 2.0, "out": 10.0},
                       "claude-haiku-5-5": {"in": 0.10, "out": 0.50}},
        }, indent=2))


def summarize(out: Path, track: str) -> None:
    rows = [json.loads(line) for line in (out / "results.jsonl").read_text().splitlines()] if (out / "results.jsonl").exists() else []
    errs = (out / "errors.jsonl").read_text().splitlines() if (out / "errors.jsonl").exists() else []
    ok = [r for r in rows if r["status"] == "ok"]
    print(f"\n{track}: {len(ok)} scored, {len(rows) - len(ok)} truncated, {len(errs)} failed attempts")
    for m in METRICS[track]:
        vals = [r["grade"][m["id"]] for r in ok]
        if vals:
            p = sum(vals) / len(vals)
            ci = 1.96 * math.sqrt(p * (1 - p) / len(vals))
            print(f"  {m['label']:<14} {p:6.1%}  ±{ci:.1%}  (n={len(vals)})")
    tok_in = sum(r["usage"]["input_tokens"] for r in rows)
    tok_out = sum(r["usage"]["output_tokens"] for r in rows)
    cr = sum(r["usage"]["cache_read_input_tokens"] for r in rows)
    cw = sum(r["usage"]["cache_creation_input_tokens"] for r in rows)
    print(f"  tokens: {tok_in:,} in, {cw:,} cache-write, {cr:,} cache-read, {tok_out:,} out")


async def main(args: argparse.Namespace) -> None:
    cases = [json.loads(line) for line in (DATA / f"{args.track}.jsonl").read_text().splitlines()]
    if args.ids:
        cases = [c for c in cases if c["id"] in set(args.ids.split(","))]
    if args.limit:
        # A stratified pilot: spread the limit across tags[0] groups.
        random.Random(0).shuffle(cases)
        by_group: dict[str, list] = {}
        for c in cases:
            by_group.setdefault(c["tags"][0], []).append(c)
        cases = [c for i in range(args.limit) for g in by_group.values() if i < len(g) for c in [g[i]]][: args.limit]

    runner = Runner(args)
    write_state(args.track)
    done = runner.done_keys()
    todo = [(c, r) for c in cases for r in range(args.reps) if (c["id"], r) not in done]
    print(f"{args.track}/{args.variant}: {len(todo)} attempts to run ({len(done)} already done)")

    env = {**os.environ, "GENOMICS_MCP_CACHE": str(DATA / "cache.sqlite"), "GENOMICS_MCP_CACHE_TTL_S": "inf"}
    if args.track == "reasoning":
        env["GENOMICS_MCP_HIDE_CLINVAR"] = "1"
    server = StdioServerParameters(command=sys.executable, args=["-m", "genomics_mcp.server"], env=env, cwd=str(ROOT))
    async with Client(server, read_timeout_seconds=300) as mcp:
        listed = (await mcp.list_tools()).tools
        if args.track == "reasoning":
            assert all(t.name != "get_clinvar_classification" for t in listed), "ClinVar must be hidden"
        tools = anthropic_tools(listed, args.track)
        system = SYSTEM[args.track] + "\n\n" + (mcp.instructions or "")
        if args.track == "reasoning":
            # Attach the criteria reference the way a host app attaches an MCP resource.
            acmg = (await mcp.read_resource("genomics://acmg-criteria")).contents[0].text
            system += "\n\n<resource uri=\"genomics://acmg-criteria\">\n" + acmg + "\n</resource>"
        await asyncio.gather(*(runner.run_case(mcp, tools, system, c, r) for c, r in todo))
    summarize(runner.out, args.track)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("track", choices=list(GRADERS))
    p.add_argument("--variant", default="baseline")
    p.add_argument("--model", default="claude-opus-5-5")
    p.add_argument("--effort", default="high", choices=["low", "medium", "high", "xhigh", "max"])
    p.add_argument("--reps", type=int, default=1)
    p.add_argument("--limit", type=int, default=0)
    p.add_argument("--ids", default="")
    p.add_argument("--concurrency", type=int, default=4)
    asyncio.run(main(p.parse_args()))
