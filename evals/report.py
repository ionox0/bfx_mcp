"""Detailed HTML report for the classification tracks (interpret / reasoning).

    uv run python evals/report.py                       # .claude/hillclimb/genomics-interpret
    uv run python evals/report.py --flow genomics-reasoning

Reads every variant directory under the flow (results.jsonl, config.json,
traces/) and writes <flow>/results.html: one static page, no scripts. Every
string from data files is HTML-escaped, since case text comes from ClinVar and
model output.
"""

import argparse
import html
import json
import math
import os
import statistics
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from grade import THREE_TIER, TIERS  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "evals" / "data"
ABBR = {"Pathogenic": "P", "Likely pathogenic": "LP", "Uncertain significance": "VUS", "Likely benign": "LB", "Benign": "B"}
# $/MTok: input, output, cache write (1.25x input), cache read
PRICES = {
    "claude-opus-5-5": (4.0, 20.0, 5.0, 0.20),
    "claude-sonnet-5-5": (2.0, 10.0, 2.5, 0.20),
    "claude-haiku-5-5": (0.10, 0.50, 0.125, 0.01),
}
FEATURE_TEXT = {
    "codon": "get_clinvar_codon_variants (PS1/PM5)",
    "pvs1": "PVS1 context in annotate_variant",
    "intronic": "intronic HGVS normalization",
}
LOF = {"frameshift_variant", "stop_gained", "splice_donor_variant", "splice_acceptor_variant", "start_lost"}
SPLICE_NEAR = {"splice_region_variant", "splice_donor_region_variant", "splice_donor_5th_base_variant",
               "splice_polypyrimidine_tract_variant", "intron_variant"}

e = html.escape


def consequence_group(term: str) -> str:
    if term == "missense_variant":
        return "missense"
    if term in LOF:
        return "LoF (frameshift / nonsense / canonical splice)"
    if term in SPLICE_NEAR:
        return "intronic / splice region"
    if term == "synonymous_variant":
        return "synonymous"
    return "other (UTR, in-frame, ...)"


def wilson(k: int, n: int) -> tuple[float, float]:
    if n == 0:
        return (0.0, 0.0)
    z, p = 1.96, k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def pct(k: int, n: int, ci: bool = True) -> str:
    if n == 0:
        return "–"
    lo, hi = wilson(k, n)
    s = f"<b>{k / n:.0%}</b> <span class=muted>{k}/{n}</span>"
    return s + (f"<br><span class=ci>{lo:.0%}–{hi:.0%}</span>" if ci else "")


def cost(rows: list[dict]) -> float:
    total = 0.0
    for r in rows:
        pin, pout, pcw, pcr = PRICES.get(r.get("model") or "", PRICES["claude-opus-5-5"])
        u = r["usage"]
        total += (u["input_tokens"] * pin + u["output_tokens"] * pout
                  + u["cache_creation_input_tokens"] * pcw + u["cache_read_input_tokens"] * pcr) / 1e6
    return total


class Variant:
    def __init__(self, path: Path):
        self.name = path.name
        self.path = path
        self.config = json.loads((path / "config.json").read_text()) if (path / "config.json").exists() else {}
        self.rows = [json.loads(l) for l in (path / "results.jsonl").read_text().splitlines()] if (path / "results.jsonl").exists() else []
        self.errors = (path / "errors.jsonl").read_text().splitlines() if (path / "errors.jsonl").exists() else []
        self.dataset = self.config.get("dataset") or ("reasoning" if not self.name.startswith("holdout") else "reasoning_holdout")
        feats = self.config.get("features", "all")
        self.features = (sorted(FEATURE_TEXT) if feats == "all" else [] if feats == "none"
                         else [f for f in feats.split(",") if f in FEATURE_TEXT])
        self.by_id = {r["prompt_id"]: r for r in self.rows}

    def trace(self, case_id: str) -> list[dict]:
        p = self.path / "traces" / f"{case_id}_rep0.json"
        return json.loads(p.read_text()) if p.exists() else []

    def stats(self) -> dict:
        rows = self.rows
        def tier_k(group):
            sub = [r for r in rows if THREE_TIER[r["expected"]["classification"]] == group]
            ok = sum(THREE_TIER.get((r["answer"] or {}).get("classification", ""), "") == group for r in sub)
            return ok, len(sub)
        return {
            "n": len(rows),
            "plp": tier_k("P/LP"), "blb": tier_k("B/LB"), "vus": tier_k("VUS"),
            "three": (sum(int(r["grade"]["three_tier"]) for r in rows), len(rows)),
            "exact": (sum(int(r["grade"]["exact"]) for r in rows), len(rows)),
            "flip": (sum(int(r["grade"]["severe_error"]) for r in rows), len(rows)),
            "selfc": (sum(int(r["grade"].get("self_consistent", 0)) for r in rows), len(rows)),
            "latency": statistics.median(r["latency_s"] for r in rows) if rows else 0,
            "tools": statistics.mean(r["tool_calls"] for r in rows) if rows else 0,
            "cost": cost(rows),
        }


# ----------------------------------------------------------------------------- sections


def section_overview(variants: list[Variant], cases_by_ds: dict[str, list[dict]]) -> str:
    out = ["<h2 id=how>How this eval works</h2>"]
    out.append(
        "<p>Each case gives the agent one germline variant (HGVS on a RefSeq transcript), its gene and a condition, "
        "through the server's <code>interpret_variant</code> MCP prompt, exactly as a user would invoke "
        "<code>/mcp__genomics__interpret_variant</code>. The agent can call the server's tools and read its resources, "
        "then must call <code>submit_answer</code> with a 5-tier classification and the ACMG/AMP criteria it applied. "
        "ClinVar's per-variant classification tool is hidden, because ClinVar is the source of the gold labels.</p>")
    out.append("<h3>Gold labels</h3><p>ClinVar germline classifications with review status <i>reviewed by expert panel</i> "
               "or <i>practice guideline</i> (3–4 stars), sampled round-robin across ClinGen VCEP genes, 12 per tier. "
               "Expert panels can use functional, segregation and case-level evidence that no public aggregate source provides, "
               "so a perfect score is not attainable from these tools. A second, disjoint <i>held-out</i> set (different seed, "
               "no shared variants) measures changes that were developed while looking at the main set's failures.</p>")
    out.append("<div class=grid2>")
    for ds, cases in cases_by_ds.items():
        tiers = Counter(c["tags"][0] for c in cases)
        groups = Counter(consequence_group(c["tags"][2]) for c in cases)
        years = Counter((c["meta"].get("date_last_evaluated") or "?")[:4] for c in cases)
        out.append(f"<div class=card><h4>{e(ds)} <span class=muted>({len(cases)} cases, {len({c['input']['gene'] for c in cases})} genes)</span></h4>")
        out.append("<table class=mini><tr><th>tier</th>" + "".join(f"<td>{ABBR[t]} {tiers[t]}</td>" for t in TIERS) + "</tr></table>")
        out.append("<table class=mini>" + "".join(f"<tr><th>{e(g)}</th><td>{n}</td></tr>" for g, n in groups.most_common()) + "</table>")
        span = sorted(y for y in years if y != "?")
        out.append(f"<p class=muted>Last evaluated by the panel: {e(span[0]) if span else '?'}–{e(span[-1]) if span else '?'}</p></div>")
    out.append("</div>")

    out.append("<h3>Variants compared</h3><table><tr><th>variant</th><th>cases</th><th>model / effort</th><th>server features</th><th>ClinVar</th></tr>")
    for v in variants:
        c = v.config
        feats = ", ".join(FEATURE_TEXT[f] for f in v.features) or "none (original server)"
        out.append(f"<tr><td><b>{e(v.name)}</b></td><td>{e(v.dataset)}</td><td>{e(c.get('model', 'claude-opus-5-5'))} / {e(c.get('effort', 'high'))}</td>"
                   f"<td>{e(feats)}</td><td>{'visible' if c.get('clinvar') else 'hidden'}</td></tr>")
    out.append("</table>")

    sample = next((v for v in variants if v.rows), None)
    if sample:
        tr = sample.trace(sample.rows[0]["prompt_id"])
        system = next((t["content"] for t in tr if t["role"] == "system"), "")
        user = next((t["content"] for t in tr if t["role"] == "user"), "")
        out.append("<details><summary>System prompt (host + MCP server instructions)</summary>"
                   f"<pre>{e(system)}</pre></details>")
        out.append(f"<details><summary>Example user message (from the <code>interpret_variant</code> prompt)</summary><pre>{e(user)}</pre></details>")

    out.append("""<h3>Metrics</h3><table>
<tr><th>P/LP sensitivity</th><td>Of variants the panel calls Pathogenic or Likely pathogenic, the share the agent puts in P/LP. A miss is a missed actionable finding. <b>Headline.</b></td></tr>
<tr><th>P&harr;B flips</th><td>Pathogenic-side gold called benign-side, or the reverse. The most harmful single error. <b>Headline</b>, lower is better.</td></tr>
<tr><th>B/LB, VUS correct</th><td>The same share for benign-side and VUS gold.</td></tr>
<tr><th>3-tier / 5-tier</th><td>Agreement after collapsing to P/LP, VUS, B/LB; and exact tier agreement.</td></tr>
<tr><th>Self-consistent</th><td>Reference-free: the submitted tier equals the tier implied by the agent's own criteria under the ACMG point system (Supporting 1, Moderate 2, Strong 4, Very strong 8; P &ge;10, LP 6–9, VUS 0–5, LB &minus;1 to &minus;6, B &le;&minus;7).</td></tr>
<tr><th>Intervals</th><td>95% Wilson score intervals. With 24 P/LP cases per set, differences under ~20 points are within noise.</td></tr>
</table>""")
    return "\n".join(out)


def section_headline(variants: list[Variant]) -> str:
    out = ["<h2 id=results>Results</h2><table class=headline><tr><th>variant</th><th>P/LP sensitivity</th><th>P&harr;B flips</th>"
           "<th>B/LB correct</th><th>VUS correct</th><th>3-tier</th><th>5-tier exact</th><th>self-consistent</th>"
           "<th>median latency</th><th>tool calls</th><th>cost</th><th>failed</th></tr>"]
    for v in variants:
        s = v.stats()
        out.append(
            f"<tr><td><b>{e(v.name)}</b><br><span class=muted>{e(v.dataset)}</span></td><td>{pct(*s['plp'])}</td><td>{pct(*s['flip'], ci=False)}</td>"
            f"<td>{pct(*s['blb'])}</td><td>{pct(*s['vus'])}</td><td>{pct(*s['three'])}</td><td>{pct(*s['exact'])}</td>"
            f"<td>{pct(*s['selfc'], ci=False)}</td><td>{s['latency']:.0f} s</td><td>{s['tools']:.1f}</td><td>${s['cost']:.2f}</td>"
            f"<td>{len(v.errors)}</td></tr>")
    out.append("</table>")
    return "\n".join(out)


def section_changes(pairs: list[tuple[Variant, Variant]], cases: dict[str, dict]) -> str:
    out = ["<h2 id=changes>What changed</h2>"]
    for base, new in pairs:
        fixed, broke, moved = [], [], []
        for cid, r0 in base.by_id.items():
            r1 = new.by_id.get(cid)
            if not r1:
                continue
            a0, a1 = (r0["answer"] or {}).get("classification"), (r1["answer"] or {}).get("classification")
            if r0["grade"]["three_tier"] < r1["grade"]["three_tier"]:
                fixed.append((cid, a0, a1))
            elif r0["grade"]["three_tier"] > r1["grade"]["three_tier"]:
                broke.append((cid, a0, a1))
            elif a0 != a1:
                moved.append((cid, a0, a1))
        out.append(f"<h3>{e(base.name)} &rarr; {e(new.name)}</h3><p><b class=good>{len(fixed)} fixed</b>, "
                   f"<b class=bad>{len(broke)} regressed</b>, {len(moved)} changed tier within the same 3-tier group.</p>")
        if fixed or broke or moved:
            out.append("<table><tr><th></th><th>gene</th><th>variant</th><th>consequence</th><th>gold</th><th>before</th><th>after</th></tr>")
            for label, items in (("fixed", fixed), ("regressed", broke), ("moved", moved)):
                for cid, a0, a1 in items:
                    c = cases[cid]
                    out.append(f"<tr class={label}><td>{label}</td><td>{e(c['input']['gene'])}</td><td><code>{e(c['input']['variant'])}</code></td>"
                               f"<td>{e(c['tags'][2])}</td><td>{e(ABBR.get(c['expected']['classification'], '?'))}</td>"
                               f"<td>{e(ABBR.get(a0 or '', '–'))}</td><td>{e(ABBR.get(a1 or '', '–'))}</td></tr>")
            out.append("</table>")
    return "\n".join(out)


def section_confusion(variants: list[Variant]) -> str:
    out = ["<h2 id=confusion>Confusion matrices</h2><p class=muted>Rows: expert-panel gold. Columns: agent. "
           "Green = correct 3-tier group; red = P&harr;B flip.</p><div class=grid3>"]
    for v in variants:
        cm = Counter((r["expected"]["classification"], (r["answer"] or {}).get("classification")) for r in v.rows)
        mx = max(cm.values(), default=1)
        out.append(f"<div class=card><h4>{e(v.name)}</h4><table class=cm><tr><th></th>" + "".join(f"<th>{ABBR[t]}</th>" for t in TIERS) + "</tr>")
        for g in TIERS:
            out.append(f"<tr><th>{ABBR[g]}</th>")
            for a in TIERS:
                n = cm[(g, a)]
                kind = "ok" if THREE_TIER[g] == THREE_TIER[a] else "flip" if {THREE_TIER[g], THREE_TIER[a]} == {"P/LP", "B/LB"} else "off"
                out.append(f"<td class='{kind}' style='--a:{n / mx:.2f}'>{n or ''}</td>")
            out.append("</tr>")
        out.append("</table></div>")
    out.append("</div>")
    return "\n".join(out)


def section_breakdown(variants: list[Variant], cases: dict[str, dict]) -> str:
    out = ["<h2 id=breakdown>By consequence and gene</h2><h3>Consequence class</h3>"
           "<table><tr><th>class</th>" + "".join(f"<th>{e(v.name)}<br><span class=muted>3-tier · P/LP sens.</span></th>" for v in variants) + "</tr>"]
    groups = sorted({consequence_group(cases[r["prompt_id"]]["tags"][2]) for v in variants for r in v.rows})
    for g in groups:
        out.append(f"<tr><th>{e(g)}</th>")
        for v in variants:
            rows = [r for r in v.rows if consequence_group(cases[r["prompt_id"]]["tags"][2]) == g]
            plp = [r for r in rows if THREE_TIER[r["expected"]["classification"]] == "P/LP"]
            k = sum(int(r["grade"]["three_tier"]) for r in rows)
            kp = sum(int(r["grade"]["three_tier"]) for r in plp)
            out.append(f"<td>{pct(k, len(rows), ci=False)} · {pct(kp, len(plp), ci=False) if plp else '–'}</td>")
        out.append("</tr>")
    out.append("</table><h3>Gene</h3><table><tr><th>gene</th>" + "".join(f"<th>{e(v.name)}</th>" for v in variants) + "</tr>")
    genes = sorted({cases[r["prompt_id"]]["input"]["gene"] for v in variants for r in v.rows})
    for gname in genes:
        out.append(f"<tr><th>{e(gname)}</th>")
        for v in variants:
            rows = [r for r in v.rows if cases[r["prompt_id"]]["input"]["gene"] == gname]
            k = sum(int(r["grade"]["three_tier"]) for r in rows)
            out.append(f"<td>{pct(k, len(rows), ci=False) if rows else ''}</td>")
        out.append("</tr>")
    out.append("</table>")
    return "\n".join(out)


def section_tools(variants: list[Variant]) -> str:
    calls: dict[str, Counter] = {}
    errs: dict[str, Counter] = {}
    for v in variants:
        c, er = Counter(), Counter()
        for r in v.rows:
            tr = v.trace(r["prompt_id"])
            for i, t in enumerate(tr):
                if t["role"] == "tool_call" and t["name"] != "submit_answer":
                    c[t["name"]] += 1
                    nxt = tr[i + 1] if i + 1 < len(tr) else {}
                    if nxt.get("role") == "tool_result" and str(nxt.get("content", "")).startswith("Error"):
                        er[t["name"]] += 1
        calls[v.name], errs[v.name] = c, er
    names = sorted({n for c in calls.values() for n in c})
    out = ["<h2 id=tools>Tool usage</h2><p class=muted>Total calls across all cases (errors in red).</p><table><tr><th>tool</th>"
           + "".join(f"<th>{e(v.name)}</th>" for v in variants) + "</tr>"]
    for n in names:
        out.append(f"<tr><th><code>{e(n)}</code></th>" + "".join(
            f"<td>{calls[v.name][n] or ''}{f' <span class=bad>({errs[v.name][n]} err)</span>' if errs[v.name][n] else ''}</td>" for v in variants) + "</tr>")
    out.append("</table>")
    return "\n".join(out)


def case_detail(v: Variant, r: dict, case: dict, out_dir: Path) -> str:
    a = r["answer"] or {}
    g = r["grade"]
    status = "ok" if g["three_tier"] else ("flip" if g["severe_error"] else "off")
    head = (f"<span class='pill {status}'>{e(ABBR.get(case['expected']['classification'], '?'))} &rarr; "
            f"{e(ABBR.get(a.get('classification', ''), '–'))}</span> <b>{e(case['input']['gene'])}</b> "
            f"<code>{e(case['input']['variant'])}</code> <span class=muted>{e(case['tags'][2])} · {r['tool_calls']} tools · {r['latency_s']:.0f}s</span>")
    body = [f"<p><b>Condition:</b> {e(case['input']['condition'])} · <b>Gold:</b> {e(case['expected']['classification'])} "
            f"({e(case['meta'].get('clinvar_accession', ''))}, {e(case['meta'].get('review_status', ''))}, evaluated {e(case['meta'].get('date_last_evaluated') or '?')})</p>"]
    crit = a.get("criteria") or []
    if crit:
        body.append("<table class=crit><tr><th>code</th><th>strength</th><th>evidence the agent cited</th></tr>" + "".join(
            f"<tr><td><b>{e(str(c.get('code', '')))}</b></td><td>{e(str(c.get('strength', '')))}</td><td>{e(str(c.get('evidence', '')))}</td></tr>" for c in crit) + "</table>")
    if a.get("not_evaluable"):
        body.append(f"<p><b>Not evaluable:</b> {e(a['not_evaluable'])}</p>")
    if a.get("summary"):
        body.append(f"<p><b>Summary:</b> {e(a['summary'])}</p>")
    steps = []
    tr = v.trace(r["prompt_id"])
    for i, t in enumerate(tr):
        if t["role"] == "tool_call" and t["name"] != "submit_answer":
            try:
                args = ", ".join(f"{k}={v_!s}" for k, v_ in json.loads(t["content"]).items())
            except (ValueError, AttributeError):
                args = t["content"]
            nxt = tr[i + 1] if i + 1 < len(tr) else {}
            err = nxt.get("role") == "tool_result" and str(nxt.get("content", "")).startswith("Error")
            steps.append(f"<li{' class=bad' if err else ''}><code>{e(t['name'])}</code>({e(args)})"
                         + (f" &rarr; {e(str(nxt['content'])[:160])}" if err else "") + "</li>")
    if steps:
        body.append("<p><b>Tool calls</b></p><ol class=steps>" + "".join(steps) + "</ol>")
    rel = os.path.relpath(v.path / "traces" / f"{r['prompt_id']}_rep0.json", out_dir)
    body.append(f"<p><a href='{e(str(rel))}'>Full transcript (JSON)</a></p>")
    return f"<details class=case><summary>{head}</summary>{''.join(body)}</details>"


def section_cases(variants: list[Variant], cases: dict[str, dict], out_dir: Path) -> str:
    out = ["<h2 id=cases>Every case</h2><p class=muted>Wrong 3-tier answers first, then by gold tier. Expand a case for "
           "the criteria the agent applied, the evidence it cited, and its tool calls.</p>"]
    for v in variants:
        rows = sorted(v.rows, key=lambda r: (r["grade"]["three_tier"], TIERS.index(r["expected"]["classification"])))
        wrong = sum(1 for r in rows if not r["grade"]["three_tier"])
        out.append(f"<details class=variant><summary><b>{e(v.name)}</b> <span class=muted>({len(rows)} cases, {wrong} wrong)</span></summary>")
        out.extend(case_detail(v, r, cases[r["prompt_id"]], out_dir) for r in rows)
        out.append("</details>")
    return "\n".join(out)


CSS = """
:root{--bg:#fbfbfa;--fg:#1d1d1b;--muted:#6b6b66;--card:#fff;--line:#e3e2dd;--good:#1d7a46;--bad:#b3261e;--warn:#9a6700;--accent:#3b5bdb;--okbg:29,122,70;--badbg:179,38,30;--offbg:154,103,0}
@media (prefers-color-scheme:dark){:root{--bg:#161615;--fg:#ecebe6;--muted:#9c9b95;--card:#1f1f1d;--line:#34332f;--good:#5cc98a;--bad:#ff8a80;--warn:#e3b341;--accent:#91a7ff}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}
main{max-width:1180px;margin:0 auto;padding:24px 16px 80px}h1{font-size:26px;margin:0 0 4px}h2{margin-top:40px;border-bottom:1px solid var(--line);padding-bottom:6px}
nav a{margin-right:14px;color:var(--accent);text-decoration:none}.muted{color:var(--muted);font-size:.9em}.ci{color:var(--muted);font-size:.8em}
table{border-collapse:collapse;width:100%;margin:8px 0 16px;font-size:14px}th,td{border-bottom:1px solid var(--line);padding:6px 8px;text-align:left;vertical-align:top}
.headline td,.headline th{text-align:center}.headline td:first-child{text-align:left}
.card{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:12px 14px}.grid2,.grid3{display:grid;gap:14px}
.grid2{grid-template-columns:repeat(auto-fit,minmax(320px,1fr))}.grid3{grid-template-columns:repeat(auto-fit,minmax(260px,1fr))}
table.mini{width:auto}table.mini td,table.mini th{padding:2px 8px;border:0}
table.cm td{text-align:center;width:16%}td.ok{background:rgba(var(--okbg),calc(.08 + .5*var(--a)))}td.flip{background:rgba(var(--badbg),calc(.08*var(--a)*10))}
td.off{background:rgba(var(--offbg),calc(.4*var(--a)))}.good{color:var(--good)}.bad{color:var(--bad)}
tr.fixed td:first-child{color:var(--good);font-weight:600}tr.regressed td:first-child{color:var(--bad);font-weight:600}
details{margin:6px 0}details.variant>summary{font-size:17px;cursor:pointer;padding:6px 0}details.case{background:var(--card);border:1px solid var(--line);border-radius:6px;padding:6px 10px;margin:6px 0 6px 12px}
details.case summary{cursor:pointer}.pill{display:inline-block;min-width:72px;text-align:center;border-radius:10px;padding:0 8px;font-size:12px;font-weight:600;color:#fff}
.pill.ok{background:var(--good)}.pill.off{background:var(--warn)}.pill.flip{background:var(--bad)}pre{white-space:pre-wrap;background:var(--card);border:1px solid var(--line);padding:10px;border-radius:6px;font-size:13px}
code{font-size:.92em}ol.steps{font-size:13px;margin:4px 0}li.bad{color:var(--bad)}table.crit td:nth-child(3){font-size:13px}
@media (max-width:640px){table{font-size:12px}th,td{padding:4px}}
"""


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--flow", default="genomics-interpret")
    p.add_argument("--order", default="baseline,v1-tools,holdout-baseline,holdout-v1")
    p.add_argument("--out", default="", help="output path (default: <flow>/results.html)")
    args = p.parse_args()
    flow = ROOT / ".claude" / "hillclimb" / args.flow
    found = {d.name: Variant(d) for d in sorted(flow.iterdir()) if d.is_dir() and (d / "results.jsonl").exists()}
    order = [n for n in args.order.split(",") if n in found] + [n for n in found if n not in args.order.split(",")]
    variants = [found[n] for n in order]
    cases_by_ds = {ds: [json.loads(l) for l in (DATA / f"{ds}.jsonl").read_text().splitlines()]
                   for ds in dict.fromkeys(v.dataset for v in variants)}
    cases = {c["id"]: c for cs in cases_by_ds.values() for c in cs}
    out = Path(args.out).resolve() if args.out else flow / "results.html"
    out.parent.mkdir(parents=True, exist_ok=True)
    pairs = [(found[a], found[b]) for a, b in (("baseline", "v1-tools"), ("holdout-baseline", "holdout-v1"))
             if a in found and b in found]
    body = "\n".join([
        f"<h1>genomics-mcp · interpret_variant eval</h1><p class=muted>Generated {datetime.now():%Y-%m-%d %H:%M} · "
        f"{sum(len(v.rows) for v in variants)} graded runs across {len(variants)} variants · "
        f"total cost ${sum(v.stats()['cost'] for v in variants):.2f}</p>",
        "<nav><a href=#how>How it works</a><a href=#results>Results</a><a href=#changes>Changes</a><a href=#confusion>Confusion</a>"
        "<a href=#breakdown>Breakdown</a><a href=#tools>Tools</a><a href=#cases>Cases</a></nav>",
        section_overview(variants, cases_by_ds),
        section_headline(variants),
        section_changes(pairs, cases),
        section_confusion(variants),
        section_breakdown(variants, cases),
        section_tools(variants),
        section_cases(variants, cases, out.parent),
    ])
    out.write_text(f"<!doctype html><html lang=en><head><meta charset=utf-8><meta name=viewport content='width=device-width,initial-scale=1'>"
                   f"<title>Genomics MCP Eval</title><style>{CSS}</style></head><body><main>{body}</main></body></html>")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
