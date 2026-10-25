"""Programmatic graders. Every track's answer is a constrained value (a 5-tier
label, a number, an id, a boolean), so no model judge is needed for the
headline metrics. Each grader returns {metric_id: 0.0-1.0}.
"""

import math
import re
from typing import Any

TIERS = ["Pathogenic", "Likely pathogenic", "Uncertain significance", "Likely benign", "Benign"]
THREE_TIER = {
    "Pathogenic": "P/LP", "Likely pathogenic": "P/LP",
    "Uncertain significance": "VUS",
    "Likely benign": "B/LB", "Benign": "B/LB",
}

_CLASSIFY_METRICS = [
    {"id": "three_tier", "label": "3-tier match", "kind": "binary"},
    {"id": "exact", "label": "5-tier exact", "kind": "binary"},
    {"id": "severe_error", "label": "P<->B flip", "kind": "binary", "better": "lower"},
    {"id": "within_one", "label": "within 1 tier", "kind": "binary"},
    {"id": "self_consistent", "label": "self-consistent", "kind": "binary"},
]

METRICS = {
    "reasoning": _CLASSIFY_METRICS,
    "interpret": _CLASSIFY_METRICS,
    "tool_use": [
        {"id": "correct", "label": "correct", "kind": "binary"},
    ],
    "retrieval": [
        {"id": "all_fields", "label": "all fields", "kind": "binary"},
        {"id": "classification", "label": "class", "kind": "binary"},
        {"id": "stars", "label": "stars", "kind": "binary"},
        {"id": "conflict", "label": "conflict", "kind": "binary"},
        {"id": "accession", "label": "accession", "kind": "binary"},
    ],
}


def _norm(s: Any) -> str:
    return re.sub(r"[\s_\-]+", " ", str(s)).strip().lower()


_POINTS = {"very strong": 8, "strong": 4, "moderate": 2, "supporting": 1}


def tier_from_criteria(criteria: list[dict]) -> str | None:
    """5-tier class implied by the stated criteria under the ACMG point system
    (Tavtigian 2018/2020). None if the criteria can't be scored."""
    total = 0
    for c in criteria:
        code = str(c.get("code", "")).strip().upper()
        strength = str(c.get("strength", "")).strip().lower()
        if code.startswith("BA1") or strength == "stand-alone":
            if code.startswith("B"):
                return "Benign"
            continue
        if code[:1] not in {"P", "B"} or strength not in _POINTS:
            return None
        total += _POINTS[strength] if code.startswith("P") else -_POINTS[strength]
    if total >= 10:
        return "Pathogenic"
    if total >= 6:
        return "Likely pathogenic"
    if total >= 0:
        return "Uncertain significance"
    return "Likely benign" if total >= -6 else "Benign"


def grade_reasoning(expected: dict, answer: dict) -> dict[str, float]:
    gold, got = expected["classification"], answer.get("classification")
    if got not in TIERS:
        return {"three_tier": 0.0, "exact": 0.0, "severe_error": 0.0, "within_one": 0.0, "self_consistent": 0.0}
    g3, a3 = THREE_TIER[gold], THREE_TIER[got]
    return {
        "three_tier": float(g3 == a3),
        "exact": float(gold == got),
        "severe_error": float({g3, a3} == {"P/LP", "B/LB"}),
        "within_one": float(abs(TIERS.index(gold) - TIERS.index(got)) <= 1),
        # Reference-free: does the submitted tier follow from the agent's own criteria?
        "self_consistent": float(tier_from_criteria(answer.get("criteria") or []) == got),
    }


def _number(s: Any) -> float | None:
    if isinstance(s, (int, float)):
        return float(s)
    m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", str(s).replace(",", ""))
    return float(m[0]) if m else None


def _hgvsp(s: str) -> str:
    s = str(s).strip()
    s = s.split(":", 1)[-1]  # drop ENSP/NP accession
    s = s.removeprefix("p.").strip("()")
    return s.replace("fsTer", "fs*").replace("Ter", "*").lower()


def _vcf_id(s: str) -> str:
    return _norm(s).replace("chr", "").replace(":", " ").replace(" ", "-")


def grade_tool_use(expected: dict, answer: dict) -> dict[str, float]:
    kind, gold = expected["type"], expected["answer"]
    got = answer.get("answer", "")
    if kind == "invalid":
        return {"correct": float(bool(answer.get("input_invalid")))}
    if answer.get("input_invalid"):
        return {"correct": 0.0}
    if kind == "number":
        x = _number(got)
        ok = x is not None and math.isclose(
            x, gold, rel_tol=expected.get("rel_tol", 0.0), abs_tol=expected.get("abs_tol", 0.0) or 1e-12
        )
        return {"correct": float(ok)}
    if kind == "integer":
        x = _number(got)
        return {"correct": float(x is not None and int(x) == int(gold))}
    if kind == "hgvsp":
        return {"correct": float(_hgvsp(got) == _hgvsp(gold))}
    if kind == "string_set":
        items = {_vcf_id(x) for x in re.split(r"[,;\s]+(?=(?:chr)?[0-9XYM]+[-:])", str(got)) if x.strip()}
        return {"correct": float(items == {_vcf_id(x) for x in gold})}
    if re.fullmatch(r"(chr)?[0-9XYM]+-\d+-[ACGT]+-[ACGT]+", str(gold)):
        return {"correct": float(_vcf_id(got) == _vcf_id(gold))}
    return {"correct": float(_norm(got) == _norm(gold))}


def grade_retrieval(expected: dict, answer: dict) -> dict[str, float]:
    scores = {
        "classification": float(_norm(answer.get("classification", "")) == _norm(expected["classification"])),
        "stars": float(answer.get("review_stars") == expected["review_stars"]),
        "conflict": float(answer.get("has_conflicting_submissions") == expected["has_conflict"]),
        "accession": float(_norm(answer.get("clinvar_accession", "")).split(".")[0] == _norm(expected["accession"])),
    }
    return {"all_fields": float(all(scores.values())), **scores}


GRADERS = {"reasoning": grade_reasoning, "interpret": grade_reasoning, "tool_use": grade_tool_use,
           "retrieval": grade_retrieval}
