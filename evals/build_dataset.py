"""Build the three eval tracks from live public data and freeze them.

    uv run python evals/build_dataset.py            # writes evals/data/*.jsonl and evals/cases.md

Every upstream response used to compute gold answers is stored in
evals/data/cache.sqlite with no expiry. The eval runner points the server at
the same cache, so the agent sees exactly the data the gold was computed from
even though ClinVar and gnomAD keep changing.

Tracks
  reasoning  ClinVar hidden. Classify variants whose gold label is a ClinVar
             expert-panel / practice-guideline classification.
  tool_use   Factual questions with deterministic answers from gnomAD / VEP /
             normalization (plus a few error traps).
  retrieval  ClinVar visible. Report the current aggregate classification,
             review stars and conflict status.
"""

import os
import sys
from pathlib import Path

DATA = Path(__file__).parent / "data"
DATA.mkdir(exist_ok=True)
os.environ["GENOMICS_MCP_CACHE"] = str(DATA / "cache.sqlite")
os.environ["GENOMICS_MCP_CACHE_TTL_S"] = "inf"

import asyncio  # noqa: E402
import json  # noqa: E402
import random  # noqa: E402
import re  # noqa: E402
from datetime import date  # noqa: E402

from genomics_mcp import clinvar, ensembl, gnomad  # noqa: E402
from genomics_mcp.http import UpstreamError, log  # noqa: E402
from genomics_mcp.variants import Allele, VariantInputError, check_reference, parse_vcf, resolve  # noqa: E402

SEED = 20261025
TIERS = ["Pathogenic", "Likely pathogenic", "Uncertain significance", "Likely benign", "Benign"]
TIER_FILTER = dict(zip(TIERS, ["clinsig_pathogenic", "clinsig_likely_pathogenic", "clinsig_vus",
                               "clinsig_likely_benign", "clinsig_benign"]))
REASONING_GENES = ["BRCA1", "BRCA2", "TP53", "PTEN", "CDH1", "MLH1", "MSH2", "MYH7", "PAH", "PTPN11",
                   "RAF1", "SOS1", "RYR1", "FBN1", "GAA", "CDKL5", "ATM", "DICER1", "RUNX1", "HNF1A", "GCK", "LDLR"]
PER_TIER = 12
MAX_PER_GENE_TIER = 2
NAME_RE = re.compile(r"^(NM_\d+\.\d+)\(([A-Za-z0-9-]+)\):(c\.[^\s]+)")
SKIP_CONDITIONS = {"not provided", "not specified", "See cases"}


async def vep(allele: Allele) -> dict | None:
    """VEP with extra patience for Ensembl outages; None if it stays down."""
    for round_ in range(3):
        try:
            return await ensembl.annotate(allele)
        except UpstreamError as e:
            log(f"VEP round {round_ + 1}/3 failed for {allele.vcf_id}: {str(e)[:80]}")
            await asyncio.sleep(30)
    return None


def hgvs_of(name: str) -> tuple[str, str] | None:
    m = NAME_RE.match(name or "")
    return (f"{m[1]}:{m[3]}", m[2]) if m else None


def condition_of(rec: dict) -> str | None:
    return next((c for c in rec["conditions"] if c not in SKIP_CONDITIONS), None)


def usable(va, rec: dict, allele: Allele | None) -> bool:
    return (
        allele is not None
        and allele.ref and allele.alt
        and max(len(allele.ref), len(allele.alt)) <= 50
        and allele.chrom != "MT"
        and hgvs_of(rec["name"]) is not None
        and condition_of(rec) is not None
    )


# ----------------------------------------------------------------------------- reasoning


async def build_reasoning(rng: random.Random) -> list[dict]:
    cases: list[dict] = []
    for tier in TIERS:
        picked: list[dict] = []
        genes = REASONING_GENES[:]
        rng.shuffle(genes)
        # Round-robin across genes so no single VCEP dominates a tier.
        for round_ in range(MAX_PER_GENE_TIER):
            for gene in genes:
                if len(picked) >= PER_TIER:
                    break
                status = '("reviewed by expert panel"[Review status] OR "practice guideline"[Review status])'
                ids = await clinvar.search_ids(f"{gene}[gene] AND {status} AND {TIER_FILTER[tier]}[Properties]", 300)
                rng_g = random.Random(f"{SEED}-{gene}-{tier}")
                rng_g.shuffle(ids)
                for va, rec in await clinvar.fetch_by_variation_ids(ids[:12]):
                    allele = clinvar.vcf_allele_of(va)
                    if (
                        rec["germline_classification"] == tier
                        and rec["review_stars"] >= 3
                        and usable(va, rec, allele)
                        and not any(c["variation_id"] == rec["variation_id"] for c in picked)
                        and sum(c["gene"] == gene for c in picked) <= round_
                    ):
                        hgvs, sym = hgvs_of(rec["name"])
                        picked.append({
                            "variation_id": rec["variation_id"], "gene": sym, "hgvs": hgvs, "allele": allele,
                            "condition": condition_of(rec), "rec": rec,
                        })
                        break
        log(f"reasoning: {tier}: {len(picked)} cases")
        for p in picked:
            ann = await vep(p["allele"])
            consequence = (
                (ann.get("primary_transcript") or {}).get("consequence_terms", [ann["most_severe_consequence"]])
                if ann else ["unannotated"]
            )
            cases.append({
                "id": f"reasoning-{p['gene']}-{p['variation_id']}",
                "track": "reasoning",
                "tags": [tier, p["gene"], consequence[0]],
                "input": {"variant": p["hgvs"], "gene": p["gene"], "condition": p["condition"]},
                "expected": {"classification": tier},
                "meta": {
                    "clinvar_accession": p["rec"]["accession"], "review_status": p["rec"]["review_status"],
                    "date_last_evaluated": p["rec"]["date_last_evaluated"], "vcf_id": p["allele"].vcf_id,
                },
            })
    return cases


# ----------------------------------------------------------------------------- tool use

TOOL_USE_VARIANTS = [
    # (input as given to the agent, why it is here)
    ("NM_000492.4:c.1521_1523del", "common recessive pathogenic deletion"),
    ("NM_007294.4:c.5266dup", "frameshift, LOFTEE HC"),
    ("17-43045712-T-C", "missense with REVEL/AlphaMissense"),
    ("NM_000059.4:c.7007G>A", "BRCA2 last-base-of-exon missense with splice impact"),
    ("NM_000546.6:c.743G>A", "TP53 hotspot missense"),
    ("NM_000277.3:c.1222C>T", "PAH recessive, observed homozygotes"),
    ("NM_000257.4:c.1208G>A", "MYH7 missense"),
    ("NM_004333.6:c.1799T>A", "BRAF V600E (somatic hotspot, rare germline)"),
    ("NM_000038.6:c.3927_3931del", "APC frameshift (p.Glu1309fs)"),
    ("NM_000051.4:c.7271T>G", "ATM missense"),
]


def q(id_, question, answer_type, gold, **grading):
    return {"id": f"tool_use-{id_}", "track": "tool_use", "tags": [grading.pop("kind")],
            "input": {"question": question}, "expected": {"answer": gold, "type": answer_type, **grading}}


async def build_tool_use() -> list[dict]:
    cases: list[dict] = []
    for i, (variant, _why) in enumerate(TOOL_USE_VARIANTS):
        [allele] = await resolve(variant)
        freq = await gnomad.variant_frequency(allele)
        ann = await vep(allele)
        if ann is None:
            log(f"tool_use: skipping annotation questions for {variant} (Ensembl down)")
        pt = (ann or {}).get("primary_transcript") or {}
        tag = f"v{i}"
        if variant != allele.vcf_id:
            cases.append(q(f"{tag}-normalize", f"What is the GRCh38 VCF-style id (chrom-pos-ref-alt) of {variant}?",
                           "string", allele.vcf_id, kind="normalization"))
        joint = freq.get("joint") if freq["found"] else None
        if joint:
            cases.append(q(f"{tag}-af", f"What is the gnomAD v4 joint (exome+genome) allele frequency of {variant}?",
                           "number", joint["af"], rel_tol=0.02, kind="frequency"))
            cases.append(q(f"{tag}-hom", f"How many homozygous individuals carry {variant} in the gnomAD v4 joint callset?",
                           "integer", joint["homozygotes"], kind="frequency"))
            if joint["grpmax_faf95_group"]:
                cases.append(q(f"{tag}-grpmax", f"In gnomAD v4 (joint), which genetic-ancestry group has the highest filtering allele frequency (grpmax FAF95) for {variant}? Answer with the gnomAD group code.",
                               "string", joint["grpmax_faf95_group"], kind="frequency"))
        else:
            cases.append(q(f"{tag}-absent", f"Is {variant} observed in gnomAD v4? Answer yes or no.",
                           "string", "no", kind="frequency"))
        if pt:
            cases.append(q(f"{tag}-csq", f"What is the Sequence Ontology consequence term of {variant} on its MANE Select transcript, per Ensembl VEP?",
                           "string", pt["consequence_terms"][0], kind="annotation"))
            if pt.get("hgvsp"):
                cases.append(q(f"{tag}-hgvsp", f"What is the HGVS protein change of {variant} on the MANE Select transcript (p. notation, three-letter amino acids)?",
                               "hgvsp", pt["hgvsp"].split(":", 1)[1], kind="annotation"))
            if pt.get("revel") is not None:
                cases.append(q(f"{tag}-revel", f"What is the REVEL score of {variant}?",
                               "number", pt["revel"], abs_tol=0.0015, kind="prediction"))
            if (sp := pt.get("spliceai")) and sp.get("max_delta_score") is not None:
                cases.append(q(f"{tag}-spliceai", f"What is the maximum SpliceAI delta score for {variant}?",
                               "number", sp["max_delta_score"], abs_tol=0.0015, kind="prediction"))
    for gene in ["BRCA1", "TP53", "PAH", "SCN1A", "MYH7"]:
        g = await gnomad.gene(gene)
        cases.append(q(f"loeuf-{gene}", f"What is the gnomAD v4 LOEUF (upper bound of the LoF observed/expected CI) for {gene}?",
                       "number", g["gnomad_constraint"]["oe_lof_upper"], abs_tol=0.005, kind="constraint"))
    # Multi-allelic rsID: the agent must report every alt allele, not the first one.
    alleles = await resolve("rs80357258")
    cases.append(q("rs-multiallelic", "List every GRCh38 VCF-style allele id (chrom-pos-ref-alt) for rs80357258.",
                   "string_set", sorted(a.vcf_id for a in alleles), kind="normalization"))
    # Traps: GRCh37 coordinates whose REF does not match GRCh38 must be reported, not answered.
    for grch37 in ["17-41197729-T-C", "7-117199644-ATCT-A", "13-32914438-T-TA"]:
        try:
            await check_reference(parse_vcf(grch37))
        except VariantInputError:
            cases.append(q(f"trap-{grch37}", f"What is the gnomAD v4 allele frequency of {grch37}? (coordinates as given)",
                           "invalid", "invalid", kind="trap"))
        else:
            log(f"tool_use: {grch37} happens to match GRCh38; not usable as a trap")
    return cases


# ----------------------------------------------------------------------------- retrieval

RETRIEVAL_STATUSES = {
    "criteria provided, conflicting classifications": 8,
    "criteria provided, single submitter": 6,
    "criteria provided, multiple submitters, no conflicts": 6,
    "reviewed by expert panel": 4,
}
RETRIEVAL_GENES = ["BRCA2", "MSH6", "LDLR", "SCN5A", "KCNH2", "MYBPC3", "TTN", "APC", "CHEK2", "ATM", "RYR1", "DSP"]


async def build_retrieval(rng: random.Random) -> list[dict]:
    cases: list[dict] = []
    for status, n in RETRIEVAL_STATUSES.items():
        genes = RETRIEVAL_GENES[:]
        rng.shuffle(genes)
        picked = 0
        for gene in genes:
            if picked >= n:
                break
            ids = await clinvar.search_ids(f'{gene}[gene] AND "{status}"[Review status]', 300)
            random.Random(f"{SEED}-{gene}-{status}").shuffle(ids)
            for va, rec in await clinvar.fetch_by_variation_ids(ids[:10]):
                allele = clinvar.vcf_allele_of(va)
                if rec["review_status"] != status or not usable(va, rec, allele):
                    continue
                # Re-run the server's own lookup so gold matches what the tool returns for this allele.
                live = await clinvar.lookup(allele)
                if len(live) != 1:
                    continue
                r = live[0]
                hgvs, sym = hgvs_of(r["name"])
                cases.append({
                    "id": f"retrieval-{sym}-{r['variation_id']}",
                    "track": "retrieval",
                    "tags": [status, sym],
                    "input": {"variant": hgvs},
                    "expected": {
                        "classification": r["germline_classification"],
                        "review_stars": r["review_stars"],
                        "has_conflict": "conflicting" in (r["review_status"] or ""),
                        "accession": r["accession"].split(".")[0],
                    },
                    "meta": {"vcf_id": allele.vcf_id, "breakdown": r["classification_breakdown"]},
                })
                picked += 1
                break
        log(f"retrieval: {status}: {picked}")
    return cases


# ----------------------------------------------------------------------------- review file


def write_markdown(tracks: dict[str, list[dict]]) -> None:
    lines = [f"# Eval cases (built {date.today().isoformat()}, seed {SEED})", ""]
    for track, cases in tracks.items():
        lines += [f"## {track} ({len(cases)} cases)", ""]
        if track == "reasoning":
            lines += ["| id | gold | gene | variant | condition | consequence | ClinVar | evaluated |", "|---|---|---|---|---|---|---|---|"]
            for c in cases:
                lines.append(f"| {c['id']} | {c['expected']['classification']} | {c['input']['gene']} | `{c['input']['variant']}` | "
                             f"{c['input']['condition']} | {c['tags'][2]} | {c['meta']['clinvar_accession']} | {c['meta']['date_last_evaluated']} |")
        elif track == "tool_use":
            lines += ["| id | kind | question | gold |", "|---|---|---|---|"]
            for c in cases:
                lines.append(f"| {c['id']} | {c['tags'][0]} | {c['input']['question']} | `{c['expected']['answer']}` |")
        else:
            lines += ["| id | variant | classification | stars | conflict | breakdown |", "|---|---|---|---|---|---|"]
            for c in cases:
                e = c["expected"]
                lines.append(f"| {c['id']} | `{c['input']['variant']}` | {e['classification']} | {e['review_stars']} | "
                             f"{e['has_conflict']} | {c['meta']['breakdown']} |")
        lines.append("")
    (Path(__file__).parent / "cases.md").write_text("\n".join(lines))


async def main(only: set[str]) -> None:
    rng = random.Random(SEED)
    builders = {"reasoning": build_reasoning, "tool_use": build_tool_use, "retrieval": build_retrieval}
    tracks: dict[str, list[dict]] = {}
    for name, build in builders.items():
        path = DATA / f"{name}.jsonl"
        if only and name not in only and path.exists():
            tracks[name] = [json.loads(line) for line in path.read_text().splitlines()]
            continue
        tracks[name] = await (build(rng) if name != "tool_use" else build())
        path.write_text("".join(json.dumps(c) + "\n" for c in tracks[name]))
        log(f"wrote {path} ({len(tracks[name])} cases)")
    write_markdown(tracks)


if __name__ == "__main__":
    asyncio.run(main(set(sys.argv[1:])))
