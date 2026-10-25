"""Ensembl VEP (REST) for consequence and in-silico predictions.

Clinical assertions that VEP attaches to co-located variants (ClinVar
significance, HGMD ids, PubMed links) are deliberately dropped: clinical
evidence comes only from the ClinVar tool, so that tool can be withheld in the
blinded eval track without leaking the answer through VEP.
"""

from typing import Any

from .http import http
from .variants import Allele

VEP_URL = "https://rest.ensembl.org/vep/human/region"
VEP_OPTIONS = {
    "canonical": 1, "mane": 1, "hgvs": 1, "numbers": 1,
    "CADD": 1, "REVEL": 1, "AlphaMissense": 1, "SpliceAI": 1, "LoF": 1,
}
MAX_OTHER_TRANSCRIPTS = 15

_PRIMARY_FIELDS = (
    "gene_symbol", "gene_id", "hgnc_id", "transcript_id", "mane_select", "biotype", "strand",
    "consequence_terms", "impact", "hgvsc", "hgvsp", "exon", "intron",
    "protein_start", "protein_end", "amino_acids", "codons",
    "sift_prediction", "sift_score", "polyphen_prediction", "polyphen_score",
    "cadd_phred", "revel", "lof", "lof_filter", "lof_flags", "lof_info",
)


def _spliceai(sp: dict | None) -> dict | None:
    if not sp:
        return None
    scores = {k: sp[k] for k in ("DS_AG", "DS_AL", "DS_DG", "DS_DL") if k in sp}
    return {
        "max_delta_score": max(scores.values(), default=None),
        **scores,
        **{k: sp[k] for k in ("DP_AG", "DP_AL", "DP_DG", "DP_DL") if k in sp},
    }


def _primary(tcs: list[dict]) -> dict | None:
    for pick in (lambda t: t.get("mane_select"), lambda t: t.get("canonical") and t.get("biotype") == "protein_coding",
                 lambda t: t.get("canonical")):
        if hit := next((t for t in tcs if pick(t)), None):
            return hit
    return tcs[0] if tcs else None


async def annotate(allele: Allele) -> dict[str, Any]:
    line = f"{allele.chrom} {allele.pos} . {allele.ref} {allele.alt} . . ."
    [res] = await http().post_json(VEP_URL, {"variants": [line], **VEP_OPTIONS}, headers={"Accept": "application/json"})
    tcs = res.get("transcript_consequences", [])
    primary = _primary(tcs)
    summary: dict[str, Any] = {
        "variant": allele.vcf_id,
        "most_severe_consequence": res.get("most_severe_consequence"),
        "dbsnp_ids": sorted({c["id"] for c in res.get("colocated_variants", []) if c.get("id", "").startswith("rs")}),
    }
    if primary:
        p = {k: primary[k] for k in _PRIMARY_FIELDS if k in primary}
        p["is_mane_select"] = bool(primary.get("mane_select"))
        if am := primary.get("alphamissense"):
            p["alphamissense"] = am
        if sp := _spliceai(primary.get("spliceai")):
            p["spliceai"] = sp
        summary["primary_transcript"] = p
    others = [
        {k: t[k] for k in ("transcript_id", "gene_symbol", "biotype", "consequence_terms", "hgvsc", "hgvsp") if k in t}
        for t in tcs
        if t is not primary and t.get("biotype") == "protein_coding"
    ]
    summary["other_protein_coding_transcripts"] = others[:MAX_OTHER_TRANSCRIPTS]
    if len(others) > MAX_OTHER_TRANSCRIPTS:
        summary["other_transcripts_truncated"] = len(others)
    if not tcs:
        summary["note"] = "No transcript overlaps this variant (intergenic or regulatory only)."
    return summary
