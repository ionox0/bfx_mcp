"""Ensembl VEP (REST) for consequence and in-silico predictions.

Clinical assertions that VEP attaches to co-located variants (ClinVar
significance, HGMD ids, PubMed links) are deliberately dropped: clinical
evidence comes only from the ClinVar tool, so that tool can be withheld in the
blinded eval track without leaking the answer through VEP.
"""

import json
from typing import Any

from . import features
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


_FLANKING = {"upstream_gene_variant", "downstream_gene_variant"}


def _primary(tcs: list[dict], most_severe: str | None, gene: str | None = None) -> dict | None:
    """The transcript to report on: MANE Select, preferring the gene the variant
    actually hits. A variant can sit inside one gene and up/downstream of a
    neighbour's MANE transcript, and VEP lists transcripts in arbitrary order."""
    if gene:
        tcs = [t for t in tcs if (t.get("gene_symbol") or "").upper() == gene.upper()] or tcs
    for pick in (
        lambda t: t.get("mane_select") and most_severe in t.get("consequence_terms", []),
        lambda t: t.get("mane_select") and not _FLANKING.issuperset(t.get("consequence_terms", [])),
        lambda t: t.get("canonical") and t.get("biotype") == "protein_coding" and most_severe in t.get("consequence_terms", []),
        lambda t: t.get("mane_select"),
        lambda t: t.get("canonical") and t.get("biotype") == "protein_coding",
        lambda t: t.get("canonical"),
    ):
        if hit := next((t for t in tcs if pick(t)), None):
            return hit
    return tcs[0] if tcs else None


VEP_BATCH = 50  # REST limit is 200; smaller batches fail less often on a slow server


def _line(allele: Allele) -> str:
    return f"{allele.chrom} {allele.pos} . {allele.ref} {allele.alt} . . ."


def _single_body(allele: Allele) -> dict[str, Any]:
    return {"variants": [_line(allele)], **VEP_OPTIONS}


async def annotate(allele: Allele, gene: str | None = None) -> dict[str, Any]:
    [res] = await http().post_json(VEP_URL, _single_body(allele), headers={"Accept": "application/json"})
    return _summarize(allele, res, gene)


async def coding_position(allele: Allele, gene: str | None = None) -> dict[str, Any] | None:
    """Raw VEP fields for the reported transcript (cds_start, strand, hgvsp...),
    or None if the variant is not in a coding sequence."""
    [res] = await http().post_json(VEP_URL, _single_body(allele), headers={"Accept": "application/json"})
    t = _primary(res.get("transcript_consequences", []), res.get("most_severe_consequence"), gene)
    if not t or t.get("cds_start") is None or t.get("protein_start") is None:
        return None
    return t


async def annotate_many(alleles: list[Allele]) -> dict[str, dict[str, Any]]:
    """Annotate in batched VEP requests, keyed by vcf_id.

    Each result is also cached under the key a single-variant `annotate` call
    would use, so later per-variant tool calls are served from cache.
    """
    out: dict[str, dict[str, Any]] = {}
    for i in range(0, len(alleles), VEP_BATCH):
        chunk = alleles[i : i + VEP_BATCH]
        by_line = {_line(a): a for a in chunk}
        results = await http().post_json(
            VEP_URL, {"variants": list(by_line), **VEP_OPTIONS}, headers={"Accept": "application/json"}
        )
        for res in results:
            allele = by_line.get(res.get("input", ""))
            if allele is None:
                continue
            http().seed("POST", VEP_URL, json.dumps([res]), json_body=_single_body(allele))
            out[allele.vcf_id] = _summarize(allele, res)
    return out


_LOF_TERMS = {"stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant", "start_lost"}


def _pvs1_context(t: dict) -> dict[str, Any] | None:
    """Facts the ClinGen PVS1 decision tree (Abou Tayoun 2018) branches on,
    pulled out of VEP/LOFTEE so the agent doesn't have to decode lof_info."""
    terms = set(t.get("consequence_terms", []))
    if not terms & _LOF_TERMS:
        return None
    info = dict(kv.split(":", 1) for kv in (t.get("lof_info") or "").split(",") if ":" in kv)
    ctx: dict[str, Any] = {
        "lof_consequence": sorted(terms & _LOF_TERMS),
        "loftee": t.get("lof"),  # HC = high-confidence LoF, LC = low-confidence
        "loftee_filter": t.get("lof_filter"),
        "loftee_flags": t.get("lof_flags"),
    }
    if exon := t.get("exon"):
        n, total = (int(x) for x in exon.split("/"))
        ctx["exon"] = exon
        ctx["in_last_exon"] = n == total
    if intron := t.get("intron"):
        ctx["intron"] = intron
    if t.get("cds_start") is not None:
        ctx["cds_position"] = t.get("cds_start")
    if "PERCENTILE" in info:
        ctx["relative_position_in_transcript"] = round(float(info["PERCENTILE"]), 3)
    if "DIST_FROM_LAST_EXON" in info:
        ctx["distance_from_last_exon_junction_nt"] = int(float(info["DIST_FROM_LAST_EXON"]))
    if "50_BP_RULE" in info:
        # FAIL = premature stop within 50 nt of the last exon-exon junction.
        ctx["fifty_nt_rule"] = info["50_BP_RULE"]
    if terms & {"stop_gained", "frameshift_variant"} and "in_last_exon" in ctx:
        ctx["nmd_predicted"] = not ctx["in_last_exon"] and info.get("50_BP_RULE", "PASS") == "PASS"
        ctx["nmd_note"] = (
            "Premature stop is predicted to trigger nonsense-mediated decay."
            if ctx["nmd_predicted"] else
            "Predicted to escape NMD (last exon or within 50 nt of the last junction): PVS1 strength depends on "
            "how much of the protein is removed and whether that region is critical."
        )
    if terms & {"splice_donor_variant", "splice_acceptor_variant"}:
        ctx["splice_note"] = (
            "Canonical +/-1,2 splice site. PVS1 strength depends on the predicted effect (exon skipping in or out "
            "of frame, cryptic site use) and whether the affected exon is in the biologically relevant transcript."
        )
    if "INTRON_SIZE" in info:
        ctx["intron_size_nt"] = int(info["INTRON_SIZE"])
    return ctx


def _summarize(allele: Allele, res: dict[str, Any], gene: str | None = None) -> dict[str, Any]:
    tcs = res.get("transcript_consequences", [])
    primary = _primary(tcs, res.get("most_severe_consequence"), gene)
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
        if features.on("pvs1") and (ctx := _pvs1_context(primary)):
            p["pvs1_context"] = ctx
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
