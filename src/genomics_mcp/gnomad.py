"""gnomAD GraphQL API: population allele frequencies and gene constraint."""

from typing import Any

from .http import UpstreamError, http
from .variants import Allele

GNOMAD_API = "https://gnomad.broadinstitute.org/api"
DATASETS = {"gnomad_r4", "gnomad_r3", "gnomad_r2_1"}  # r2.1 is GRCh37-only and is rejected below

_SEQ_FIELDS = "ac an homozygote_count hemizygote_count filters faf95 { popmax popmax_population } populations { id ac an homozygote_count }"

VARIANT_QUERY = f"""
query($id: String!, $dataset: DatasetId!) {{
  variant(variantId: $id, dataset: $dataset) {{
    variant_id rsids flags
    exome {{ {_SEQ_FIELDS} }}
    genome {{ {_SEQ_FIELDS} }}
    joint {{ {_SEQ_FIELDS} }}
  }}
}}"""

GENE_QUERY = """
query($symbol: String!) {
  gene(gene_symbol: $symbol, reference_genome: GRCh38) {
    gene_id symbol name chrom start stop strand hgnc_id omim_id
    canonical_transcript_id
    mane_select_transcript { ensembl_id ensembl_version refseq_id refseq_version }
    gnomad_constraint { exp_lof obs_lof oe_lof oe_lof_lower oe_lof_upper pli lof_z exp_mis obs_mis oe_mis oe_mis_upper mis_z syn_z }
  }
}"""

# Genetic-ancestry groups reported in gnomAD v4 (sex-stratified rows are dropped).
ANCESTRY_LABELS = {
    "afr": "African/African American", "amr": "Admixed American", "asj": "Ashkenazi Jewish",
    "eas": "East Asian", "fin": "Finnish", "mid": "Middle Eastern", "nfe": "Non-Finnish European",
    "sas": "South Asian", "ami": "Amish", "remaining": "Remaining", "oth": "Other",
}


class NotInGnomad(Exception):
    pass


async def _query(query: str, variables: dict) -> dict:
    res = await http().post_json(GNOMAD_API, {"query": query, "variables": variables})
    if errs := res.get("errors"):
        msg = "; ".join(e.get("message", "") for e in errs)
        if "not found" in msg.lower():
            raise NotInGnomad(msg)
        raise UpstreamError(f"gnomAD: {msg}")
    return res["data"]


def _summarize(seq: dict | None) -> dict | None:
    if not seq or not seq.get("an"):
        return None
    pops = []
    for p in seq.get("populations") or []:
        if not p["id"] or "_" in p["id"] or p["id"] in {"XX", "XY"} or not p["an"]:
            continue
        pops.append({
            "id": p["id"], "label": ANCESTRY_LABELS.get(p["id"], p["id"]),
            "ac": p["ac"], "an": p["an"], "af": p["ac"] / p["an"], "homozygotes": p["homozygote_count"],
        })
    return {
        "ac": seq["ac"], "an": seq["an"], "af": seq["ac"] / seq["an"],
        "homozygotes": seq["homozygote_count"], "hemizygotes": seq.get("hemizygote_count"),
        "filters": seq.get("filters") or [],
        "grpmax_faf95": (seq.get("faf95") or {}).get("popmax"),
        "grpmax_faf95_group": (seq.get("faf95") or {}).get("popmax_population"),
        "by_ancestry": sorted(pops, key=lambda p: -p["af"]),
    }


async def variant_frequency(allele: Allele, dataset: str = "gnomad_r4") -> dict[str, Any]:
    if dataset not in DATASETS or dataset == "gnomad_r2_1":
        raise ValueError(f"dataset must be gnomad_r4 or gnomad_r3 (GRCh38); got {dataset!r}")
    vid = f"{allele.gnomad_chrom}-{allele.pos}-{allele.ref}-{allele.alt}"
    try:
        v = (await _query(VARIANT_QUERY, {"id": vid, "dataset": dataset}))["variant"]
    except NotInGnomad:
        return {
            "variant": allele.vcf_id, "dataset": dataset, "found": False,
            "note": "Not observed in gnomAD. Absence is evidence of rarity only where the site is well covered.",
        }
    out: dict[str, Any] = {"variant": v["variant_id"], "dataset": dataset, "found": True,
                           "rsids": v.get("rsids") or [], "flags": v.get("flags") or []}
    for key in ("joint", "exome", "genome"):
        out[key] = _summarize(v.get(key))
    return out


async def gene(symbol: str) -> dict[str, Any] | None:
    try:
        return (await _query(GENE_QUERY, {"symbol": symbol.upper()}))["gene"]
    except NotInGnomad:
        return None
