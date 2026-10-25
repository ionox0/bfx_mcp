"""genomics-mcp: public variant-interpretation sources for AI agents, over MCP.

Tools (model-invoked)
  normalize_variant           any notation -> GRCh38 left-aligned VCF allele(s)
  annotate_variant            Ensembl VEP: consequence, HGVS, REVEL/AlphaMissense/CADD/SpliceAI, LOFTEE
  get_population_frequency    gnomAD v4/v3: AC/AN/AF, homozygotes, grpmax FAF95, per-ancestry AF
  get_gene_constraint         gnomAD: pLI, LOEUF, missense Z, MANE transcript
  get_clinvar_classification  ClinVar: aggregate classification, review status, submissions
  get_clinvar_codon_variants  ClinVar: other alleles at the same codon (PS1/PM5), own record excluded

Resources (app/user-attached)
  genomics://sources                 provenance, versions and limitations
  genomics://acmg-criteria           ACMG/AMP criteria summary mapped to these tools
  genomics://variant/{variant}       one-shot evidence bundle for a variant

Prompts (user-picked)
  interpret_variant                  ACMG-style interpretation workflow

Set GENOMICS_MCP_HIDE_CLINVAR=1 to withhold the variant's own ClinVar record
(used by the blinded eval tracks, where ClinVar classifications are the gold
labels). GENOMICS_MCP_FEATURES selects server improvements (see features.py).
"""

import asyncio
import json
import math
import os
from importlib.resources import files
from typing import Any, Literal

from mcp.server.mcpserver import MCPServer
from mcp.server.mcpserver.exceptions import ToolError
from mcp.types import ToolAnnotations

from . import clinvar, ensembl, features, gnomad
from .http import CACHE_TTL_S, UpstreamError, log
from .variants import Allele, VariantInputError, classify_input, resolve, resolve_one

HIDE_CLINVAR = os.environ.get("GENOMICS_MCP_HIDE_CLINVAR", "").lower() in {"1", "true", "yes"}

INSTRUCTIONS = """\
Public genomics evidence for germline variant interpretation (GRCh38 only).

Variants can be given as VCF-style ids (17-43045712-T-C), rsIDs, or HGVS with a
versioned RefSeq accession (NM_007294.4:c.5266dup). Every tool normalizes its
input; call normalize_variant first only when an rsID may be multi-allelic or
you need the canonical id. Upstream APIs are rate-limited and Ensembl can be
slow (10-60 s); call tools for one variant at a time and reuse results.

Read genomics://acmg-criteria for how each tool's output maps to ACMG/AMP
criteria. State which criteria you apply, at what strength, and which ones
cannot be evaluated from public aggregate data.

For research and education only. Output is not a clinical interpretation.
"""

if HIDE_CLINVAR:
    INSTRUCTIONS += "\nClinVar is not available in this session: interpret from primary evidence.\n"

mcp = MCPServer(name="genomics", instructions=INSTRUCTIONS, version="0.1.0")

READ_ONLY = ToolAnnotations(read_only_hint=True, idempotent_hint=True, open_world_hint=True)


async def _allele(variant: str) -> Allele:
    try:
        return await resolve_one(variant)
    except VariantInputError as e:
        raise ToolError(str(e)) from e
    except UpstreamError as e:
        raise ToolError(f"Could not normalize {variant!r}: {e}") from e


def _upstream(source: str, e: Exception) -> ToolError:
    log(f"{source} failed: {e}")
    return ToolError(f"{source} is unavailable right now ({e}). Retry later or continue without it.")


# --------------------------------------------------------------------------- tools


@mcp.tool(annotations=READ_ONLY)
async def normalize_variant(variant: str) -> dict[str, Any]:
    """Resolve a variant to GRCh38, left-aligned VCF-style allele(s).

    Accepts VCF-style ids (17-43045712-T-C, chr17:43045712 T>C), rsIDs
    (rs80357258) and HGVS with versioned RefSeq accessions (NM_007294.4:c.5266dup,
    NC_000017.11:g.43045712T>C). A multi-allelic rsID returns several alleles;
    pass one of their `vcf_id`s to the other tools.
    """
    try:
        alleles = await resolve(variant)
    except VariantInputError as e:
        raise ToolError(str(e)) from e
    except UpstreamError as e:
        raise _upstream("NCBI Variation Services", e) from e
    return {"input": variant, "input_type": classify_input(variant), "alleles": [a.as_dict() for a in alleles]}


@mcp.tool(annotations=READ_ONLY)
async def annotate_variant(variant: str, gene: str | None = None) -> dict[str, Any]:
    """Predicted molecular consequence and in-silico scores from Ensembl VEP.

    Returns the most severe consequence; for the MANE Select transcript the
    variant falls in (else canonical): HGVS c./p., exon/intron, impact, SIFT,
    PolyPhen, CADD, REVEL, AlphaMissense, SpliceAI delta scores and LOFTEE
    loss-of-function calls; plus consequences on other protein-coding
    transcripts. Pass `gene` (HGNC symbol) to report on that gene's transcript
    where genes overlap. Contains no clinical assertions. Ensembl can take 10-60 s.
    """
    allele = await _allele(variant)
    try:
        return await ensembl.annotate(allele, gene)
    except UpstreamError as e:
        raise _upstream("Ensembl VEP", e) from e


@mcp.tool(annotations=READ_ONLY)
async def get_population_frequency(
    variant: str, dataset: Literal["gnomad_r4", "gnomad_r3"] = "gnomad_r4"
) -> dict[str, Any]:
    """Population allele frequency from gnomAD (default v4, GRCh38).

    Returns allele count/number/frequency and homozygote counts for the joint
    (exome+genome), exome and genome callsets, the grpmax filtering allele
    frequency (FAF95, the value ACMG BA1/BS1 thresholds should be compared to),
    QC filters, and per-genetic-ancestry frequencies. `found: false` means the
    allele was not observed.
    """
    allele = await _allele(variant)
    try:
        return await gnomad.variant_frequency(allele, dataset)
    except UpstreamError as e:
        raise _upstream("gnomAD", e) from e


@mcp.tool(annotations=READ_ONLY)
async def get_gene_constraint(gene_symbol: str) -> dict[str, Any]:
    """Gene-level constraint and transcript info from gnomAD v4.

    Returns pLI, LOEUF (oe_lof_upper; < 0.6 suggests loss-of-function
    intolerance), observed/expected LoF and missense counts, missense and
    synonymous Z scores, coordinates, and the MANE Select transcript. Constraint
    is context for disease mechanism, not proof of it.
    """
    try:
        g = await gnomad.gene(gene_symbol.strip())
    except UpstreamError as e:
        raise _upstream("gnomAD", e) from e
    if g is None:
        raise ToolError(f"gnomAD has no GRCh38 gene with symbol {gene_symbol!r}. Use the HGNC-approved symbol.")
    return g


async def get_clinvar_classification(variant: str) -> dict[str, Any]:
    """ClinVar germline classification for an exact allele.

    Returns each matching ClinVar variation record (VCV) with the aggregate
    germline classification, review status and star rating (0-4), date last
    evaluated, conditions, and per-submitter classifications (first 25). Only
    records whose GRCh38 VCF allele matches exactly are returned; an empty list
    means ClinVar has no record of this allele.
    """
    allele = await _allele(variant)
    try:
        records = await clinvar.lookup(allele)
    except UpstreamError as e:
        raise _upstream("ClinVar", e) from e
    return {"variant": allele.vcf_id, "records": records}


if not HIDE_CLINVAR:
    mcp.tool(annotations=READ_ONLY)(get_clinvar_classification)


async def get_clinvar_codon_variants(variant: str, gene: str | None = None) -> dict[str, Any]:
    """Other ClinVar variants at the same codon, for ACMG PS1 and PM5.

    Returns ClinVar records (classification, review stars) for other alleles in
    the codon this coding variant changes, each labelled same_amino_acid_change
    (PS1 candidate), different_missense (PM5 candidate), nonsense, other_protein_change or synonymous.
    The queried variant's own ClinVar record is always excluded, so this tool is
    available even when get_clinvar_classification is not. Pass `gene` where
    genes overlap.
    """
    allele = await _allele(variant)
    try:
        transcript = await ensembl.coding_position(allele, gene)
        if transcript is None:
            raise ToolError(f"{allele.vcf_id} is not in a coding sequence, so PS1/PM5 do not apply.")
        return await clinvar.codon_variants(allele, transcript)
    except UpstreamError as e:
        raise _upstream("Ensembl VEP / ClinVar", e) from e


if features.on("codon"):
    mcp.tool(annotations=READ_ONLY)(get_clinvar_codon_variants)


# --------------------------------------------------------------------------- resources

SOURCES_MD = f"""\
# Data sources

| Source | What | Access | Notes |
|---|---|---|---|
| NCBI Variation Services | HGVS / rsID / SPDI -> VCF normalization | api.ncbi.nlm.nih.gov/variation/v0 | GRCh38 primary assembly only |
| Ensembl VEP | Consequences, HGVS, SIFT, PolyPhen, CADD, REVEL, AlphaMissense, SpliceAI, LOFTEE | rest.ensembl.org | Current Ensembl release; clinical co-located annotations stripped |
| gnomAD | Allele frequencies (v4.1 default, v3.1), gene constraint | gnomad.broadinstitute.org/api | Strict rate limit; ~1 request/s |
| ClinVar | Aggregate and per-submitter germline classifications | NCBI E-utilities | {"**Withheld in this session** (GENOMICS_MCP_HIDE_CLINVAR)" if HIDE_CLINVAR else "Live; reclassifications appear as ClinVar updates"} |

Responses are cached on disk ({"no expiry: frozen snapshot" if math.isinf(CACHE_TTL_S) else f"{CACHE_TTL_S / 3600:.0f} h"}).
Set NCBI_API_KEY to raise the NCBI rate limit from 3 to 10 requests/s.

Not covered: structural/copy-number variants, somatic interpretation,
mitochondrial heteroplasmy, GRCh37 coordinates, literature and functional
data, case-level segregation or de novo data.
"""


@mcp.resource("genomics://sources", mime_type="text/markdown")
def sources() -> str:
    """Provenance, versions and limitations of every upstream source."""
    return SOURCES_MD


@mcp.resource("genomics://acmg-criteria", mime_type="text/markdown")
def acmg_criteria() -> str:
    """ACMG/AMP criteria summary, combining rules, and which tool supplies each criterion's evidence."""
    text = files("genomics_mcp").joinpath("acmg_criteria.md").read_text()
    if features.on("codon"):
        text = text.replace("ClinVar (other alleles at the codon)", "`get_clinvar_codon_variants`")
    if features.on("pvs1"):
        text = text.replace("`annotate_variant` (consequence, LOFTEE", "`annotate_variant` (`pvs1_context`: NMD prediction, last exon; consequence, LOFTEE")
    return text


@mcp.resource("genomics://variant/{variant}", mime_type="application/json")
async def variant_bundle(variant: str) -> str:
    """Every available source for one variant in a single document (attach with @)."""
    allele = await _allele(variant)
    jobs: dict[str, Any] = {
        "consequence": ensembl.annotate(allele),
        "population_frequency": gnomad.variant_frequency(allele),
    }
    if not HIDE_CLINVAR:
        jobs["clinvar"] = clinvar.lookup(allele)
    results = await asyncio.gather(*jobs.values(), return_exceptions=True)
    bundle: dict[str, Any] = {"variant": allele.as_dict()}
    for key, res in zip(jobs, results):
        bundle[key] = {"error": str(res)} if isinstance(res, Exception) else res
    return json.dumps(bundle, indent=2)


# --------------------------------------------------------------------------- prompts


@mcp.prompt()
def interpret_variant(variant: str, condition: str = "", inheritance: str = "") -> str:
    """Walk through an ACMG/AMP-style interpretation of one germline variant."""
    context = "".join([
        f"\nCondition of interest: {condition}" if condition else "",
        f"\nMode of inheritance: {inheritance}" if inheritance else "",
    ])
    clinvar_step = (
        "" if HIDE_CLINVAR else
        "5. Check ClinVar last, and compare: does the existing classification agree with your evidence? "
        "Note review status and any conflicts.\n"
    )
    return f"""Interpret the germline variant {variant}.{context}

1. Normalize the variant and confirm the gene and MANE Select transcript.
2. Get the predicted consequence and in-silico scores.
3. Get gnomAD frequency (use grpmax FAF95 for frequency criteria) and gene constraint.
4. Read genomics://acmg-criteria. For each criterion you can evaluate, state the code, strength and the specific value that triggers it. List the criteria that need data not available here.
{clinvar_step}
Finish with a 5-tier classification (Pathogenic, Likely pathogenic, Uncertain significance, Likely benign, Benign), the criteria combination that supports it, and your main uncertainty. This is a research exercise, not a clinical report."""


def main() -> None:
    log(f"starting on stdio (ClinVar {'hidden' if HIDE_CLINVAR else 'enabled'})")
    mcp.run("stdio")


if __name__ == "__main__":
    main()
