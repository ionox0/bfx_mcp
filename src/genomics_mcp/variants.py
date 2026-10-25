"""Variant identifiers: parse what the agent typed and normalize it to one
GRCh38, left-aligned, VCF-style allele that every upstream understands.

Accepted inputs:
  * VCF-style      17-43045712-T-C, chr17:43045712:T:C, 17:43045712 T>C
  * dbSNP          rs80357258        (may expand to several alt alleles)
  * HGVS           NM_007294.4:c.5266dup, NC_000017.11:g.43045712T>C

HGVS and rsIDs are resolved with NCBI Variation Services, which return SPDI
and convert it to left-aligned VCF fields (the representation gnomAD uses).
VCF-style input is checked against the GRCh38 reference via the UCSC API.
"""

import re
from dataclasses import asdict, dataclass
from urllib.parse import quote

from .http import UpstreamError, http

NCBI_VARIATION = "https://api.ncbi.nlm.nih.gov/variation/v0"
UCSC_SEQUENCE = "https://api.genome.ucsc.edu/getData/sequence"

# GRCh38 RefSeq accessions for the primary assembly.
GRCH38_ACCESSIONS = {
    "1": "NC_000001.11", "2": "NC_000002.12", "3": "NC_000003.12", "4": "NC_000004.12",
    "5": "NC_000005.10", "6": "NC_000006.12", "7": "NC_000007.14", "8": "NC_000008.11",
    "9": "NC_000009.12", "10": "NC_000010.11", "11": "NC_000011.10", "12": "NC_000012.12",
    "13": "NC_000013.11", "14": "NC_000014.9", "15": "NC_000015.10", "16": "NC_000016.10",
    "17": "NC_000017.11", "18": "NC_000018.10", "19": "NC_000019.10", "20": "NC_000020.11",
    "21": "NC_000021.9", "22": "NC_000022.11", "X": "NC_000023.11", "Y": "NC_000024.10",
    "MT": "NC_012920.1",
}
_ACCESSION_TO_CHROM = {acc: chrom for chrom, acc in GRCH38_ACCESSIONS.items()}

_VCF_RE = re.compile(
    r"^(?:chr)?(?P<chrom>[0-9]{1,2}|X|Y|MT?)[-:_\s]+(?P<pos>\d+)[-:_\s]+(?P<ref>[ACGTN]+)(?:[-:_\s]+|\s*>\s*)(?P<alt>[ACGTN]+)$",
    re.IGNORECASE,
)
_RSID_RE = re.compile(r"^rs\d+$", re.IGNORECASE)
_HGVS_RE = re.compile(r"^[A-Z]{2}_\d+(\.\d+)?(\([A-Za-z0-9-]+\))?:[cgnmr]\.", re.IGNORECASE)


@dataclass(frozen=True)
class Allele:
    """One GRCh38 allele, left-aligned VCF representation."""

    chrom: str  # 1-22, X, Y, MT
    pos: int
    ref: str
    alt: str

    @property
    def vcf_id(self) -> str:
        """gnomAD-style id, e.g. 17-43045712-T-C."""
        return f"{self.chrom}-{self.pos}-{self.ref}-{self.alt}"

    @property
    def gnomad_chrom(self) -> str:
        return "M" if self.chrom == "MT" else self.chrom

    @property
    def spdi(self) -> str:
        # SPDI is 0-based. For a VCF allele we can write it without the shared
        # anchor base trimmed; NCBI accepts either form.
        return f"{GRCH38_ACCESSIONS[self.chrom]}:{self.pos - 1}:{self.ref}:{self.alt}"

    def as_dict(self) -> dict:
        return asdict(self) | {"vcf_id": self.vcf_id, "assembly": "GRCh38"}


class VariantInputError(ValueError):
    """The input could not be parsed or resolved to a GRCh38 allele."""


def classify_input(text: str) -> str:
    t = text.strip()
    if _VCF_RE.match(t):
        return "vcf"
    if _RSID_RE.match(t):
        return "rsid"
    if _HGVS_RE.match(t):
        return "hgvs"
    return "unknown"


def parse_vcf(text: str) -> Allele:
    m = _VCF_RE.match(text.strip())
    if not m:
        raise VariantInputError(f"Not a VCF-style variant: {text!r}")
    chrom = m["chrom"].upper()
    chrom = "MT" if chrom in {"M", "MT"} else chrom
    return Allele(chrom, int(m["pos"]), m["ref"].upper(), m["alt"].upper())


def _spdi_str(x: dict) -> str:
    return f"{x['seq_id']}:{x['position']}:{x['deleted_sequence']}:{x['inserted_sequence']}"


async def _spdi_to_allele(spdi: str) -> Allele:
    data = (await http().get_json(f"{NCBI_VARIATION}/spdi/{quote(spdi, safe='')}/vcf_fields"))["data"]
    chrom = _ACCESSION_TO_CHROM.get(data["chrom"])
    if chrom is None:
        raise VariantInputError(f"{spdi} is not on a GRCh38 primary-assembly chromosome")
    return Allele(chrom, int(data["pos"]), data["ref"], data["alt"])


async def _from_hgvs(hgvs: str) -> list[Allele]:
    try:
        ctx = await http().get_json(f"{NCBI_VARIATION}/hgvs/{quote(hgvs, safe='')}/contextuals")
    except UpstreamError as e:
        raise VariantInputError(f"NCBI could not parse HGVS {hgvs!r}: {e}") from e
    spdis = ctx.get("data", {}).get("spdis") or []
    if not spdis:
        raise VariantInputError(f"NCBI could not resolve HGVS {hgvs!r}")
    first = _spdi_str(spdis[0])
    if spdis[0]["seq_id"] in _ACCESSION_TO_CHROM:
        return [await _spdi_to_allele(first)]
    # Transcript-level SPDI: project it onto the GRCh38 chromosome.
    eq = await http().get_json(f"{NCBI_VARIATION}/spdi/{quote(first, safe='')}/all_equivalent_contextual")
    genomic = [x for x in eq["data"]["spdis"] if x["seq_id"] in _ACCESSION_TO_CHROM]
    if not genomic:
        raise VariantInputError(f"{hgvs!r} has no GRCh38 genomic placement")
    return [await _spdi_to_allele(_spdi_str(genomic[0]))]


async def _from_rsid(rsid: str) -> list[Allele]:
    try:
        snp = await http().get_json(f"{NCBI_VARIATION}/refsnp/{rsid.lower().removeprefix('rs')}")
    except UpstreamError as e:
        raise VariantInputError(f"dbSNP lookup failed for {rsid}: {e}") from e
    alleles: list[Allele] = []
    for placement in snp.get("primary_snapshot_data", {}).get("placements_with_allele", []):
        if not placement.get("is_ptlp"):  # primary top-level placement = GRCh38 chromosome
            continue
        for a in placement["alleles"]:
            s = a["allele"]["spdi"]
            if s["deleted_sequence"] != s["inserted_sequence"]:
                alleles.append(await _spdi_to_allele(_spdi_str(s)))
    if not alleles:
        raise VariantInputError(f"{rsid} has no GRCh38 alternate alleles in dbSNP (merged or withdrawn?)")
    return alleles


async def check_reference(allele: Allele) -> None:
    """Raise if REF does not match the GRCh38 reference (catches mistyped coordinates)."""
    chrom = "chrM" if allele.chrom == "MT" else f"chr{allele.chrom}"
    res = await http().get_json(
        UCSC_SEQUENCE,
        params={"genome": "hg38", "chrom": chrom, "start": allele.pos - 1, "end": allele.pos - 1 + len(allele.ref)},
    )
    actual = (res.get("dna") or "").upper()
    if actual != allele.ref:
        raise VariantInputError(
            f"REF mismatch for {allele.vcf_id}: GRCh38 has {actual or '?'!r} at {allele.chrom}:{allele.pos}, "
            f"not {allele.ref!r}. Check the coordinate and assembly (GRCh37 positions differ)."
        )


async def resolve(text: str) -> list[Allele]:
    """Resolve any supported input to one or more GRCh38 alleles."""
    kind = classify_input(text)
    if kind == "vcf":
        allele = parse_vcf(text)
        await check_reference(allele)
        return [allele]
    if kind == "hgvs":
        return await _from_hgvs(text.strip())
    if kind == "rsid":
        return await _from_rsid(text.strip())
    raise VariantInputError(
        f"Unrecognized variant {text!r}. Use VCF-style (17-43045712-T-C), an rsID (rs80357258), "
        "or HGVS with a versioned RefSeq accession (NM_007294.4:c.5266dup)."
    )


async def resolve_one(text: str) -> Allele:
    """Resolve to exactly one allele; multi-allelic rsIDs must be disambiguated."""
    alleles = await resolve(text)
    if len(alleles) > 1:
        raise VariantInputError(
            f"{text} matches {len(alleles)} alleles: {', '.join(a.vcf_id for a in alleles)}. "
            "Pass one of these VCF-style ids instead."
        )
    return alleles[0]
