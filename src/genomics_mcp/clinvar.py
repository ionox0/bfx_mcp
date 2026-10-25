"""ClinVar via NCBI E-utilities.

Lookup is by GRCh38 position (esearch on chrpos38), then exact allele matching
against the VCF-style coordinates in each record's VCV XML. Searching ClinVar
by HGVS text directly is fuzzy - one dup query returns ~10 unrelated records -
so we never trust it for identity.
"""

import os
import xml.etree.ElementTree as ET
from typing import Any

from .http import http
from .variants import Allele

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

REVIEW_STARS = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,
    "criteria provided, single submitter": 1,
    "no assertion criteria provided": 0,
    "no classification provided": 0,
    "no classification for the individual variant": 0,
}

MAX_SUBMISSIONS = 25


def _params(**kw: Any) -> dict[str, Any]:
    p = {"db": "clinvar", "tool": "genomics-mcp", **kw}
    if key := os.environ.get("NCBI_API_KEY"):
        p["api_key"] = key
    return p


def _text(e: ET.Element | None, path: str) -> str | None:
    if e is None:
        return None
    found = e.find(path)
    return found.text.strip() if found is not None and found.text else None


def _conditions(classification: ET.Element) -> list[str]:
    names = []
    for trait in classification.iterfind("ConditionList/TraitSet/Trait"):
        name = _text(trait, "Name/ElementValue[@Type='Preferred']")
        if name and name not in names:
            names.append(name)
    return names


def _vcf_match(record: ET.Element, allele: Allele) -> bool:
    # Only the record's own allele: haplotype and compound records also contain
    # SequenceLocations for each member allele, which must not match.
    for loc in record.iterfind("ClassifiedRecord/SimpleAllele/Location/SequenceLocation"):
        if (
            loc.get("Assembly") == "GRCh38"
            and loc.get("positionVCF") == str(allele.pos)
            and loc.get("referenceAlleleVCF") == allele.ref
            and loc.get("alternateAlleleVCF") == allele.alt
            and (loc.get("Chr") or "").upper() == allele.chrom
        ):
            return True
    return False


def _parse_record(va: ET.Element) -> dict[str, Any]:
    cr = va.find("ClassifiedRecord")
    germ = cr.find("Classifications/GermlineClassification") if cr is not None else None
    review = _text(germ, "ReviewStatus")
    out: dict[str, Any] = {
        "accession": f"{va.get('Accession')}.{va.get('Version')}",
        "variation_id": va.get("VariationID"),
        "name": va.get("VariationName"),
        "variation_type": va.get("VariationType"),
        "record_last_updated": va.get("DateLastUpdated"),
        "germline_classification": _text(germ, "Description"),
        "review_status": review,
        "review_stars": REVIEW_STARS.get(review or "", 0),
        "classification_breakdown": _text(germ, "Explanation"),
        "date_last_evaluated": germ.get("DateLastEvaluated") if germ is not None else None,
        "conditions": _conditions(germ) if germ is not None else [],
        "url": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{va.get('VariationID')}/",
    }
    # Somatic classifications live alongside germline ones on newer records.
    for tag, key in (("SomaticClinicalImpact", "somatic_clinical_impact"), ("OncogenicityClassification", "oncogenicity")):
        node = cr.find(f"Classifications/{tag}") if cr is not None else None
        if node is not None:
            out[key] = {"description": _text(node, "Description"), "review_status": _text(node, "ReviewStatus")}

    submissions = []
    for ca in cr.iterfind("ClinicalAssertionList/ClinicalAssertion") if cr is not None else []:
        acc = ca.find("ClinVarAccession")
        cls = ca.find("Classification")
        submissions.append({
            "scv": f"{acc.get('Accession')}.{acc.get('Version')}" if acc is not None else None,
            "submitter": acc.get("SubmitterName") if acc is not None else None,
            "classification": _text(cls, "GermlineClassification"),
            "review_status": _text(cls, "ReviewStatus"),
            "date_last_evaluated": cls.get("DateLastEvaluated") if cls is not None else None,
        })
    out["submission_count"] = len(submissions)
    out["submissions"] = submissions[:MAX_SUBMISSIONS]
    if len(submissions) > MAX_SUBMISSIONS:
        out["submissions_truncated"] = True
    return out


async def lookup(allele: Allele) -> list[dict[str, Any]]:
    """ClinVar records whose GRCh38 VCF allele exactly matches `allele`."""
    end = allele.pos + max(len(allele.ref), len(allele.alt))
    term = f"{allele.chrom}[chr] AND {allele.pos}:{end}[chrpos38]"
    search = await http().get_json(f"{EUTILS}/esearch.fcgi", params=_params(term=term, retmode="json", retmax=200))
    ids = search["esearchresult"]["idlist"]
    if not ids:
        return []
    records = []
    for i in range(0, len(ids), 50):
        xml = await http().request(
            "GET", f"{EUTILS}/efetch.fcgi",
            params=_params(id=",".join(ids[i : i + 50]), rettype="vcv", is_variationid="true"),
        )
        for va in ET.fromstring(xml).iter("VariationArchive"):
            if _vcf_match(va, allele):
                records.append(_parse_record(va))
    return records


async def fetch_by_variation_ids(ids: list[str]) -> list[tuple[ET.Element, dict[str, Any]]]:
    """Raw VCV elements plus parsed summaries, for dataset building."""
    out = []
    for i in range(0, len(ids), 50):
        xml = await http().request(
            "GET", f"{EUTILS}/efetch.fcgi",
            params=_params(id=",".join(ids[i : i + 50]), rettype="vcv", is_variationid="true"),
        )
        out.extend((va, _parse_record(va)) for va in ET.fromstring(xml).iter("VariationArchive"))
    return out


def vcf_allele_of(va: ET.Element) -> Allele | None:
    """The GRCh38 VCF allele of a simple-allele VCV record, if ClinVar provides one."""
    for loc in va.iterfind("ClassifiedRecord/SimpleAllele/Location/SequenceLocation"):
        if loc.get("Assembly") == "GRCh38" and loc.get("positionVCF"):
            chrom = (loc.get("Chr") or "").upper()
            return Allele("MT" if chrom == "MT" else chrom, int(loc.get("positionVCF")),
                          loc.get("referenceAlleleVCF"), loc.get("alternateAlleleVCF"))
    return None


async def search_ids(term: str, retmax: int = 500) -> list[str]:
    res = await http().get_json(f"{EUTILS}/esearch.fcgi", params=_params(term=term, retmode="json", retmax=retmax))
    return res["esearchresult"]["idlist"]
