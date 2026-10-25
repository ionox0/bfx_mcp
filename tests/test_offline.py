"""Tests that need no network."""

import asyncio
import json
import os
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from genomics_mcp import clinvar
from genomics_mcp.variants import Allele, VariantInputError, classify_input, parse_vcf

sys.path.insert(0, str(Path(__file__).parent.parent / "evals"))
from grade import grade_reasoning, grade_retrieval, grade_tool_use  # noqa: E402


@pytest.mark.parametrize("text,expected", [
    ("17-43045712-T-C", Allele("17", 43045712, "T", "C")),
    ("chr17:43045712:T:C", Allele("17", 43045712, "T", "C")),
    ("17:43045712 T>C", Allele("17", 43045712, "T", "C")),
    ("chrX-100-a-g", Allele("X", 100, "A", "G")),
    ("chrM-3243-A-G", Allele("MT", 3243, "A", "G")),
    ("7-117559590-ATCT-A", Allele("7", 117559590, "ATCT", "A")),
])
def test_parse_vcf(text, expected):
    assert parse_vcf(text) == expected


@pytest.mark.parametrize("text,kind", [
    ("rs80357258", "rsid"),
    ("NM_007294.4:c.5266dup", "hgvs"),
    ("NC_000017.11:g.43045712T>C", "hgvs"),
    ("NM_007294.4(BRCA1):c.5266dup", "hgvs"),
    ("17-43045712-T-C", "vcf"),
    ("BRCA1 c.5266dup", "unknown"),
])
def test_classify_input(text, kind):
    assert classify_input(text) == kind


def test_parse_vcf_rejects_garbage():
    with pytest.raises(VariantInputError):
        parse_vcf("17-abc-T-C")


def _vcv(record_type: str, chrom="7", pos="117559590", ref="ATCT", alt="A") -> ET.Element:
    loc = (f'<Location><SequenceLocation Assembly="GRCh38" Chr="{chrom}" positionVCF="{pos}" '
           f'referenceAlleleVCF="{ref}" alternateAlleleVCF="{alt}"/></Location>')
    body = (f"<SimpleAllele>{loc}</SimpleAllele>" if record_type == "simple"
            else f"<Haplotype><SimpleAllele>{loc}</SimpleAllele></Haplotype>")
    return ET.fromstring(f'<VariationArchive Accession="VCV1" Version="1"><ClassifiedRecord>{body}</ClassifiedRecord></VariationArchive>')


def test_clinvar_matches_only_the_records_own_allele():
    f508del = Allele("7", 117559590, "ATCT", "A")
    assert clinvar._vcf_match(_vcv("simple"), f508del)
    assert not clinvar._vcf_match(_vcv("haplotype"), f508del)
    assert not clinvar._vcf_match(_vcv("simple", alt="G"), f508del)


def test_grade_reasoning():
    g = grade_reasoning({"classification": "Pathogenic"}, {"classification": "Likely pathogenic"})
    assert g == {"three_tier": 1.0, "exact": 0.0, "severe_error": 0.0, "within_one": 1.0}
    g = grade_reasoning({"classification": "Benign"}, {"classification": "Likely pathogenic"})
    assert g["severe_error"] == 1.0 and g["three_tier"] == 0.0
    assert grade_reasoning({"classification": "Benign"}, {})["three_tier"] == 0.0


@pytest.mark.parametrize("expected,answer,ok", [
    ({"type": "number", "answer": 0.0119312, "rel_tol": 0.02}, {"answer": "1.19e-2"}, True),
    ({"type": "number", "answer": 0.0119312, "rel_tol": 0.02}, {"answer": "0.0125"}, False),
    ({"type": "number", "answer": 0.836, "abs_tol": 0.0015}, {"answer": "REVEL = 0.836"}, True),
    ({"type": "integer", "answer": 58}, {"answer": "58"}, True),
    ({"type": "hgvsp", "answer": "p.Gln1756ProfsTer74"}, {"answer": "p.(Gln1756Profs*74)"}, True),
    ({"type": "hgvsp", "answer": "p.Phe508del"}, {"answer": "ENSP00000003084.6:p.Phe508del"}, True),
    ({"type": "string", "answer": "7-117559590-ATCT-A"}, {"answer": "chr7:117559590:ATCT:A"}, True),
    ({"type": "string", "answer": "frameshift_variant"}, {"answer": "Frameshift variant"}, True),
    ({"type": "string_set", "answer": ["17-43045712-T-A", "17-43045712-T-C"]},
     {"answer": "17-43045712-T-C, 17-43045712-T-A"}, True),
    ({"type": "string_set", "answer": ["17-43045712-T-A", "17-43045712-T-C"]}, {"answer": "17-43045712-T-C"}, False),
    ({"type": "invalid", "answer": "invalid"}, {"answer": "", "input_invalid": True}, True),
    ({"type": "invalid", "answer": "invalid"}, {"answer": "0.001", "input_invalid": False}, False),
    ({"type": "number", "answer": 0.5, "abs_tol": 0.01}, {"answer": "0.5", "input_invalid": True}, False),
])
def test_grade_tool_use(expected, answer, ok):
    assert grade_tool_use(expected, answer)["correct"] == float(ok)


def test_grade_retrieval():
    exp = {"classification": "Conflicting classifications of pathogenicity", "review_stars": 1,
           "has_conflict": True, "accession": "VCV000012345"}
    ans = {"classification": "conflicting classifications of pathogenicity", "review_stars": 1,
           "has_conflicting_submissions": True, "clinvar_accession": "VCV000012345.3"}
    assert grade_retrieval(exp, ans)["all_fields"] == 1.0
    assert grade_retrieval(exp, {**ans, "review_stars": 2})["all_fields"] == 0.0


def test_hide_clinvar_removes_the_tool():
    code = (
        "import asyncio, json; from mcp import Client; from genomics_mcp.server import mcp\n"
        "async def m():\n"
        "    async with Client(mcp) as c: print(json.dumps([t.name for t in (await c.list_tools()).tools]))\n"
        "asyncio.run(m())"
    )
    def tools(hide: str) -> list[str]:
        env = {**os.environ, "GENOMICS_MCP_HIDE_CLINVAR": hide}
        out = subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True, check=True)
        return json.loads(out.stdout.strip().splitlines()[-1])
    assert "get_clinvar_classification" in tools("")
    assert "get_clinvar_classification" not in tools("1")
