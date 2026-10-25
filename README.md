# genomics-mcp

An MCP server that gives AI agents ClinVar, gnomAD and Ensembl VEP, plus an
eval for variant interpretation in three tracks. GRCh38 only. For research and
education: output is not a clinical interpretation.

The original teaching notebook server is in git history (`git show main:server.py`).

## What the server exposes

| Kind | Name | Source | Returns |
|---|---|---|---|
| tool | `normalize_variant` | NCBI Variation Services, UCSC | GRCh38 left-aligned VCF allele(s) from VCF-style ids, rsIDs or HGVS. Checks REF against the genome |
| tool | `annotate_variant` | Ensembl VEP | Consequence, HGVS c./p. on MANE Select, SIFT, PolyPhen, CADD, REVEL, AlphaMissense, SpliceAI, LOFTEE. ClinVar/HGMD annotations are stripped |
| tool | `get_population_frequency` | gnomAD v4 (or v3) | Joint / exome / genome AC, AN, AF, homozygotes, grpmax FAF95, per-ancestry AF, filters |
| tool | `get_gene_constraint` | gnomAD v4 | pLI, LOEUF, o/e and Z scores, MANE transcript |
| tool | `get_clinvar_classification` | ClinVar (E-utilities) | Aggregate classification, review stars, conditions, per-submitter calls. Exact allele match only |
| resource | `genomics://acmg-criteria` | — | ACMG/AMP criteria, combining rules, calibrated REVEL/SpliceAI thresholds, mapped to the tools above |
| resource | `genomics://sources` | — | Provenance and limitations |
| resource template | `genomics://variant/{variant}` | all | One JSON evidence bundle, e.g. `@genomics:genomics://variant/17-43045712-T-C` |
| prompt | `interpret_variant` | — | ACMG-style workflow (`/mcp__genomics__interpret_variant`) |

Every tool accepts any supported notation and normalizes it first, so an agent
can pass `NM_000492.4:c.1521_1523del`, `rs113993960` or `7-117559590-ATCT-A`.
A multi-allelic rsID is rejected with the list of alleles to choose from.

| Environment variable | Effect |
|---|---|
| `GENOMICS_MCP_HIDE_CLINVAR=1` | Remove the ClinVar tool (used by the blinded eval track) |
| `NCBI_API_KEY` | Raise the NCBI limit from 3 to 10 requests/s |
| `GENOMICS_MCP_CACHE` | SQLite response cache path (default `~/.cache/genomics-mcp/http.sqlite`) |
| `GENOMICS_MCP_CACHE_TTL_S` | Cache lifetime, default 7 days; `inf` never expires |

## Use it from Claude Code

```bash
claude mcp add genomics -- uv run --directory /path/to/mcp_server genomics-mcp
```

Then ask, for example, "Is NM_000059.4:c.7007G>A likely to affect splicing?",
or attach `@genomics:genomics://acmg-criteria`.

## Upstream behaviour worth knowing

- **Ensembl REST is slow and flaky** (10-60 s per VEP call, frequent 500/503).
  The HTTP layer retries with backoff and caches every success on disk.
- **gnomAD** rate-limits aggressively; requests are paced to 1/s per process.
- **ClinVar text search is fuzzy**, so lookups search by GRCh38 position and then
  match the VCF allele of each record's own SimpleAllele. Haplotype and compound
  records that merely contain the allele are not returned.
- HGVS and rsIDs are normalized with NCBI because Ensembl's variant recoder was
  returning 500s during development.

## Eval

Three tracks, one runner, programmatic grading throughout.

| Track | ClinVar | Task | Gold | Headline metric |
|---|---|---|---|---|
| `reasoning` | hidden | Classify a variant (5-tier) with ACMG criteria from primary evidence | ClinVar expert-panel / practice-guideline classification | 3-tier match (P/LP, VUS, B/LB); also 5-tier exact, P↔B flips |
| `tool_use` | visible | Factual questions: AF, homozygotes, grpmax group, consequence, HGVS p., REVEL, SpliceAI, LOEUF, normalization, GRCh37 traps | Computed from the APIs at build time | correct |
| `retrieval` | visible | Report ClinVar's classification, stars, conflict status, accession | ClinVar at build time | all fields correct |

```bash
uv run python evals/build_dataset.py              # rebuild all tracks (slow: Ensembl)
uv run python evals/build_dataset.py retrieval    # rebuild one track
uv run --extra eval python evals/run.py reasoning --limit 5      # pilot
uv run --extra eval python evals/run.py tool_use --model claude-sonnet-5-5 --variant sonnet
node <claude-api skill>/shared/evals/report/build-report-lite.mjs .claude/hillclimb/genomics-reasoning/
```

`build_dataset.py` writes the cases to `evals/data/*.jsonl` and a review table
to `evals/cases.md`. It also stores every upstream response it used in
`evals/data/cache.sqlite` with no expiry. The runner starts the server against
that same cache, so the agent sees the data the gold answers came from even
after ClinVar and gnomAD update. Delete the cache and rebuild to refresh.

Known limits of the reasoning track:

- Expert panels use case, segregation and functional data that public
  aggregate sources do not have, so 100% is not attainable from these tools.
  Many expert-panel VUS and LP calls depend on that data.
- Well-known variants may be memorized. `meta.date_last_evaluated` is recorded
  per case so results can be split by evaluation date.
- Gene-specific VCEP rules (for example BRCA1/2 ENIGMA frequency thresholds)
  are not given to the agent unless you add them.

## Tests

```bash
uv run pytest tests        # offline: parsing, graders, ClinVar allele matching, ClinVar hiding
uv run python scripts/smoke.py [variant ...]   # live: every tool through an MCP client session
```
