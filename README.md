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

## Results: `interpret_variant` eval

The `interpret` track tests the server's `interpret_variant` MCP prompt end to
end. Each case's user message is the prompt text a user would get from
`/mcp__genomics__interpret_variant`, the agent reads `genomics://acmg-criteria`
itself through a resource tool, and ClinVar's per-variant tool is hidden. It
uses the reasoning cases: 60 variants with ClinVar expert-panel or
practice-guideline classifications, 12 per tier across 22 genes. A disjoint
60-case held-out set (`reasoning_holdout`, different seed) checks changes
developed against the main set's failures. Agent: `claude-opus-5-5` at high
effort, one run per case.

Only the classification label is graded against gold. The criteria evidence
text, summary and "not evaluable" notes are not scored, so a correct tier
could rest on wrong reasoning. `self-consistent` checks only that the tier
follows from the stated criteria codes and strengths under the ACMG point
system.

| | Original server, main set | Improved server, main set | Original server, held-out |
|---|---|---|---|
| **P/LP sensitivity** (P/LP gold called P/LP) | 14/24 (58%) | **17/24 (71%)** | 14/24 (58%) |
| **P↔B flips** | 0/60 | 0/60 | 0/57 |
| B/LB correct | 21/24 | 20/24 | 19/21 |
| VUS correct | 11/12 | 11/12 | 9/12 |
| 3-tier match | 46/60 | 48/60 | 42/57 |
| 5-tier exact | 30/60 | 31/60 | 28/57 |
| Self-consistent | 60/60 | 59/60 | 54/57 |
| Median latency / tool calls | 32 s / 5.7 | 26 s / 5.6 | 29 s / 5.3 |
| Cost | $6.76 | $6.66 | $6.32 |

"Improved server" adds three changes, each switchable with
`GENOMICS_MCP_FEATURES` (see `src/genomics_mcp/features.py`):
`get_clinvar_codon_variants` (other ClinVar alleles at the same codon, for
PS1/PM5, excluding the variant's own record), a `pvs1_context` block in
`annotate_variant` (NMD prediction, last exon, LOFTEE detail), and intronic
HGVS normalization through Ensembl when NCBI declines.

What the results show:

- **The main failure is missed pathogenic variants, not flips.** With the
  original server the agent called 10 of 24 P/LP variants VUS and never
  called a variant Pathogenic. It never placed a pathogenic variant in a benign
  tier or the reverse. Most misses were missense variants where the expert
  panel likely relied on functional, segregation or case data. In its notes
  the agent most often asked for same-codon ClinVar evidence (PS1/PM5) and
  gene-specific hotspot lists (PM1).
- **The codon tool fixed 3 misses on the main set.** CDKL5, DICER1 and LDLR
  missense variants moved from VUS to Likely pathogenic after the agent found
  other pathogenic missense changes at the same residue and applied PM5. One
  benign TP53 variant regressed from Likely benign to VUS; it is the only answer
  whose tier contradicts its own criteria.
- **Intronic normalization removed workarounds.** With the original server,
  11 of 60 cases (mostly splice variants) failed normalization and the agent
  guessed genomic coordinates, taking up to 14 tool calls and 251 s. With the
  fix there were no tool errors and median latency fell from 32 s to 26 s.
- **PVS1 context changed no outcome.** For the FBN1 frameshift it confirmed
  the agent's moderate PVS1: the stop is in the last exon and escapes NMD.
- **The gain is within noise.** +3 of 24 P/LP cases sits inside the 95%
  intervals, and the main set is the one the changes were developed against.
  The held-out run of the improved server is pending: the original-server
  held-out run stopped at 57/60 cases and the improved-server run did not
  start, because the API account ran out of credit.

The full report, with method, confusion matrices, breakdowns by consequence
and gene, tool usage and every case's criteria and tool calls, is
`evals/results/interpret.html`. Regenerate it with
`uv run python evals/report.py --out evals/results/interpret.html`. Its
transcript links point into `.claude/hillclimb/`, which is gitignored, so they
only work locally.

## Tests

```bash
uv run pytest tests        # offline: parsing, graders, ClinVar allele matching, ClinVar hiding
uv run python scripts/smoke.py [variant ...]   # live: every tool through an MCP client session
```
