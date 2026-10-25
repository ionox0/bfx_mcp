# ACMG/AMP variant classification criteria (summary)

Source: Richards et al. 2015, *Genet Med* 17:405-424, with ClinGen Sequence
Variant Interpretation (SVI) working group refinements noted inline. This is a
condensed reference for reasoning, not a substitute for the primary
publications or gene-specific ClinGen VCEP specifications, which override the
generic rules (often with different frequency cut-offs and strength levels).

The "Evidence source" column says which tool in this server, if any, supplies
the data. Many criteria need case-level or functional data that no public
aggregate source provides; say so rather than guessing.

## Pathogenic criteria

| Code | Strength | Meaning | Evidence source |
|---|---|---|---|
| PVS1 | Very strong | Null variant (nonsense, frameshift, canonical +/-1,2 splice, initiation codon, single/multi-exon deletion) in a gene where loss of function is a known disease mechanism. Apply via the ClinGen PVS1 decision tree (Abou Tayoun 2018): downgrade for NMD escape (last exon / within 50 nt of last junction), truncations removing <10% of protein, alternative transcripts, rescue. | `annotate_variant` (consequence, LOFTEE `lof`/`lof_filter`/`lof_flags`, exon number); `get_gene_constraint` (LoF intolerance as supporting context for mechanism, not proof) |
| PS1 | Strong | Same amino-acid change as an established pathogenic variant, different nucleotide change. | ClinVar (other alleles at the codon) |
| PS2 | Strong | De novo, maternity and paternity confirmed, in a patient with the disease and no family history. | Case data - not available |
| PS3 | Strong | Well-established functional studies show a damaging effect. | Literature - not available |
| PS4 | Strong | Prevalence in affected individuals significantly increased vs controls. | Case-control data - not available |
| PM1 | Moderate | Located in a mutational hot spot or well-established functional domain without benign variation. | Partly: protein position from `annotate_variant`; domain knowledge |
| PM2 | Moderate (SVI 2020: apply at **Supporting**) | Absent from controls, or extremely low frequency for recessive disorders. | `get_population_frequency` |
| PM3 | Moderate | For recessive disorders, detected in trans with a pathogenic variant. | Case data - not available |
| PM4 | Moderate | Protein length change from in-frame indel in a non-repeat region, or stop-loss. | `annotate_variant` |
| PM5 | Moderate | Novel missense change at a residue where a different missense change is established pathogenic. | ClinVar (other alleles at the codon) |
| PM6 | Moderate | Assumed de novo without confirmation of parentage. | Case data - not available |
| PP1 | Supporting | Co-segregation with disease in multiple affected family members. | Family data - not available |
| PP2 | Supporting | Missense in a gene with a low rate of benign missense variation where missense is a common disease mechanism. | `get_gene_constraint` (missense Z / o/e as context) |
| PP3 | Supporting (up to Strong with calibrated predictors) | Multiple lines of computational evidence support a deleterious effect. | `annotate_variant` (REVEL, AlphaMissense, CADD, SpliceAI) |
| PP4 | Supporting | Patient phenotype or family history highly specific for a single-gene disease. | Case data - not available |
| PP5 | Supporting | Reputable source reports pathogenic without shareable evidence. **SVI (2018) recommends not using PP5.** | - |

## Benign criteria

| Code | Strength | Meaning | Evidence source |
|---|---|---|---|
| BA1 | Stand-alone | Allele frequency > 5% in a large population database (SVI: use grpmax filtering allele frequency, with a small list of exceptions). | `get_population_frequency` (`grpmax_faf95`) |
| BS1 | Strong | Allele frequency greater than expected for the disorder (threshold depends on prevalence, penetrance and allelic heterogeneity; usually gene-specific). | `get_population_frequency` |
| BS2 | Strong | Observed in a healthy adult with full penetrance expected at an early age (recessive: homozygous; dominant: heterozygous; X-linked: hemizygous). | `get_population_frequency` (homozygote / hemizygote counts as partial evidence) |
| BS3 | Strong | Well-established functional studies show no damaging effect. | Literature - not available |
| BS4 | Strong | Lack of segregation in affected family members. | Family data - not available |
| BP1 | Supporting | Missense in a gene where primarily truncating variants cause disease. | Gene knowledge |
| BP2 | Supporting | Observed in trans with a pathogenic variant for a fully penetrant dominant disorder, or in cis with a pathogenic variant in any inheritance pattern. | Case data - not available |
| BP3 | Supporting | In-frame indel in a repetitive region without known function. | `annotate_variant` |
| BP4 | Supporting (up to Strong with calibrated predictors) | Multiple lines of computational evidence suggest no impact. | `annotate_variant` |
| BP5 | Supporting | Found in a case with an alternate molecular basis for disease. | Case data - not available |
| BP6 | Supporting | Reputable source reports benign without shareable evidence. **SVI (2018) recommends not using BP6.** | - |
| BP7 | Supporting | Synonymous (or deep intronic, per SVI splicing guidance) variant with no predicted splice impact, at a non-conserved position. | `annotate_variant` (consequence, SpliceAI) |

## Calibrated computational thresholds (ClinGen SVI)

- **REVEL** (Pejaver et al. 2022, single predictor, missense only):
  PP3 Supporting >= 0.644, Moderate >= 0.773, Strong >= 0.932;
  BP4 Supporting <= 0.290, Moderate <= 0.183, Strong <= 0.016.
  Values between 0.290 and 0.644 are uninformative.
- **SpliceAI** (Walker et al. 2023, splicing): max delta score >= 0.2 supports
  PP3 for splicing; <= 0.1 supports BP4/BP7.
- Use one calibrated predictor for PP3/BP4 rather than counting several
  correlated tools as independent lines of evidence. Do not combine PP3 with
  PVS1 for the same mechanism.

## Combining criteria (Richards 2015, Table 5)

**Pathogenic** if any of:
1. 1 Very strong AND (>=1 Strong OR >=2 Moderate OR 1 Moderate + 1 Supporting OR >=2 Supporting)
2. >=2 Strong
3. 1 Strong AND (>=3 Moderate OR 2 Moderate + >=2 Supporting OR 1 Moderate + >=4 Supporting)

**Likely pathogenic** if any of:
1. 1 Very strong AND 1 Moderate
2. 1 Strong AND 1-2 Moderate
3. 1 Strong AND >=2 Supporting
4. >=3 Moderate
5. 2 Moderate AND >=2 Supporting
6. 1 Moderate AND >=4 Supporting

**Benign**: BA1, or >=2 Strong benign.

**Likely benign**: 1 Strong + 1 Supporting benign, or >=2 Supporting benign.

**Uncertain significance**: criteria for none of the above are met, or
benign and pathogenic criteria contradict.

The Bayesian point system (Tavtigian et al. 2018/2020: Supporting = 1,
Moderate = 2, Strong = 4, Very strong = 8 points; benign criteria negative;
P >= 10, LP 6-9, VUS 0-5, LB -1 to -6, B <= -7) is an equivalent and
often clearer way to combine evidence at non-default strengths.
