"""Feature switches, so evals can compare server versions side by side.

GENOMICS_MCP_FEATURES is a comma-separated list; unset means all features,
"none" means the original v0.1 behaviour.

  codon     get_clinvar_codon_variants tool (PS1/PM5 evidence)
  pvs1      PVS1 context (NMD prediction, last exon, LOFTEE detail) in annotate_variant
  intronic  intronic HGVS normalization via Ensembl when NCBI declines
"""

import os

ALL = frozenset({"codon", "pvs1", "intronic"})

_raw = os.environ.get("GENOMICS_MCP_FEATURES")
if _raw is None:
    ENABLED = ALL
elif _raw.strip().lower() in {"", "none"}:
    ENABLED = frozenset()
else:
    ENABLED = frozenset(f.strip() for f in _raw.split(",")) & ALL


def on(name: str) -> bool:
    return name in ENABLED
