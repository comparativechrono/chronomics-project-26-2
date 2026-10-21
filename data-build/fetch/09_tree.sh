#!/bin/bash
# 09_tree.sh OUT – the dated tree of the mammals of Upham, Esselstyn & Jetz (2019): the DNA-only, node-dated
# (exponential prior) maximum clade credibility tree of 4,098 species, from the authors' repository on GitHub.
#   OUT/tree/MamPhy_DNAonly_4098sp_NDexp_MCC.tre   (NEXUS; tips named Genus_species_FAMILY_ORDER)
set -euo pipefail
mkdir -p "$1/tree"
curl -sSf --retry 3 -o "$1/tree/MamPhy_DNAonly_4098sp_NDexp_MCC.tre" \
  "https://raw.githubusercontent.com/n8upham/MamPhy_v1/master/_DATA/MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre"
ls -l "$1/tree"
