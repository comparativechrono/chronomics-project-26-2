#!/bin/bash
# 00_inputs.sh OUT – put the hand-made input files (../inputs) where the analysis expects them:
#   OUT/human/variants.tsv             the seven FNSS mutations and NPSR1 N107I, with their papers
#   OUT/sleep/pos_species_map.tsv      names in the Phylogeny of Sleep database -> species
#   OUT/sleep/msleep_species_map.tsv   names in msleep -> species
#   OUT/tree/tree_name_synonyms.tsv    species names of the genomes -> names in the tree of Upham et al. 2019
set -euo pipefail
OUT=$1
IN=$(cd "$(dirname "$0")/../inputs" && pwd)
mkdir -p "$OUT/human" "$OUT/sleep" "$OUT/tree"
cp "$IN/variants.tsv" "$OUT/human/"
cp "$IN/pos_species_map.tsv" "$IN/msleep_species_map.tsv" "$OUT/sleep/"
cp "$IN/tree_name_synonyms.tsv" "$OUT/tree/"
