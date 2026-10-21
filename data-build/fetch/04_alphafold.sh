#!/bin/bash
# 04_alphafold.sh OUT – the AlphaFold models (AlphaFold DB, version 6) of the five human proteins, by UniProt accession:
#   OUT/alphafold/AF-Q9C0J9-F1-model_v6.pdb (DEC2), AF-P08588-… (ADRB1), AF-Q6W5P4-… (NPSR1),
#   AF-Q13255-… (GRM1, isoform alpha), AF-Q9Y2K2-… (SIK3, UniProt's canonical isoform of 1,321 amino acids)
# (The server refuses curl's own user agent, hence -A.)
set -euo pipefail
mkdir -p "$1/alphafold"
for a in Q9C0J9 P08588 Q6W5P4 Q13255 Q9Y2K2; do
  curl -sSf --retry 3 -A "Mozilla/5.0" -o "$1/alphafold/AF-$a-F1-model_v6.pdb" "https://alphafold.ebi.ac.uk/files/AF-$a-F1-model_v6.pdb"
done
ls -l "$1/alphafold"
