#!/bin/bash
# 01_toga.sh OUT – TOGA (Kirilenko et al. 2023), human hg38 reference: the multiple codon alignments of the five
# genes (the reference transcript TOGA used for each), the table of assemblies, and TOGA's README.
#   OUT/toga/ENST00000242728.BHLHE41.fasta.gz  ENST00000369295.ADRB1.fasta.gz  ENST00000381539.NPSR1.fasta.gz
#   OUT/toga/ENST00000361719.GRM1.fasta.gz     ENST00000445177.SIK3.fasta.gz
#   OUT/toga/overview.table.tsv  README.txt  README.multipleCodonAlignments.txt
set -euo pipefail
OUT=$1/toga
mkdir -p "$OUT"
BASE=https://genome.senckenberg.de/download/TOGA
for t in ENST00000242728.BHLHE41 ENST00000369295.ADRB1 ENST00000381539.NPSR1 ENST00000361719.GRM1 ENST00000445177.SIK3; do
  curl -sSf --retry 3 -o "$OUT/$t.fasta.gz" "$BASE/human_hg38_reference/MultipleCodonAlignments/$t.fasta.gz"
done
curl -sSf --retry 3 -o "$OUT/overview.table.tsv" "$BASE/human_hg38_reference/overview.table.tsv"
curl -sSf --retry 3 -o "$OUT/README.txt" "$BASE/README.txt"
curl -sSf --retry 3 -o "$OUT/README.multipleCodonAlignments.txt" "$BASE/README.multipleCodonAlignments.txt"
ls -l "$OUT"
