#!/bin/bash
# 03_human.sh OUT – the human reference proteins (NCBI RefSeq, E-utilities):
#   OUT/human/human_refseq.fasta  DEC2 NP_110389.1, ADRB1 NP_000675.1, NPSR1 isoform A NP_997055.1,
#                                 GRM1 isoforms alpha NP_001264993.1 and beta NP_001264994.1, SIK3 NP_001353615.1
#                                 (the MANE Select proteins, and mGluR1b), and NPSR1 isoform G NP_001287864.1
#                                 (the protein of the transcript TOGA used for NPSR1)
set -euo pipefail
mkdir -p "$1/human"
IDS=NP_110389.1,NP_000675.1,NP_997055.1,NP_001264993.1,NP_001264994.1,NP_001353615.1,NP_001287864.1
curl -sSf --retry 3 -o "$1/human/human_refseq.fasta" \
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=$IDS&rettype=fasta&retmode=text"
grep -c ">" "$1/human/human_refseq.fasta"
