#!/bin/bash
# 05_alphamissense.sh OUT – AlphaMissense (Cheng et al. 2023; Zenodo record 8208688; CC BY-NC-SA 4.0), cut down to
# the five genes. The two files are streamed and filtered on the way, never stored whole (1.2 and 2.5 GB, compressed):
#   OUT/alphamissense/AlphaMissense_aa_substitutions.fnss.tsv           the five canonical UniProt proteins
#   OUT/alphamissense/AlphaMissense_isoforms_aa_substitutions.fnss.tsv  every Ensembl transcript of the five genes
#                                                                       that the isoform file has (GRM1 isoform beta,
#                                                                       for R889W, and SIK3's MANE transcript, for N783Y)
#   OUT/ensembl/lookup_GENE.json   the Ensembl genes, with their transcripts (Ensembl REST API), to choose them
# The comment lines of the files (copyright and licence) are kept.
set -euo pipefail
OUT=$1
mkdir -p "$OUT/alphamissense" "$OUT/ensembl"
Z=https://zenodo.org/records/8208688/files
curl -sSfL --retry 3 "$Z/AlphaMissense_aa_substitutions.tsv.gz" | gzip -dc |
  awk -F'\t' '/^#/ || $1 == "uniprot_id" || $1 == "Q9C0J9" || $1 == "P08588" || $1 == "Q6W5P4" || $1 == "Q13255" || $1 == "Q9Y2K2"' \
  > "$OUT/alphamissense/AlphaMissense_aa_substitutions.fnss.tsv"
for g in BHLHE41 ADRB1 NPSR1 GRM1 SIK3; do
  curl -sSf --retry 3 -o "$OUT/ensembl/lookup_$g.json" \
    "https://rest.ensembl.org/lookup/symbol/homo_sapiens/$g?expand=1;content-type=application/json"
  sleep 1
done
python3 -c 'import json, sys
for f in sys.argv[1:]:
    for t in json.load(open(f))["Transcript"]:
        print(t["id"])' "$OUT"/ensembl/lookup_*.json > "$OUT/ensembl/transcripts.txt"
curl -sSfL --retry 3 "$Z/AlphaMissense_isoforms_aa_substitutions.tsv.gz" | gzip -dc |
  awk -F'\t' 'NR == FNR {want[$1]; next} /^#/ || $1 == "transcript_id" || substr($1, 1, index($1, ".") - 1) in want' \
  "$OUT/ensembl/transcripts.txt" - > "$OUT/alphamissense/AlphaMissense_isoforms_aa_substitutions.fnss.tsv"
wc -l "$OUT"/alphamissense/*.tsv
