#!/bin/bash
# make_bundle.sh – build the student data of the short-sleep project (data/examples/fnss/) from the raw inputs.
#
#   bash make_bundle.sh DATA WORK
#
# DATA is a folder made by the scripts of ../fetch (see ../README.md): toga/, ncbi/, human/, sleep/, tree/,
# alphafold/, alphamissense/, gnomad/, genebass/, reads/. WORK is a new folder for the work: the scripts are
# copied into WORK/scripts, DATA is linked as WORK/data, the tables go to WORK/outputs/tables and the student
# data to WORK/outputs/student/fnss (everything of data/examples/fnss/ but its README.md).
#
# Needs Python 3 with Biopython, and R with the packages ape and jsonlite (for example
# R_LIBS_USER=/path/to/library). Takes about four minutes, most of it the alignments of steps 02 and 03.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
DATA=$(cd "$1" && pwd)
WORK=$2
mkdir -p "$WORK/scripts"
WORK=$(cd "$WORK" && pwd)
cp "$HERE"/*.py "$HERE"/*.R "$WORK/scripts/"
ln -sfn "$DATA" "$WORK/data"
cd "$WORK"
python3 scripts/02_toga_calls.py          # outputs/tables/toga_site_calls.tsv
python3 scripts/03_ncbi_calls.py          # outputs/tables/ncbi_site_calls.tsv
python3 scripts/04_species_table.py       # outputs/tables/species_site_calls.tsv
python3 scripts/06_human_variation.py     # outputs/tables/gnomad_*.tsv (and AlphaMissense summaries)
python3 scripts/07_sleep_data.py          # outputs/tables/sleep_species.tsv
# the orders and common names of the sleep species (NCBI taxonomy), if not fetched yet
if [ ! -f data/sleep/sleep_species_taxonomy.json ]; then
  python3 "$HERE/../fetch/10_sleep_taxonomy.py" data outputs/tables/sleep_species.tsv
fi
python3 scripts/09_ukb_exomes.py          # outputs/tables/ukb_exome_*.tsv
python3 scripts/14_student_bundle.py      # outputs/student/fnss/ (all but the tree)
Rscript scripts/15_student_tree.R         # outputs/student/fnss/tree/mammals.nwk
echo "The student data are in $WORK/outputs/student/fnss"
