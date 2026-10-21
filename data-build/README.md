# How the project's data were made

The files in `data/examples/fnss/` are small pieces of large public data sets. These scripts made them, on
21 October 2026, from the sources named in each. They are here so that the data can be checked and remade (for
instance when the sources are updated) – students do not need them.

There are two steps. The scripts in `fetch/` download the raw inputs into a folder; `analysis/make_bundle.sh`
then makes the student data from that folder, with the scripts of the project's analysis. They need a Linux or
macOS computer with network access, `curl`, Python 3 with Biopython, R with the packages ape, jsonlite and
ggplot2, and, for the porcupine's reads, `vdb-dump` of the SRA Toolkit. Each fetch script writes into the folder
it is given:

```bash
mkdir raw
bash fetch/00_inputs.sh raw          # human/variants.tsv, sleep/*_species_map.tsv, tree/tree_name_synonyms.tsv (from inputs/)
bash fetch/01_toga.sh raw            # toga/: the codon alignments of the five genes, overview.table.tsv, READMEs
python3 fetch/02_ncbi.py raw         # ncbi/: orthologs_*.json, proteins_*.faa, taxonomy.json (after 01)
bash fetch/03_human.sh raw           # human/human_refseq.fasta
bash fetch/04_alphafold.sh raw       # alphafold/AF-*-F1-model_v6.pdb
bash fetch/05_alphamissense.sh raw   # alphamissense/*.fnss.tsv, ensembl/lookup_*.json (streams 3.7 GB)
python3 fetch/06_gnomad.py raw       # gnomad/gnomad_v4.1_*.json
python3 fetch/07_genebass.py raw     # genebass/*.json
python3 fetch/08_sleep.py raw        # sleep/Complete_Sleep_Dataset_FINAL.xml, phylogeny_of_sleep_db.tsv, msleep.tsv
bash fetch/09_tree.sh raw            # tree/MamPhy_DNAonly_4098sp_NDexp_MCC.tre
python3 fetch/11_scan_sra.py raw     # reads/coendou_SRR10331604_NPSR1_reads.txt (some hours)
bash analysis/make_bundle.sh raw work   # work/outputs/student/fnss/: data/examples/fnss/ but its README.md
```

`fetch/10_sleep_taxonomy.py` looks up the order, family and common name of each sleep species in NCBI's taxonomy
(`sleep/sleep_species_taxonomy.json`). It needs the list of sleep species that the analysis makes, so
`make_bundle.sh` runs it when the file is missing. `inputs/` holds the four files made by hand: the table of the
mutations, from the papers, and three tables of names – the names of the two sleep compilations, and of a few
genomes in the tree, matched to species.

## What `make_bundle.sh` does

`analysis/` has copies of the scripts of the project's full analysis that the student data need, run in this
order (about four minutes):

| Script | Makes |
| --- | --- |
| `02_toga_calls.py` | the amino acid of every TOGA projection at each mutated position, from a pairwise alignment (BLOSUM62) of its protein to the human one, with the identity of the 21 positions around it (`fnsslib.py` has the shared code) |
| `03_ncbi_calls.py` | the same for every NCBI protein of every ortholog, keeping for each species the isoform best aligned there |
| `04_species_table.py` | one row per species and position, both sources together |
| `06_human_variation.py` | the mutations in gnomAD v4.1: counts, by ancestry, and the other changes seen at the same codons |
| `07_sleep_data.py` | the sleep of each species: the median of the records of the Phylogeny of Sleep database that meet the criteria of Capellini et al. 2008 (EEG, at least 12 hours, adults, on land), and msleep |
| `09_ukb_exomes.py` | the mutations in Genebass: carriers, effect and an approximate 95% confidence interval for six sleep questions; gene-level tests |
| `14_student_bundle.py` | the student files: per gene and source, one protein per mammal species – the human first; of NCBI's isoforms whose alignment around the mutation passes the check, the one covering most of the human protein, and of TOGA's projections the best-aligned one; names `Genus_species` throughout; the sleep table, with orders and common names from NCBI's taxonomy; the AlphaFold models; the population tables; AlphaMissense at the seven positions; the porcupine's reads, each once |
| `15_student_tree.R` | `tree/mammals.nwk`: the tree of Upham et al. 2019 cut down to the species of the files, tips renamed to match them |

A species' amino acid at a position counts only where the alignment around it is good: an amino acid, not a gap;
at least half of the 21 positions around it identical to the human ones; at most four of them unaligned, and at
most four amino acids inserted among them.

## Checked

On 21 October 2026 the fetch scripts 00–04 and 06–10 were run again and gave files identical, byte for byte, to
those the student data were made from: TOGA's alignments and tables; NCBI's orthologs, proteins and taxonomy;
the human proteins; the AlphaFold models; gnomAD; all 41 Genebass files; the sleep workbook, its conversion and
msleep; the tree; and the taxonomy of the sleep species. `05_alphamissense.sh` was not run in full (it streams
3.7 GB): its filters were run on the start of both files followed by the rows kept, and gave exactly the rows
kept. `11_scan_sra.py` was run on two ranges of the run (rows 1–200,000, and 230,164,641–231,959,781, which holds
one of the reads) and found what the original scan found there; reading the whole run takes some hours. A full
scan finds the same reads, but may write them in another order (the original file was put together from three
passes over parts of the run), and the order of the lines of `reads/coendou_NPSR1_reads.txt` follows it.
`make_bundle.sh`, run on the original raw files and on the newly fetched ones, made student data identical to
`data/examples/fnss/` (all but its README.md).

NCBI, Ensembl and Genebass are updated from time to time: run later, the scripts may fetch slightly different
records, and the student data may then differ a little.

## Sources

| Files | Source | Version / date |
| --- | --- | --- |
| `toga/*.fasta.gz`, `overview.table.tsv`, `README*.txt` | TOGA (Kirilenko et al. 2023), human hg38 reference: `genome.senckenberg.de/download/TOGA/human_hg38_reference/` (`MultipleCodonAlignments/`) | the published release; README of 17 January 2025 |
| `ncbi/orthologs_*.json`, `proteins_*.faa` | NCBI Datasets API v2: `gene/id/{id}/orthologs` and `gene/download` (FASTA_PROTEIN), for genes 79365 (BHLHE41), 153 (ADRB1), 387129 (NPSR1), 2911 (GRM1), 23387 (SIK3) | 21 October 2026 |
| `ncbi/taxonomy.json` | NCBI E-utilities, `efetch db=taxonomy` (old ids of merged taxa kept, with `merged_into`) | 21 October 2026 |
| `human/human_refseq.fasta` | NCBI E-utilities, `efetch db=protein` (RefSeq) | 21 October 2026 |
| `alphafold/*.pdb` | AlphaFold DB, `alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v6.pdb` | version 6 |
| `alphamissense/*.fnss.tsv` | AlphaMissense (Cheng et al. 2023), Zenodo record 8208688: `AlphaMissense_aa_substitutions.tsv.gz`, `AlphaMissense_isoforms_aa_substitutions.tsv.gz`; CC BY-NC-SA 4.0 on that record (DeepMind's repository now gives CC BY 4.0) | 2023 |
| `ensembl/lookup_*.json` | Ensembl REST API, `lookup/symbol/homo_sapiens/{gene}?expand=1` | release 116 |
| `gnomad/gnomad_v4.1_*.json` | gnomAD GraphQL API (`gnomad.broadinstitute.org/api`), dataset `gnomad_r4` | v4.1 |
| `genebass/*.json` | Genebass API (`main.genebass.org/api`; Karczewski et al. 2022: 394,841 UK Biobank exomes) | 21 October 2026 |
| `sleep/Complete_Sleep_Dataset_FINAL.xml`, `phylogeny_of_sleep_db.tsv` | Phylogeny of Sleep database (McNamara et al. 2008), `www.bu.edu/phylogeny/datasets/Complete_Sleep_Dataset_FINAL.xml` | workbook of 4 December 2007 |
| `sleep/msleep.tsv` | R package ggplot2, data set `msleep` (Savage & West 2007) | ggplot2 4.0.3 |
| `sleep/sleep_species_taxonomy.json` | NCBI E-utilities, `esearch` and `efetch db=taxonomy` | 21 October 2026 |
| `tree/MamPhy_DNAonly_4098sp_NDexp_MCC.tre` | Upham, Esselstyn & Jetz 2019: `github.com/n8upham/MamPhy_v1`, `_DATA/MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre` | v1 |
| `reads/coendou_SRR10331604_NPSR1_reads.txt` | NCBI SRA, run SRR10331604 (DNA Zoo, BioProject PRJNA512907, BioSample SAMN13072004; Illumina NovaSeq 6000, 463,919,562 read pairs), read with `vdb-dump` of the SRA Toolkit 3.4.1 | – |
| `inputs/*` | made for the project: the mutations from the papers (PubMed IDs in the table); names matched by hand, with NCBI's taxonomy | 21 October 2026 |
