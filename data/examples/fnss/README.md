# The data of the natural short sleep project

Small pieces of large public data sets – the parts about seven positions in five proteins – so that the project
runs in the browser. The guide's page *The project data* says the same, with the papers to cite.

## The mutations and the human proteins

| File | What it is |
| --- | --- |
| `variants.tsv` | The seven familial natural short sleep (FNSS) mutations, and for comparison NPSR1 N107I, a common variant. One row each: `variant`, `gene`, `protein` (the RefSeq protein whose numbering the mutation uses), `position`, `human_residue`, `fnss_residue`, `change`, `rsid`, `study`, `pmid`, `people`, `animal_model`, `where_in_protein`, `conservation_claim` (what the paper said), `kind`. A dot means "nothing to say". |
| `human/human_proteins.fasta` | The six human proteins: BHLHE41 (DEC2, NP_110389.1), ADRB1 (NP_000675.1), NPSR1 isoform A (NP_997055.1), GRM1 isoform alpha (mGluR1a, NP_001264993.1), GRM1 isoform beta (mGluR1b, NP_001264994.1; named `GRM1b`), SIK3 (NP_001353615.1). All but mGluR1b are the MANE Select proteins. |

## `orthologs/` – the same proteins in other mammals

One protein per species. Every name line starts with the species as `Genus_species`, the same names as in the
sleep table and the tree; the human protein is first.

| Files | What they are |
| --- | --- |
| `GENE_ncbi_mammals.fasta` (BHLHE41, ADRB1, NPSR1, GRM1, GRM1b, SIK3) | The RefSeq proteins of NCBI's orthologs of the human gene in mammals (NCBI Datasets, October 2026). Of a species' isoforms whose alignment to the human protein is good around the mutation (window identity ≥ 0.5 over ±10 positions), the one that covers most of the human protein. `GRM1b` has only isoforms with the mGluR1b tail. Name line: `>Genus_species ACCESSION common name`. 203–266 species per gene. |
| `GENE_toga_mammals.fasta` (BHLHE41, ADRB1, NPSR1, GRM1, SIK3) | The same genes predicted by TOGA (Kirilenko et al. 2023) in placental mammal genome assemblies, from TOGA's multiple codon alignments with the human genome (hg38): one protein per species (the projection with the best-aligned window around the mutation). Name line: `>Genus_species ASSEMBLY common name`. 299–401 species per gene, and the human protein of the transcript TOGA used (for NPSR1 isoform G, NP_001287864.1, 390 amino acids, identical to isoform A up to position 342). `X`: a codon TOGA could not read; `*`: a stop codon. |

**A known fault.** In TOGA's alignment of *BHLHE41* the human reference goes out of frame at codon 432, and so do
all the species' predicted proteins from there to the end: the last 50 amino acids in
`BHLHE41_toga_mammals.fasta` are not real. Positions 362 and 384 are not affected.

## `sleep/`

`mammal_sleep.tsv` – 101 mammal species with a daily sleep time in one or both of two compilations:

| Column | What it is |
| --- | --- |
| `species`, `common_name`, `order` | the species (as the sleep studies name it), its common name (msleep, or NCBI Taxonomy), its order (NCBI Taxonomy) |
| `environment` | `land` or `marine` |
| `sleep_pos_h`, `pos_records` | total sleep (hours per 24 h) in the Phylogeny of Sleep database (McNamara et al. 2008), from the records that meet the criteria of Capellini et al. (2008) – EEG, at least 12 h of recording, adults, not marine – the median of the species' records, and how many records (66 species) |
| `sleep_msleep_h`, `body_mass_kg` | total sleep and body mass from msleep (Savage & West 2007), as shipped with the R package ggplot2 (77 species) |
| `genome_species`, `match` | the species of the sequence files that stands for this one, and whether it is the same species (`species`), one of the same genus (`genus`) or none (`none`) |

## `tree/`

`mammals.nwk` – the dated mammal tree of Upham, Esselstyn & Jetz (2019; DNA-only, node-dated, maximum clade
credibility, 4,098 species), cut down to the 485 species of the files above, with the tips named as in them.
Branch lengths in millions of years.

## `structures/`

The AlphaFold models (AlphaFold DB v6) of the human proteins, by UniProt accession: Q9C0J9 (DEC2), P08588 (ADRB1),
Q6W5P4 (NPSR1), Q13255 (GRM1, isoform alpha), Q9Y2K2 (SIK3, 1,321 amino acids: N783 of the MANE protein is N735
here). The B-factor column is the pLDDT.

## `population/`

| File | What it is |
| --- | --- |
| `gnomad_fnss_variants.tsv` | each mutation in gnomAD v4.1 (exomes + genomes): allele count `ac`, number of alleles `an`, frequency `af`, homozygotes `hom`, carrier frequency |
| `gnomad_fnss_by_ancestry.tsv` | the same by genetic ancestry group |
| `gnomad_other_changes_same_codon.tsv` | the other changes gnomAD has at the same codons |
| `ukb_exomes_fnss_variants.tsv` | each mutation in the 394,841 UK Biobank exomes of Genebass (Karczewski et al. 2022): carriers' alleles `ac`, effect `beta` with an approximate 95% confidence interval (from the effect and its P-value) and `p`, for sleep duration (hours), chronotype (0–3, higher = later), ease of getting up (0–3), insomnia (0–2), daytime napping and daytime dozing (log odds) |
| `ukb_exomes_gene_burden.tsv` | Genebass gene-level tests of the five genes for sleep duration |
| `alphamissense_fnss_sites.tsv` | AlphaMissense scores (Cheng et al. 2023) of all 19 possible changes at each of the seven positions (for R889W and N783Y, of the transcripts with that numbering) |

## `reads/`

`coendou_NPSR1_reads.txt` – the 57 read pairs of SRR10331604 (DNA Zoo, whole-genome Illumina sequencing of a
Brazilian porcupine, *Coendou prehensilis*; 463,919,562 pairs) that contain the 20 bases before or after *NPSR1*
codon 206 on either strand. One pair per line: both 151-base reads, joined.

## Licences and citation

NCBI's data (Datasets, Taxonomy, the SRA) are free to use – NCBI places no restrictions on them; TOGA's data are free to
use with citation of Kirilenko et al. (2023); the Phylogeny of Sleep database and msleep are published research data (msleep
as shipped with ggplot2, MIT licence); the Upham et al. tree comes from the authors' MamPhy_v1 repository (the same trees are
on Dryad, doi:10.5061/dryad.tb03d03, under CC0); AlphaFold DB models are CC BY 4.0; AlphaMissense predictions are CC
BY-NC-SA 4.0 (non-commercial use) on the Zenodo record they come from, though DeepMind's repository now gives CC BY 4.0; gnomAD's summary data are released without restriction on use (gnomAD's terms of use);
Genebass's results are public summary statistics of UK Biobank. Cite the papers listed on the guide's page *The project
data*. How the files were made is in `data-build/` of the
repository.
