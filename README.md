# Seven mutations, four hundred mammals – the natural short sleep project

**A guided undergraduate research project of the Chronomic Medicine Group, in a workbench that runs in the
browser.** Seven rare mutations in five genes – *DEC2*/*BHLHE41* (P384R, Y362H), *ADRB1* (A187V), *NPSR1*
(Y206H), *GRM1* (S458A, R889W) and *SIK3* (N783Y) – have been reported to cause familial natural short sleep, and
every paper that reported one called the changed amino acid "conserved". The student tests that claim in the
proteins of 200–450 mammals from two independent sources (NCBI RefSeq orthologs and TOGA), finds where on the
dated mammal tree the variable positions changed, asks whether those changes go with how long species sleep
(with and without phylogenetic regression), places the mutations in AlphaFold models, looks at them in gnomAD and
the UK Biobank exomes, and checks, in the raw sequencing reads, a Brazilian porcupine that carries one of the
mutations. It needs no programming experience: the guide teaches the Linux terminal, `kalign`, R (with `ape` and
`phylolm`) and the 3D viewer as the project needs them.

The site is built from the [Bioinformatics workbench](https://github.com/comparativechrono/bioinformatics-workbench)
of the same group: the same notebooks (R and Python), terminal, Files tab, Databases tab, 3D viewer, IDEs and AI
assistant – with the programs and tabs the project does not need left out, and a guide in the left pane in place
of the workbench's help (the help is still there, behind **Help**).

---

## Contents

1. [The project in brief](#the-project-in-brief)
2. [The guide](#the-guide)
3. [Put it online](#put-it-online)
4. [Before the student starts](#before-the-student-starts)
5. [Editing the guide](#editing-the-guide)
6. [The data](#the-data)
7. [How it was tested](#how-it-was-tested)
8. [What is different from the Bioinformatics workbench](#what-is-different-from-the-bioinformatics-workbench)
9. [Files](#files)
10. [Credits and licences](#credits-and-licences)

---

## The project in brief

**Question.** The evidence that these mutations shorten sleep comes from one or two families each (sometimes one
person) and from mice made to carry them. Each paper supports its mutation with a claim of conservation drawn
from a handful of species. Mammals sleep from about 3 to 20 hours a day: if one amino acid can shorten sleep, has
evolution changed the same positions, and do the species that changed them sleep differently?

**Plan.** Eleven stages, about one a week (10–12 weeks part-time). Stages 1–7 are the core; 8, 9 and 10 each add a
separate line of evidence; each stage ends with "going further" ideas, so the student can make the project her own.

| Stage | Title | The student learns |
| --- | --- | --- |
| 0 | Start here; the science; the data; the reading list; glossary | how sleep is controlled, the seven mutations and what conservation means; how to keep a logbook |
| 1 | Seven mutations | shell commands, paths, pipes; tables in text files; FASTA; `seqtk`, `paste`, `awk`; a `while read` loop; UniProt and Europe PMC |
| 2 | Hundreds of mammals | orthologs from NCBI and TOGA; searching and counting FASTA files; `comm`; `for` loops; why numbering differs between species |
| 3 | Line them up | multiple alignment with `kalign`; first R notebook: functions, matrices, `which`, `table`, `apply`; checking an alignment locally |
| 4 | Conserved or not? | an R function run in a loop; `merge`; comparing two sources; conservation profiles; finding a fault in the data |
| 5 | Where the changes happened | trees in R (`ape`); ancestral states with `ace`; counting changes along branches |
| 6 | How long mammals sleep | exploring a table; histograms, scatter and box plots; correlation; where published numbers come from |
| 7 | Do the amino acids track sleep? | `aggregate`; t-tests; non-independence of species; phylogenetic regression (`phylolm`, Pagel's λ, Brownian motion); reading a null result |
| 8 | Inside the proteins | AlphaFold and pLDDT; the 3D viewer; fixed-width PDB columns; isoform numbering (SIK3 N783 = N735 in UniProt) |
| 9 | The mutations in people | gnomAD frequencies and carriers; UK Biobank exome associations (Genebass) and a forest plot; multiple testing; AlphaMissense |
| 10 | A porcupine with a human mutation | sequencing reads, strands and reverse complements; a binomial test for heterozygosity; ruling out a duplicated gene |
| 11 | Write it up | figures, methods with versions, results with uncertainty, discussion; turning notebooks into scripts and a `run_all.sh` |

As far as we know the analysis is new: the conservation claims of all seven mutations tested across hundreds of
mammals and two sources, the history of the variable positions on a dated tree, a phylogenetic test against
species' sleep, and a heterozygous NPSR1 Y206H porcupine confirmed in its raw reads. **A supervisor's answer key is
not in this repository** (it is public, and the student would find it): the results that the guide's *Check* boxes
give are there to confirm each step, and the interpretation is left to the student.

## The guide

The left pane has a row of tabs: **Start**, the stages **1–11**, and **Help** (the manual of the tools). Each
stage page has: what the student will learn; the science where needed; numbered steps with grey code boxes
(**Run in the terminal** or **Add to notebook**); green **Check** boxes with the numbers she should get; *Your
turn* tasks with hidden hints and answers; *Think* questions; a logbook box; *going further*; and a checklist,
whose ticks are kept in the browser and fill in the stage's tab. Each stage with R has its own notebook
(`~/notebooks/stageNN-….ipynb`, made by a button), and the student keeps a logbook, `~/project/logbook.md`.

Early stages give the code and explain it line by line; later ones ask her to write more of it herself, with the
code behind *Show the answer*. Her results go into `~/project/results` and `~/project/figures` (alignments into
`~/project/alignments`), and stage 11 turns the work into scripts that rebuild everything.

## Put it online

1. In this repository on GitHub: *Settings → Pages → Build and deployment*: Source **Deploy from a branch**,
   branch **main**, folder **/ (root)**. Save.
2. After a minute the site is at **https://comparativechrono.github.io/chronomics-project-26-2/**. The first time
   a browser opens it, the page loads once more by itself: it installs a small service worker (`isolation-sw.js`)
   that makes the page cross-origin isolated (needed only by PLINK 2, which this project does not use, but kept
   from the workbench).

The repository must be public for GitHub Pages on a free plan. To try it on your own computer: `python3 serve.py
8000`, then open `http://localhost:8000` (opening `index.html` as a file does not work).

## Before the student starts

- **Computer.** A laptop or desktop with a current Chrome, Edge, Firefox or Safari; 8 GB of memory is comfortable.
- **Network.** Besides this site, the browser must reach `webr.r-wasm.org`, `repo.r-wasm.org` and
  `cran.r-universe.dev` (R and its packages – `ape` and `phylolm` are fetched the first time stage 5 or 7 loads
  them), `cdn.jsdelivr.net` (only if she opens a Python notebook), and the services the guide uses:
  `rest.uniprot.org`, `www.ebi.ac.uk` (Europe PMC, ENA, InterPro, BLAST), `alphafold.ebi.ac.uk`, `files.rcsb.org`
  and `data.rcsb.org` (the 3D viewer), `rest.ensembl.org` and `eutils.ncbi.nlm.nih.gov`.
- **Her files live in her browser** – not on a server: the same computer and browser each time, and **My files →
  Download all my files** once a week (the guide says so). A private window, or clearing the browser's data, loses
  them. The site keeps its files apart from the Bioinformatics workbench's and from the PER3 project's
  (`storePrefix: 'chrono26b'` in `assets/js/config.js`), although all are on `comparativechrono.github.io`.
- **The AI assistant** (optional) needs a key of her own (Google Gemini has a free tier); the guide tells her to
  check what it says and to declare its use. Decide whether you want her to use it.
- **Public data only.** All the data are public; the human data are summary counts and association results, not
  anyone's genome. The guide tells her never to put personal or patient data into the workbench or an AI tool.

## Editing the guide

The page is built from `build/`. After any change there, run:

```bash
python3 build/make_page.py        # writes index.html; warns about links to pages that do not exist
```

- **Guide pages**: `build/help/00-home.html` (Start), `00a-science.html`, `00b-data.html`, `00c-reading.html`,
  `00d-glossary.html` and the stages `00s01-mutations.html` … `00s11-writeup.html`. The first line of each holds
  its id, title, group (`The project`), `stage: N`, and `tab: yes` for the page that a stage's tab opens.
  Checklists are `<ul class="checklist"><li data-check="unique-id">…</li></ul>`; the classes of the boxes
  (`callout check`, `think`, `log`, `science`, `task`, `details.hint`) are styled in `assets/css/guide.css`, and
  the tabs, checklists and stage-notebook buttons are run by `assets/js/guide.js`.
- **Tool help**: the other pages of `build/help/` (the workbench's, trimmed to this project).
- **Starter notebooks**: `python3 build/make_starters.py` writes `data/starter/*.ipynb` (`start-here.ipynb` opens
  first, with the sleep table).
- **Settings**: `assets/js/config.js` (title, the cards of the Databases tab, R and AI settings).
- **Programs**: after a change in `assets/vendor/`, run `python3 build/vendor_sums.py`.

## The data

`data/examples/fnss/` (about 5 MB) holds everything the project needs, described in
[data/examples/fnss/README.md](data/examples/fnss/README.md) and in the guide (*The project data*): the seven
mutations and the six human proteins; the orthologous proteins of the five genes in mammals from NCBI (203–266
species per gene) and from TOGA (299–401); the sleep times of 101 mammals from the Phylogeny of Sleep database and
msleep; the dated mammal tree of Upham et al. (2019) for the 485 species; AlphaFold models of the five proteins;
the mutations in gnomAD v4.1, the UK Biobank exomes (Genebass) and AlphaMissense; and the 57 read pairs of a
Brazilian porcupine around *NPSR1* codon 206. How each file was made – the scripts and the addresses – is in
[data-build/](data-build/). The other folders of `data/examples/` are small examples kept from the workbench for
its help pages.

## How it was tested

On 21 October 2026, in Chromium (Playwright): every code box of the guide (stages 1–11, the first session, and the
help page on alignments) was run in order in the workbench – terminal boxes in the terminal as the *Run in the
terminal* button runs them, R boxes with `Rscript` and, for stages 3–5, through **Add to notebook** in the stage's
notebook (tables, plots) – and gave the numbers of the *Check* boxes. The same code was run natively (kalign 3.3.5,
R with `ape` and `phylolm`): the results are the same, except where the browser's older kalign (3.3.1) aligns a few
sequences differently – the DEC2 alignment is 820 columns wide instead of 777, and five more NPSR1 sequences from
TOGA are misaligned at position 206 (the window check sets them aside, and stage 5 shows one of them, the North
American porcupine, by another route). The *Check* boxes give the browser's numbers. The stage tabs, checklists,
stage notebooks, the 3D viewer's buttons (AlphaFold models, labels, `site`) and the Databases buttons (UniProt,
ENA, Europe PMC) were checked. The analysis behind the guide was done first, independently, in Python and R, with a
pairwise alignment of every protein to the human one: in all 2,869 species and positions that both routes call,
they give the same amino acid. The AI assistant was not tested with a real key.

## What is different from the Bioinformatics workbench

- **Left out**: the Galaxy, Networks and Genome browser tabs, and the programs for read QC, mapping, RNA-seq,
  single-cell, somatic and metagenomic analyses, peak calling, gene finding and HMMER (their files are not in
  `assets/vendor/`, and `assets/js/project-tools.js` takes them out of the terminal); the example data that went
  with them; most Databases cards (`databases:` in `config.js`).
- **Kept**: seqtk, kalign, MUSCLE, FastTree, seq-align, LASTZ, samtools, bcftools, bgzip/tabix, minimap2, bedtools,
  gffread, jq, wgsim, PLINK 2, the GNU tools, Python and R.
- **Added**: the guide (stage tabs, checklists, stage notebooks, logbook), the project's data and its
  `data-build/`, and short-sleep examples in the Databases tab.

## Files

```
index.html             the page (built by build/make_page.py)
isolation-sw.js        the service worker that makes the page cross-origin isolated (keep it next to index.html)
assets/js/             the workbench; guide.js (the guide), project-tools.js (the programs of this project), config.js
assets/css/            styles; guide.css
assets/vendor/         bundled programs and libraries, with their licences (see THIRD_PARTY.md)
data/examples/fnss/    the project's data, and its README.md
data/examples/         small example files kept from the workbench
data/starter/          the starter notebooks
data-build/            how the project's data were made
build/                 make_page.py, help/ (the guide and the help pages), template.html, citations.json,
                       make_starters.py, vendor_sums.py
serve.py               a small local web server (optional)
THIRD_PARTY.md         the third-party software and data, their licences and sources
```

## Credits and licences

Built from the Bioinformatics workbench of the Chronomic Medicine Group. The bundled programs and libraries keep
their own licences – see [THIRD_PARTY.md](THIRD_PARTY.md), and keep `assets/vendor/licenses/` with the site. The
data are public (sources and terms in THIRD_PARTY.md and the data README); note that the Zenodo record of the AlphaMissense
predictions gives them as CC BY-NC-SA 4.0 (non-commercial use), although DeepMind's repository now gives CC BY 4.0. The studies and resources to cite are listed in the guide (*The project
data*). The guide and the site's own code and text are by the Chronomic Medicine Group.
