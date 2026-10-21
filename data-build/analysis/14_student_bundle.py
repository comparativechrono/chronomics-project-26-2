"""Step 14: the data of the student project (the workbench chronomics-project-26-2), in outputs/student/fnss/.

 orthologs/<GENE>_ncbi_mammals.fasta   one RefSeq protein per mammal species (the isoform that aligns best
                                      at the FNSS site; GRM1: isoform alpha; GRM1b: isoform beta)
 orthologs/<GENE>_toga_mammals.fasta   the same genes predicted from genome assemblies by TOGA (one per species)
 human/human_proteins.fasta           the human proteins whose numbering the mutations use
 variants.tsv                         the seven FNSS mutations (and N107I for comparison)
 sleep/mammal_sleep.tsv               daily sleep of mammal species, two compilations, linked to the sequences
 tree/mammals.nwk                     the dated mammal tree (Upham et al. 2019) for the species of the files
 structures/                          AlphaFold models of the five proteins
 population/                          gnomAD v4.1, UK Biobank exomes (Genebass), AlphaMissense
 reads/                               reads of the Brazilian porcupine around NPSR1 codon 206
FASTA names are Genus_species, the same in every file, the sleep table and the tree.
"""
import csv
import json
import os
import re
import shutil
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fnsslib as L

OUTD = os.path.join(L.OUT, "student", "fnss")
if os.path.exists(OUTD):
    shutil.rmtree(OUTD)
for d in ("orthologs", "human", "sleep", "tree", "structures", "population", "reads"):
    os.makedirs(os.path.join(OUTD, d), exist_ok=True)

tax = {r["tax_id"]: r for r in json.load(open(os.path.join(L.DATA, "ncbi", "taxonomy.json")))}


def sid(name):
    w = name.replace("[", "").replace("]", "").split()
    return "_".join(w[:2])


def wrap(s, n=60):
    return "\n".join(s[i:i + n] for i in range(0, len(s), n))


# ------------------------------------------------------------------ NCBI orthologs
ncbi = L.read_tsv(os.path.join(L.TAB, "ncbi_site_calls.tsv"))
SITE_OF = {"BHLHE41": "BHLHE41_P384R", "ADRB1": "ADRB1_A187V", "NPSR1": "NPSR1_Y206H", "GRM1": "GRM1_S458A",
           "GRM1b": "GRM1_R889W", "SIK3": "SIK3_N783Y"}
GID = {"BHLHE41": "79365", "ADRB1": "153", "NPSR1": "387129", "GRM1": "2911", "GRM1b": "2911", "SIK3": "23387"}
PREFER = {"Canis_lupus": "Canis lupus familiaris"}
counts = {}
ncbi_species = defaultdict(set)
REF_OF = {"BHLHE41": "NP_110389.1", "ADRB1": "NP_000675.1", "NPSR1": "NP_997055.1", "GRM1": "NP_001264993.1",
          "GRM1b": "NP_001264994.1", "SIK3": "NP_001353615.1"}
POS_OF = {"BHLHE41": 384, "ADRB1": 187, "NPSR1": 206, "GRM1": 458, "GRM1b": 889, "SIK3": 783}
hp_all = L.human_proteins()
for gene, site in SITE_OF.items():
    fa = L.readfa(os.path.join(L.DATA, "ncbi", f"proteins_{GID[gene]}.faa"))
    human = hp_all[REF_OF[gene]]
    rows = [r for r in ncbi if r["variant_id"] == site and r["class"] == "Mammalia" and r["confident"] == "1"]
    best = {}
    for r in rows:
        k = sid(r["species"])
        # the isoform for the student's file: of the isoforms whose call at the site is confident (they never
        # disagree there), the one that covers most of the human protein - the full-length protein rather
        # than a short isoform (for GRM1, isoform alpha where the species has one)
        accs = [x.split(":")[0] for x in r["all_isoforms"].split(";") if x]
        cands = []
        for acc in accs:
            seq = fa[acc]
            m, ins, _ = L.map_positions(human, seq)
            c = L.call_site(human, seq, m, ins, POS_OF[gene])
            if L.confident(c):
                cands.append(((round(float(c["coverage"]), 2), float(c["win_id"]), float(c["global_id"]),
                               acc.startswith("NP_")), acc, c["aa"]))
        cands.sort(reverse=True)
        assert cands and cands[0][2] == r["aa"], (gene, k, r["aa"], cands[:2])
        score = (r["species"] == PREFER.get(k, r["species"]),) + cands[0][0]
        cur = best.get(k)
        if cur is None or score > cur[0]:
            best[k] = (score, dict(r, accession=cands[0][1]))
    with open(os.path.join(OUTD, "orthologs", f"{gene}_ncbi_mammals.fasta"), "w") as fh:
        for k in sorted(best, key=lambda x: (x != "Homo_sapiens", x)):
            r = best[k][1]
            fh.write(f">{k} {r['accession']} {r['common_name']}\n{wrap(fa[r['accession']])}\n")
            ncbi_species[gene].add(k)
    counts[f"{gene}_ncbi"] = len(best)

# ------------------------------------------------------------------ TOGA
ov = L.toga_overview()
toga_species = defaultdict(set)
for gene, (fname, ref_acc) in L.TOGA_FILES.items():
    human = L.human_proteins()[ref_acc]
    seqs = L.toga_sequences(gene)
    calls = [r for r in L.read_tsv(os.path.join(L.TAB, "toga_site_calls.tsv")) if r["gene"] == gene]
    conf = defaultdict(int)
    gid = {}
    for r in calls:
        if r["confident"] == "1":
            conf[(r["assembly"], r["projection"])] += 1
        gid[(r["assembly"], r["projection"])] = float(r["global_id"])
    best = {}
    for asm, projs in seqs.items():
        if asm == "REFERENCE":
            continue
        info = ov.get(asm, {})
        t = tax.get(info.get("Species Taxonomy ID", ""), {})
        k = sid(t.get("name") or info.get("Species", asm))
        for prot, cods, proj in projs:
            score = (conf.get((asm, proj), 0), gid.get((asm, proj), 0))
            if score[0] == 0:
                continue
            if k not in best or score > best[k][0]:
                best[k] = (score, asm, prot, info.get("Common name", ""))
    with open(os.path.join(OUTD, "orthologs", f"{gene}_toga_mammals.fasta"), "w") as fh:
        fh.write(f">Homo_sapiens hg38 human\n{wrap(human)}\n")
        for k in sorted(best):
            _, asm, prot, common = best[k]
            if k == "Homo_sapiens":
                continue
            fh.write(f">{k} {asm} {common}\n{wrap(prot)}\n")
            toga_species[gene].add(k)
    counts[f"{gene}_toga"] = len(best)

# ------------------------------------------------------------------ human proteins and the variants
hp = L.readfa(os.path.join(L.DATA, "human", "human_refseq.fasta"))
names = {"NP_110389.1": "BHLHE41 (DEC2), the MANE Select protein", "NP_000675.1": "ADRB1, the MANE Select protein",
         "NP_997055.1": "NPSR1 isoform A, the MANE Select protein",
         "NP_001264993.1": "GRM1 isoform alpha (mGluR1a), the MANE Select protein",
         "NP_001264994.1": "GRM1 isoform beta (mGluR1b)", "NP_001353615.1": "SIK3, the MANE Select protein (isoform 4)"}
genes_of = {"NP_110389.1": "BHLHE41", "NP_000675.1": "ADRB1", "NP_997055.1": "NPSR1", "NP_001264993.1": "GRM1",
            "NP_001264994.1": "GRM1b", "NP_001353615.1": "SIK3"}
with open(os.path.join(OUTD, "human", "human_proteins.fasta"), "w") as fh:
    for acc in names:
        fh.write(f">{genes_of[acc]} {acc} {names[acc]}\n{wrap(hp[acc])}\n")

gn = {r["variant_id"]: r for r in L.read_tsv(os.path.join(L.TAB, "gnomad_fnss_variants.tsv"))}
vrows = []
for v in L.load_variants():
    vrows.append({"variant": v["variant_id"], "gene": v["gene"] if v["variant_id"] != "GRM1_R889W" else "GRM1b",
                  "protein": v["refseq_protein"], "position": v["pos"], "human_residue": v["ref_aa"],
                  "fnss_residue": v["alt_aa"], "change": v["protein_change"],
                  "rsid": gn.get(v["variant_id"], {}).get("rsid", ""), "study": v["study"], "pmid": v["pmid"],
                  "people": v["human_evidence"], "animal_model": v["animal_model"],
                  "where_in_protein": v["location_reported"], "conservation_claim": v["conservation_claim"],
                  "kind": "familial natural short sleep" if v["category"] == "FNSS" else "common variant, for comparison"})
L.write_tsv(os.path.join(OUTD, "variants.tsv"), vrows)

# ------------------------------------------------------------------ sleep
sl = L.read_tsv(os.path.join(L.TAB, "sleep_species.tsv"))
msm = {r["species"]: r for r in L.read_tsv(os.path.join(L.DATA, "sleep", "msleep_species_map.tsv")) if r["species"] != "."}
pmap = L.read_tsv(os.path.join(L.DATA, "sleep", "pos_species_map.tsv"))
marine = {r["species"] for r in pmap if "marine" in (r["note"] or "")}
marine |= {"Globicephala macrorhynchus", "Phocoena phocoena", "Tursiops truncatus", "Callorhinus ursinus",
           "Halichoerus grypus", "Pusa caspica"}
allseq = set().union(*ncbi_species.values(), *toga_species.values())
srows = []
for r in sl:
    sp = r["species"]
    g = r["species_genome_name"]
    k = sid(g)
    if k in allseq:
        gs, match = k, "species"
    else:
        cands = sorted(x for x in allseq if x.split("_")[0] == k.split("_")[0])
        gs, match = (cands[0], "genus") if cands else ("", "none")
    t = None
    order = r["msleep_order"]
    srows.append({"species": sp, "common_name": msm.get(sp, {}).get("name", "") or "",
                  "order": order, "environment": "marine" if sp in marine else "land",
                  "sleep_pos_h": r["pos_total_sleep_h"], "pos_records": r["pos_n_records"],
                  "sleep_msleep_h": r["msleep_total_sleep_h"], "body_mass_kg": r["msleep_bodywt_kg"],
                  "genome_species": gs, "match": match})
# orders from NCBI taxonomy for every species (msleep uses older orders such as Soricomorpha and Cetacea);
# common names from msleep, or from NCBI where msleep has none
# (data/sleep/sleep_species_taxonomy.json: NCBI esearch + efetch of each sleep species, see data-build)
stax = json.load(open(os.path.join(L.DATA, "sleep", "sleep_species_taxonomy.json")))
for r in srows:
    t = stax.get(r["species"], {})
    if t.get("order"):
        r["order"] = t["order"]
    if not r["common_name"]:
        r["common_name"] = t.get("common", "")
    if r["order"] == "Soricomorpha":
        r["order"] = "Eulipotyphla"
    c = r["common_name"].replace("opposum", "opossum").replace(" american ", " American ")
    r["common_name"] = c[:1].upper() + c[1:]
missing = [r["species"] for r in srows if not r["order"] or not r["common_name"]]
if missing:
    print("sleep species without order or common name:", missing)
L.write_tsv(os.path.join(OUTD, "sleep", "mammal_sleep.tsv"), srows)

# ------------------------------------------------------------------ tree
syn = {r["name"]: r["tree_binomial"] for r in L.read_tsv(os.path.join(L.DATA, "tree", "tree_name_synonyms.tsv"))}
txt = open(os.path.join(L.DATA, "tree", "MamPhy_DNAonly_4098sp_NDexp_MCC.tre")).read()
want = set(allseq) | {r["genome_species"] for r in srows if r["genome_species"]}
to_tip = {}
for k in sorted(want):
    sp = k.replace("_", " ")
    tip = syn.get(sp, k)
    to_tip.setdefault(tip, k)  # first data name wins if two map to one tip
json.dump(to_tip, open(os.path.join(L.WORK, "tree_tip_map.json"), "w"))
counts["tree_wanted"] = len(to_tip)

# ------------------------------------------------------------------ structures, population, reads
for f in os.listdir(os.path.join(L.DATA, "alphafold")):
    if f.endswith(".pdb"):
        shutil.copy(os.path.join(L.DATA, "alphafold", f), os.path.join(OUTD, "structures", f))
for f, new in [("gnomad_fnss_variants.tsv", "gnomad_fnss_variants.tsv"),
               ("gnomad_fnss_by_ancestry.tsv", "gnomad_fnss_by_ancestry.tsv"),
               ("gnomad_same_codon.tsv", "gnomad_other_changes_same_codon.tsv"),
               ("ukb_exome_fnss_variants.tsv", "ukb_exomes_fnss_variants.tsv"),
               ("ukb_exome_gene_burden.tsv", "ukb_exomes_gene_burden.tsv")]:
    shutil.copy(os.path.join(L.TAB, f), os.path.join(OUTD, "population", new))
# AlphaMissense: every substitution at the seven sites
am_main = defaultdict(dict)
for line in open(os.path.join(L.DATA, "alphamissense", "AlphaMissense_aa_substitutions.fnss.tsv")):
    if line.startswith("#") or line.startswith("uniprot_id"):
        continue
    u, pv, s, c = line.rstrip("\n").split("\t")
    am_main[u][pv] = (s, c)
am_iso = defaultdict(dict)
for line in open(os.path.join(L.DATA, "alphamissense", "AlphaMissense_isoforms_aa_substitutions.fnss.tsv")):
    if line.startswith("#") or line.startswith("transcript_id"):
        continue
    f = line.rstrip("\n").split("\t")
    am_iso[f[0].split(".")[0]][f[1]] = (f[2], f[3] if len(f) > 3 else "")
CLS = {"likely_benign": "benign", "likely_pathogenic": "pathogenic"}
arows = []
for v in L.load_variants(include_comparison=False):
    if v["variant_id"] == "GRM1_R889W":
        tab, pos, src = am_iso["ENST00000507907"], v["pos"], "ENST00000507907 (mGluR1b)"
    elif v["variant_id"] == "SIK3_N783Y":
        tab, pos, src = am_iso["ENST00000445177"], v["pos"], "ENST00000445177 (SIK3 MANE)"
    else:
        tab, pos, src = am_main[v["uniprot"]], int(v["uniprot_pos"]), v["uniprot"]
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        key = f"{v['ref_aa']}{pos}{aa}"
        if key in tab:
            s, c = tab[key]
            arows.append({"variant": v["variant_id"], "protein_used": src, "substitution": f"{v['ref_aa']}{v['pos']}{aa}",
                          "new_residue": aa, "am_score": s, "am_class": CLS.get(c, c),
                          "is_fnss_mutation": int(aa == v["alt_aa"])})
L.write_tsv(os.path.join(OUTD, "population", "alphamissense_fnss_sites.tsv"), arows)
rd = os.path.join(L.DATA, "reads", "coendou_SRR10331604_NPSR1_reads.txt")
seen = []
for line in open(rd):
    s = line.strip().split("\t")[-1]
    if s and s not in seen:
        seen.append(s)
with open(os.path.join(OUTD, "reads", "coendou_NPSR1_reads.txt"), "w") as fh:
    for s in seen:
        fh.write(s + "\n")
counts["reads"] = len(seen)
json.dump(counts, open(os.path.join(L.WORK, "student_bundle_counts.json"), "w"), indent=1)
print(counts)
