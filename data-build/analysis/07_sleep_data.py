"""Step 07: daily sleep time of mammal species, linked to the species with FNSS-site calls and to the
mammal tree.

Two compilations:
 - Phylogeny of Sleep database (McNamara et al. 2008), as used by Capellini et al. (2008): records of
   terrestrial species with EEG recordings of at least 12 h, excluding studies of juveniles or infants
   (the 'restricted' criteria of Capellini et al.); species value = median of the qualifying records.
 - msleep (Savage & West 2007, as distributed with ggplot2): one value per species.

A sleep species is linked to the sequence data by its binomial (with a few synonyms); if the species
itself has no genome, a congener is used only when all sequenced congeners carry the same residue
(recorded as match = 'genus').

Outputs: outputs/tables/sleep_species.tsv, outputs/tables/sleep_site_residues.tsv
"""
import os
import re
import statistics
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fnsslib as L

S = os.path.join(L.DATA, "sleep")
# names used by the sleep compilations -> names used by NCBI taxonomy for the genomes
SYN = {"Canis lupus familiaris": "Canis lupus", "Mustela putorius": "Mustela putorius",
       "Chlorocebus aethiops": "Chlorocebus aethiops", "Osphranter rufus": "Osphranter rufus",
       "Sapajus apella": "Sapajus apella", "Pusa caspica": "Pusa caspica",
       "Equus asinus": "Equus asinus", "Nannospalax leucodon": "Nannospalax leucodon"}

# ---------------------------------------------------------------- Phylogeny of Sleep database
pmap = {r["reported_name"]: r for r in L.read_tsv(os.path.join(S, "pos_species_map.tsv"))}
db = L.read_tsv(os.path.join(S, "phylogeny_of_sleep_db.tsv"))
recs = defaultdict(list)
n_q = 0
for r in db:
    name = r["SpeciesName_Reported"].strip()
    sp = pmap.get(name, {}).get("species", ".")
    if sp == ".":
        continue
    try:
        tst = float(r["Total_daily_sleep"])
    except ValueError:
        continue
    q = (r["type"] == "Land" and r["EEG"] == "Yes"
         and r["Twenty_four_hour"] in ("2", "3", "At least 24 hours")
         and r["Summary_age_class"] not in ("Juvenile", "Infant", "Infant & Juvenile", "Adult & Juvenile",
                                            "Infant Juvenile & Adult", "Adult & Infant"))
    if q:
        n_q += 1
        recs[sp].append(tst)
pos = {sp: (statistics.median(v), len(v), min(v), max(v)) for sp, v in recs.items()}

# ---------------------------------------------------------------- msleep
mmap = {r["name"].strip(): r for r in L.read_tsv(os.path.join(S, "msleep_species_map.tsv"))}
ms = {}
for r in L.read_tsv(os.path.join(S, "msleep.tsv")):
    m = mmap[r["name"].strip()]
    if m["species"] == ".":
        continue
    ms[m["species"]] = (float(r["sleep_total"]), r["name"].strip(), r["order"], r["bodywt"], r["vore"])

species = sorted(set(pos) | set(ms))
rows = []
for sp in species:
    p = pos.get(sp)
    m = ms.get(sp)
    rows.append({
        "species": sp, "species_genome_name": SYN.get(sp, sp),
        "pos_total_sleep_h": round(p[0], 2) if p else "", "pos_n_records": p[1] if p else 0,
        "pos_min": round(p[2], 2) if p else "", "pos_max": round(p[3], 2) if p else "",
        "msleep_total_sleep_h": m[0] if m else "", "msleep_name": m[1] if m else "",
        "msleep_order": m[2] if m else "", "msleep_bodywt_kg": m[3] if m else "", "msleep_diet": m[4] if m else "",
    })
L.write_tsv(os.path.join(L.TAB, "sleep_species.tsv"), rows)
both = [r for r in rows if r["pos_total_sleep_h"] != "" and r["msleep_total_sleep_h"] != ""]
print(f"Phylogeny of Sleep: {n_q} qualifying records, {len(pos)} species; msleep: {len(ms)} species; both: {len(both)}")
if len(both) > 3:
    x = [float(r["pos_total_sleep_h"]) for r in both]
    y = [float(r["msleep_total_sleep_h"]) for r in both]
    print("correlation of the two compilations (Pearson r):", round(statistics.correlation(x, y), 3))

# ---------------------------------------------------------------- link to residues
calls = L.read_tsv(os.path.join(L.TAB, "species_site_calls.tsv"))
by_site_sp = defaultdict(dict)
by_site_genus = defaultdict(lambda: defaultdict(set))
for c in calls:
    if not c["call"] or c["evidence"] == "conflict":
        continue
    by_site_sp[c["variant_id"]][c["species"]] = c
    by_site_genus[c["variant_id"]][c["species"].split()[0]].add(c["call"])
out = []
for vid in sorted(by_site_sp):
    for r in rows:
        g = r["species_genome_name"]
        c = by_site_sp[vid].get(g)
        if c:
            res, match, src = c["call"], "species", c["species"]
        else:
            gen = g.split()[0]
            s = by_site_genus[vid].get(gen, set())
            if len(s) == 1:
                res, match = list(s)[0], "genus"
                src = ";".join(sorted(k for k in by_site_sp[vid] if k.split()[0] == gen))
            else:
                continue
        out.append({"variant_id": vid, "species": r["species"], "residue": res, "match": match,
                    "sequence_species": src, "pos_total_sleep_h": r["pos_total_sleep_h"],
                    "msleep_total_sleep_h": r["msleep_total_sleep_h"], "msleep_bodywt_kg": r["msleep_bodywt_kg"]})
L.write_tsv(os.path.join(L.TAB, "sleep_site_residues.tsv"), out)
for vid in sorted(by_site_sp):
    o = [x for x in out if x["variant_id"] == vid]
    from collections import Counter
    print(vid, "sleep species with a residue:", len(o), "(species-level", sum(1 for x in o if x["match"] == "species"),
          ")", dict(Counter(x["residue"] for x in o)))
