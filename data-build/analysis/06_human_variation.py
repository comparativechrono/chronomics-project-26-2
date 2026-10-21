"""Step 06: the FNSS variants in human populations (gnomAD v4.1) and their AlphaMissense scores.

For each FNSS variant: allele counts in gnomAD v4.1 exomes + genomes, overall and by genetic ancestry
group, homozygotes, and the implied number of carriers; every other missense change seen at the same
codon. AlphaMissense: the score of the FNSS substitution, its rank among all possible substitutions in
the protein, and the scores of the residues that other mammals carry at the same position.

Outputs: outputs/tables/gnomad_fnss_variants.tsv, gnomad_fnss_by_ancestry.tsv, gnomad_same_codon.tsv,
         alphamissense_fnss.tsv
"""
import csv
import json
import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fnsslib as L

variants = L.load_variants()
three = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
         "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
         "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*"}
# Ensembl transcript that carries each variant's numbering
TRANSCRIPT = {"BHLHE41_P384R": "ENST00000242728", "BHLHE41_Y362H": "ENST00000242728",
              "ADRB1_A187V": "ENST00000369295", "NPSR1_Y206H": "ENST00000360581",
              "NPSR1_N107I": "ENST00000360581", "GRM1_S458A": "ENST00000282753",
              "GRM1_R889W": "ENST00000492807", "SIK3_N783Y": "ENST00000445177"}
ANC = {"afr": "African/African American", "amr": "Admixed American", "asj": "Ashkenazi Jewish",
       "eas": "East Asian", "fin": "Finnish", "mid": "Middle Eastern", "nfe": "Non-Finnish European",
       "sas": "South Asian", "remaining": "Remaining"}


def parse_hgvsp(h):
    m = re.match(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$", h or "")
    if not m:
        return None
    return three.get(m.group(1)), int(m.group(2)), three.get(m.group(3))


rows, anc_rows, codon_rows = [], [], []
for v in variants:
    g = json.load(open(os.path.join(L.DATA, "gnomad", f"gnomad_v4.1_{v['gene']}.json")))["data"]["gene"]
    tx = TRANSCRIPT[v["variant_id"]]
    hit = None
    for x in g["variants"]:
        p = parse_hgvsp(x["hgvsp"])
        if x["transcript_id"] != tx or p is None:
            continue
        ref, pos, alt = p
        if pos != v["pos"] or ref != v["ref_aa"]:
            continue
        ex, ge = x["exome"] or {}, x["genome"] or {}
        ac = (ex.get("ac") or 0) + (ge.get("ac") or 0)
        an = (ex.get("an") or 0) + (ge.get("an") or 0)
        hom = (ex.get("ac_hom") or 0) + (ge.get("ac_hom") or 0)
        rec = {"variant_id": v["variant_id"], "gnomad_variant": x["variant_id"], "rsid": ",".join(x["rsids"] or []),
               "hgvsp": x["hgvsp"], "transcript": tx, "alt_aa_observed": alt, "ac": ac, "an": an, "hom": hom,
               "ac_exome": ex.get("ac"), "an_exome": ex.get("an"), "ac_genome": ge.get("ac"), "an_genome": ge.get("an"),
               "filters_exome": ",".join(ex.get("filters") or []), "filters_genome": ",".join(ge.get("filters") or []),
               "af": ac / an if an else 0}
        if alt == v["alt_aa"]:
            hit = rec
            pops = defaultdict(lambda: [0, 0, 0])
            for src in (ex, ge):
                for pp in src.get("populations") or []:
                    if pp["id"] in ANC:
                        pops[pp["id"]][0] += pp["ac"]
                        pops[pp["id"]][1] += pp["an"]
                        pops[pp["id"]][2] += pp.get("ac_hom") or 0
            for pid, (pac, pan, phom) in sorted(pops.items()):
                anc_rows.append({"variant_id": v["variant_id"], "ancestry": pid, "ancestry_label": ANC[pid],
                                 "ac": pac, "an": pan, "hom": phom,
                                 "af": f"{pac / pan:.3g}" if pan else ""})
        else:
            codon_rows.append(dict(rec, reported_variant=v["variant_id"]))
    if hit is None:
        rows.append({"variant_id": v["variant_id"], "ac": 0, "an": "", "af": 0, "note": "not in gnomAD v4.1"})
        continue
    af = hit["af"]
    hit["af"] = f"{af:.3g}"
    hit["carrier_frequency"] = f"{2 * af:.3g}"
    hit["one_in_n_carriers"] = int(round(1 / (2 * af))) if af > 0 else ""
    hit["note"] = ""
    rows.append(hit)

L.write_tsv(os.path.join(L.TAB, "gnomad_fnss_variants.tsv"), rows,
            ["variant_id", "gnomad_variant", "rsid", "hgvsp", "transcript", "ac", "an", "af", "hom",
             "carrier_frequency", "one_in_n_carriers", "ac_exome", "an_exome", "ac_genome", "an_genome",
             "filters_exome", "filters_genome", "note"])
L.write_tsv(os.path.join(L.TAB, "gnomad_fnss_by_ancestry.tsv"), anc_rows)
L.write_tsv(os.path.join(L.TAB, "gnomad_same_codon.tsv"), codon_rows,
            ["reported_variant", "gnomad_variant", "rsid", "hgvsp", "transcript", "alt_aa_observed", "ac", "an", "hom",
             "filters_exome", "filters_genome"])

# ---------------------------------------------------------------- AlphaMissense
am = defaultdict(dict)
with open(os.path.join(L.DATA, "alphamissense", "AlphaMissense_aa_substitutions.fnss.tsv")) as fh:
    for line in fh:
        if line.startswith("#") or line.startswith("uniprot_id"):
            continue
        u, pv, s, c = line.rstrip("\n").split("\t")
        am[u][pv] = (float(s), c)
amiso = defaultdict(dict)
isof = os.path.join(L.DATA, "alphamissense", "AlphaMissense_isoforms_aa_substitutions.fnss.tsv")
if os.path.exists(isof):
    with open(isof) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("transcript_id"):
                continue
            f = line.rstrip("\n").split("\t")
            amiso[f[0].split(".")[0]][f[1]] = (float(f[2]), f[3] if len(f) > 3 else "")

species = L.read_tsv(os.path.join(L.TAB, "species_site_calls.tsv"))
other_res = defaultdict(set)
for r in species:
    if r["call"] and r["evidence"] != "conflict" and r["call"] != r["ref_aa"] and r["species"] != "Homo sapiens":
        other_res[r["variant_id"]].add(r["call"])

am_rows = []
# the isoform that carries each variant's numbering: UniProt canonical, except GRM1 R889W (mGluR1b,
# ENST00000507907) and SIK3 N783Y (MANE isoform ENST00000445177; N735 of the UniProt canonical isoform)
ISO = {"GRM1_R889W": "ENST00000507907", "SIK3_N783Y": "ENST00000445177"}
CLS = {"likely_benign": "benign", "likely_pathogenic": "pathogenic", "ambiguous": "ambiguous",
       "benign": "benign", "pathogenic": "pathogenic"}
for v in variants:
    canon_score = ""
    if v["uniprot"] in am and f"{v['ref_aa']}{v['uniprot_pos']}{v['alt_aa']}" in am[v["uniprot"]]:
        canon_score = am[v["uniprot"]][f"{v['ref_aa']}{v['uniprot_pos']}{v['alt_aa']}"][0]
    if v["variant_id"] in ISO:
        table = amiso.get(ISO[v["variant_id"]], {})
        src, upos = f"AlphaMissense isoforms ({ISO[v['variant_id']]})", v["pos"]
    else:
        table = am.get(v["uniprot"], {})
        src, upos = f"AlphaMissense ({v['uniprot']})", int(v["uniprot_pos"])
    key = f"{v['ref_aa']}{upos}{v['alt_aa']}"
    if key not in table:
        am_rows.append({"variant_id": v["variant_id"], "source": src, "substitution": key, "note": "not found"})
        continue
    s, c = table[key]
    c = CLS.get(c, c)
    allscores = [x[0] for x in table.values()]
    pct = 100.0 * sum(1 for x in allscores if x < s) / len(allscores)
    site = {k: x for k, x in table.items() if re.match(rf"^{v['ref_aa']}{upos}[A-Z]$", k)}
    others = sorted(other_res[v["variant_id"]])
    am_rows.append({
        "variant_id": v["variant_id"], "source": src, "substitution": key, "am_score": s, "am_class": c,
        "am_score_uniprot_canonical": canon_score,
        "percentile_in_protein": round(pct, 1),
        "site_mean_score": round(sum(x[0] for x in site.values()) / len(site), 3),
        "site_n_pathogenic": sum(1 for x in site.values() if CLS.get(x[1], x[1]) == "pathogenic"),
        "residues_in_other_mammals": ",".join(others),
        "scores_of_residues_in_other_mammals": ";".join(
            f"{a}:{site.get(v['ref_aa'] + str(upos) + a, ('',))[0]}" for a in others),
        "note": "",
    })
L.write_tsv(os.path.join(L.TAB, "alphamissense_fnss.tsv"), am_rows)

for r in rows:
    print(r.get("variant_id"), r.get("rsid"), "AC", r.get("ac"), "AN", r.get("an"), "AF", r.get("af"),
          "hom", r.get("hom"), "1 in", r.get("one_in_n_carriers"), r.get("note"))
for r in am_rows:
    print(r)
