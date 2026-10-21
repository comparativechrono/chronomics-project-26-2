"""Step 04: one residue per mammal species and FNSS site, from the two independent sources.

TOGA (genome alignments; several assemblies per species possible) and NCBI (annotated proteins) are
merged by species (binomial). A species' call is
  'both'   - both sources have a confident call and agree,
  'TOGA' / 'NCBI' - only one source has a confident call,
  'conflict' - the sources (or two assemblies of the species) disagree.
Also writes the vertebrate (non-mammal) calls of NCBI for context.

Outputs: outputs/tables/species_site_calls.tsv (mammals), outputs/tables/vertebrate_site_calls.tsv
"""
import json
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fnsslib as L

tax = {r["tax_id"]: r for r in json.load(open(os.path.join(L.DATA, "ncbi", "taxonomy.json")))}
variants = {v["variant_id"]: v for v in L.load_variants()}


def binomial(name):
    w = name.replace("[", "").replace("]", "").split()
    return " ".join(w[:2]) if len(w) >= 2 else name


toga = L.read_tsv(os.path.join(L.TAB, "toga_site_calls.tsv"))
ncbi = L.read_tsv(os.path.join(L.TAB, "ncbi_site_calls.tsv"))

rec = defaultdict(lambda: {"toga": [], "ncbi": [], "info": {}})
for r in toga:
    t = tax.get(r["taxid"], {})
    sp = binomial(t.get("name") or r["species"])
    k = (r["variant_id"], sp)
    rec[k]["toga"].append(r)
    rec[k]["info"].update(order=t.get("order", ""), family=t.get("family", ""), genus=t.get("genus", ""),
                          common=rec[k]["info"].get("common") or r["common_name"], taxid=r["taxid"])
for r in ncbi:
    if r["class"] != "Mammalia":
        continue
    sp = binomial(r["species"])
    k = (r["variant_id"], sp)
    rec[k]["ncbi"].append(r)
    if not rec[k]["info"].get("order"):
        t = tax.get(r["taxid"], {})
        rec[k]["info"].update(order=r["order"] or t.get("order", ""), family=r["family"] or t.get("family", ""),
                              genus=t.get("genus", ""), taxid=r["taxid"])
    rec[k]["info"].setdefault("common", r["common_name"])

rows = []
for (vid, sp), d in sorted(rec.items()):
    v = variants[vid]
    tconf = [r for r in d["toga"] if r["confident"] == "1"]
    nconf = [r for r in d["ncbi"] if r["confident"] == "1"]
    tres = sorted({r["aa"] for r in tconf})
    nres = sorted({r["aa"] for r in nconf})
    allres = sorted(set(tres) | set(nres))
    if len(allres) == 0:
        call, evidence = "", "none"
    elif len(allres) > 1:
        call, evidence = "/".join(allres), "conflict"
    else:
        call = allres[0]
        evidence = "both" if (tres and nres) else ("TOGA" if tres else "NCBI")
    codons = sorted({r["codon"] for r in tconf})
    rows.append({
        "variant_id": vid, "gene": v["gene"], "species": sp, "common_name": d["info"].get("common", ""),
        "order": d["info"].get("order", ""), "family": d["info"].get("family", ""),
        "genus": d["info"].get("genus", "") or sp.split()[0], "taxid": d["info"].get("taxid", ""),
        "ref_aa": v["ref_aa"], "alt_aa": v["alt_aa"], "call": call, "evidence": evidence,
        "is_human_residue": int(call == v["ref_aa"]), "is_fnss_residue": int(call == v["alt_aa"]),
        "toga_residues": ",".join(tres), "toga_codons": ",".join(codons),
        "toga_assemblies": ",".join(sorted({r["assembly"] for r in d["toga"]})),
        "toga_n_confident": len(tconf), "toga_n": len(d["toga"]),
        "ncbi_residue": ",".join(nres), "ncbi_accession": ",".join(r["accession"] for r in d["ncbi"]),
        "toga_context": (tconf[0]["context_other"] if tconf else (d["toga"][0]["context_other"] if d["toga"] else "")),
        "ncbi_context": (nconf[0]["context_other"] if nconf else (d["ncbi"][0]["context_other"] if d["ncbi"] else "")),
        "min_win_id": min([float(r["win_id"]) for r in tconf + nconf], default=""),
    })
L.write_tsv(os.path.join(L.TAB, "species_site_calls.tsv"), rows)

# non-mammal vertebrates (NCBI only)
vrows = []
for r in ncbi:
    if r["class"] == "Mammalia":
        continue
    v = variants[r["variant_id"]]
    vrows.append({
        "variant_id": r["variant_id"], "gene": r["gene"], "species": binomial(r["species"]),
        "common_name": r["common_name"], "class": r["class"] or "other", "order": r["order"],
        "family": r["family"], "ref_aa": v["ref_aa"], "alt_aa": v["alt_aa"],
        "call": r["aa"] if r["confident"] == "1" else "", "accession": r["accession"],
        "context": r["context_other"], "win_id": r["win_id"],
    })
L.write_tsv(os.path.join(L.TAB, "vertebrate_site_calls.tsv"), vrows)

# summary
from collections import Counter
for vid in variants:
    rr = [r for r in rows if r["variant_id"] == vid and r["call"]]
    c = Counter(r["call"] for r in rr)
    e = Counter(r["evidence"] for r in rr)
    print(f"{vid:14s} species with a call {len(rr):4d}  {dict(c.most_common())}  evidence {dict(e)}")
