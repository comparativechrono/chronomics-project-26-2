"""Step 09: the FNSS variants in the UK Biobank exomes (Genebass; 394,841 people; SAIGE single-variant
and gene-burden tests; Karczewski et al. 2022).

For each FNSS variant and sleep phenotype: allele count, effect (BETA) and P, with an approximate 95%
confidence interval (standard error recovered from BETA and the two-sided P). Sleep duration is in
hours (UK Biobank field 1160, adjusted for smoking, alcohol, BMI and BMI squared); chronotype is on a
0-3 scale (0 = definitely a morning person); napping is a case-control trait (log odds).
Also: gene-level burden tests of the five genes for sleep duration.

Outputs: outputs/tables/ukb_exome_fnss_variants.tsv, outputs/tables/ukb_exome_gene_burden.tsv
"""
import json
import os
import sys
from statistics import NormalDist

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fnsslib as L

G = os.path.join(L.DATA, "genebass")
variants = L.load_variants()
gn = {r["variant_id"]: r for r in L.read_tsv(os.path.join(L.TAB, "gnomad_fnss_variants.tsv"))}
PHENO = {"continuous-sleep_duration_1160-both_sexes--custom": ("sleep duration", "hours"),
         "continuous-chronotype_1180_quant-both_sexes--custom": ("chronotype", "scale 0-3, higher = later"),
         "continuous-easy_wakeup_1170_quant-both_sexes--custom": ("ease of waking", "scale 0-3"),
         "continuous-insomnia_1200_quant-both_sexes--custom": ("insomnia", "scale 0-2"),
         "categorical-daytime_nap_1190_binary-both_sexes--custom": ("daytime napping", "log odds"),
         "categorical-narcolepsy_1220_binary-both_sexes--custom": ("daytime dozing", "log odds")}
nd = NormalDist()
rows = []
for v in variants:
    vid = gn.get(v["variant_id"], {}).get("gnomad_variant", "")
    for pid, (label, unit) in PHENO.items():
        f = os.path.join(G, f"variants_{v['gene']}_{pid}.json")
        if not os.path.exists(f):
            continue
        hit = [x for x in json.load(open(f)) if x["variant_id"] == vid]
        if not hit:
            rows.append({"variant_id": v["variant_id"], "genebass_variant": vid, "phenotype": label, "unit": unit,
                         "ac": 0, "note": "not tested in Genebass"})
            continue
        x = hit[0]
        b, p = x["BETA"], x["Pvalue"]
        z = nd.inv_cdf(1 - p / 2) if 0 < p < 1 else float("nan")
        se = abs(b) / z if z and z == z and z > 0 else float("nan")
        rows.append({"variant_id": v["variant_id"], "genebass_variant": vid, "phenotype": label, "unit": unit,
                     "ac": int(x["AC"]), "af": f"{x['AF']:.3g}", "beta": round(b, 4), "p": f"{p:.3g}",
                     "se_approx": round(se, 4) if se == se else "",
                     "ci95_low": round(b - 1.96 * se, 3) if se == se else "",
                     "ci95_high": round(b + 1.96 * se, 3) if se == se else "", "note": ""})
L.write_tsv(os.path.join(L.TAB, "ukb_exome_fnss_variants.tsv"), rows,
            ["variant_id", "genebass_variant", "phenotype", "unit", "ac", "af", "beta", "se_approx", "ci95_low",
             "ci95_high", "p", "note"])

grows = []
for g in sorted({v["gene"] for v in variants}):
    f = os.path.join(G, f"geneanalysis_{g}_continuous-sleep_duration_1160-both_sexes--custom.json")
    if not os.path.exists(f):
        continue
    for x in json.load(open(f)):
        grows.append({"gene": g, "phenotype": "sleep duration", "burden_set": x.get("annotation"),
                      "n_variants": x.get("total_variants"), "p_skato": x.get("Pvalue"),
                      "p_burden": x.get("Pvalue_Burden"), "p_skat": x.get("Pvalue_SKAT"),
                      "beta_burden_h": x.get("BETA_Burden")})
L.write_tsv(os.path.join(L.TAB, "ukb_exome_gene_burden.tsv"), grows)
for r in rows:
    if r["phenotype"] in ("sleep duration", "chronotype"):
        print(r["variant_id"], r["phenotype"], "AC", r.get("ac"), "beta", r.get("beta"), "95% CI",
              r.get("ci95_low"), r.get("ci95_high"), "P", r.get("p"), r.get("note"))
for r in grows:
    print(r)
