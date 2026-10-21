"""07_genebass.py OUT [GENE ...] – results for the five genes in Genebass (Karczewski et al. 2022): the exomes of
394,841 people of UK Biobank tested against their traits; here six sleep questions. From Genebass's API, as it came.

  OUT/genebass/phenotypes.json                    the list of traits
  OUT/genebass/variants_GENE_PHENOTYPE.json       every variant of the gene tested against the trait (single variants)
  OUT/genebass/gene_GENE.json                     the gene
  OUT/genebass/geneanalysis_GENE_PHENOTYPE.json   the gene-level (burden) tests, for sleep duration
Give gene names after OUT to fetch only those. Genebass sometimes answers 500 to a request it has not cached:
the script waits and asks again, and keeps files it already has, so running it again fills the gaps.
"""
import gzip
import json
import os
import sys
import time
import urllib.request

OUT = os.path.join(sys.argv[1], "genebass")
os.makedirs(OUT, exist_ok=True)
API = "https://main.genebass.org/api"
GENES = {"BHLHE41": "ENSG00000123095", "ADRB1": "ENSG00000043591", "NPSR1": "ENSG00000187258",
         "GRM1": "ENSG00000152822", "SIK3": "ENSG00000160584"}
PHENO = ["continuous-sleep_duration_1160-both_sexes--custom", "continuous-chronotype_1180_quant-both_sexes--custom",
         "continuous-easy_wakeup_1170_quant-both_sexes--custom", "continuous-insomnia_1200_quant-both_sexes--custom",
         "categorical-daytime_nap_1190_binary-both_sexes--custom", "categorical-narcolepsy_1220_binary-both_sexes--custom"]


def get(path, tries=4):
    for k in range(tries):
        try:
            req = urllib.request.Request(API + path, headers={"User-Agent": "Mozilla/5.0", "Accept-Encoding": "gzip"})
            raw = urllib.request.urlopen(req, timeout=600).read()
            break
        except Exception:
            if k == tries - 1:
                raise
            time.sleep(20 * (k + 1))
    try:
        return json.loads(raw)
    except Exception:
        return json.loads(gzip.decompress(raw))


def save(name, path):
    f = os.path.join(OUT, name)
    if os.path.exists(f):
        return
    try:
        json.dump(get(path), open(f, "w"))
    except Exception as ex:
        print("failed", name, ex, file=sys.stderr)
    time.sleep(1)


save("phenotypes.json", "/phenotypes")
for g in sys.argv[2:] or GENES:
    e = GENES[g]
    for p in PHENO:
        save(f"variants_{g}_{p}.json", f"/gene/{e}/variant-associations/{p}")
    save(f"gene_{g}.json", f"/gene/{e}")
    save(f"geneanalysis_{g}_{PHENO[0]}.json", f"/gene/{e}/analysis/{PHENO[0]}")
    print(g, "done")
