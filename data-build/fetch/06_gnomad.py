"""06_gnomad.py OUT – every variant of the five genes in gnomAD v4.1 (exomes and genomes, GRCh38), with allele
counts overall and by genetic ancestry group, and each gene's constraint, from gnomAD's GraphQL API.

  OUT/gnomad/gnomad_v4.1_GENE.json   (GENE = BHLHE41, ADRB1, NPSR1, GRM1, SIK3) the API's answer, as it came
"""
import json
import os
import sys
import time
import urllib.request

OUT = os.path.join(sys.argv[1], "gnomad")
os.makedirs(OUT, exist_ok=True)
Q = """query($g:String!){ gene(gene_symbol:$g, reference_genome: GRCh38) { gene_id symbol canonical_transcript_id mane_select_transcript { ensembl_id refseq_id }
  gnomad_constraint { exp_lof obs_lof oe_lof oe_lof_upper exp_mis obs_mis oe_mis oe_mis_upper mis_z syn_z pLI }
  variants(dataset: gnomad_r4) { variant_id pos ref alt rsids consequence hgvsc hgvsp transcript_id
    exome { ac an ac_hom faf95 { popmax popmax_population } populations { id ac an ac_hom } filters }
    genome { ac an ac_hom faf95 { popmax popmax_population } populations { id ac an ac_hom } filters } } } }"""


def gql(query, variables, tries=5):
    body = json.dumps({"query": query, "variables": variables}).encode()
    for k in range(tries):
        try:
            req = urllib.request.Request("https://gnomad.broadinstitute.org/api", data=body,
                                         headers={"Content-Type": "application/json", "User-Agent": "Mozilla/5.0"})
            return json.loads(urllib.request.urlopen(req, timeout=300).read())
        except Exception as e:  # the API limits how often it may be asked: wait and retry
            if k == tries - 1:
                raise
            print("  retrying after", e, file=sys.stderr)
            time.sleep(30 * (k + 1))


for g in ["BHLHE41", "ADRB1", "NPSR1", "GRM1", "SIK3"]:
    d = gql(Q, {"g": g})
    if "errors" in d:
        sys.exit(f"{g}: {d['errors'][:2]}")
    json.dump(d, open(os.path.join(OUT, f"gnomad_v4.1_{g}.json"), "w"))
    print(g, len(d["data"]["gene"]["variants"]), "variants")
    time.sleep(10)
