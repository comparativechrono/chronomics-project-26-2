"""02_ncbi.py OUT – the orthologs of the five genes in NCBI (NCBI Datasets v2), their proteins, and the taxonomy
of every species in them and in TOGA's table of assemblies (NCBI E-utilities).

  OUT/ncbi/orthologs_GENEID.json  the ortholog report of each human gene (79365 BHLHE41, 153 ADRB1, 387129 NPSR1,
                                  2911 GRM1, 23387 SIK3): one gene per species
  OUT/ncbi/proteins_GENEID.faa    the RefSeq proteins of all those genes, every isoform
  OUT/ncbi/taxonomy.json          name, rank, common name, lineage, class, order, family and genus of each taxon:
                                  first those of the orthologs, then those of TOGA's assemblies (needs
                                  OUT/toga/overview.table.tsv from 01_toga.sh)
NCBI updates these records from time to time; a later download may differ a little.
"""
import csv
import io
import json
import os
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
import zipfile

OUT = os.path.join(sys.argv[1], "ncbi")
os.makedirs(OUT, exist_ok=True)
API = "https://api.ncbi.nlm.nih.gov/datasets/v2"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GENES = ["79365", "153", "387129", "2911", "23387"]


def fetch(url, data=None, headers=None, tries=6):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers=headers or {})
            return urllib.request.urlopen(req, timeout=600).read()
        except Exception as e:  # rate limits (429) and the odd time-out: wait and retry
            if k == tries - 1:
                raise
            print("  retrying after", e, file=sys.stderr)
            time.sleep(3 + 5 * k)


tax_ids = set()
for g in GENES:
    # the orthologs, one page of up to 1,000 (each gene has fewer)
    raw = fetch(f"{API}/gene/id/{g}/orthologs?page_size=1000")
    d = json.loads(raw)
    assert not d.get("next_page_token"), "more than one page: merge the pages"
    open(os.path.join(OUT, f"orthologs_{g}.json"), "wb").write(raw)
    ids = [r["gene"]["gene_id"] for r in d["reports"]]
    tax_ids |= {r["gene"]["tax_id"] for r in d["reports"]}
    # their proteins: a data package with protein.faa
    body = json.dumps({"gene_ids": [int(i) for i in ids], "include_annotation_type": ["FASTA_PROTEIN"]}).encode()
    z = zipfile.ZipFile(io.BytesIO(fetch(f"{API}/gene/download", data=body,
                                         headers={"Content-Type": "application/json"})))
    open(os.path.join(OUT, f"proteins_{g}.faa"), "wb").write(z.read("ncbi_dataset/data/protein.faa"))
    print(g, len(ids), "orthologs")
    time.sleep(1)


def taxonomy(ids, batch=200):
    """efetch db=taxonomy, in batches: one record per id asked for. An id that NCBI has merged into another taxon
    comes back as that taxon; it is then kept twice, under the new id and under the old one (its AkaTaxIds),
    so that the old ids of TOGA's table can still be looked up."""
    recs = []
    for i in range(0, len(ids), batch):
        asked = ids[i:i + batch]
        body = urllib.parse.urlencode({"db": "taxonomy", "id": ",".join(asked), "retmode": "xml"}).encode()
        root = ET.fromstring(fetch(f"{EUTILS}/efetch.fcgi", data=body))
        for t in root.findall("Taxon"):
            lex = t.find("LineageEx")
            ranks = {x.findtext("Rank"): x.findtext("ScientificName") for x in (lex if lex is not None else [])}
            other = t.find("OtherNames")
            common = ""
            if other is not None:
                common = other.findtext("GenbankCommonName") or other.findtext("CommonName") or ""
            rec = {"tax_id": t.findtext("TaxId"), "name": t.findtext("ScientificName"),
                   "rank": t.findtext("Rank"), "common": common, "lineage": t.findtext("Lineage"),
                   "cls": ranks.get("class", ""), "order": ranks.get("order", ""),
                   "family": ranks.get("family", ""), "genus": ranks.get("genus", "")}
            recs.append(rec)
            for old in t.findall("AkaTaxIds/TaxId"):
                if old.text in asked:
                    recs.append(dict(rec, tax_id=old.text, merged_into=rec["tax_id"]))
        time.sleep(0.5)
    return recs


first = sorted(tax_ids, key=int)
toga = [r["Species Taxonomy ID"] for r in
        csv.DictReader(open(os.path.join(sys.argv[1], "toga", "overview.table.tsv")), delimiter="\t")]
second = sorted(set(toga) - set(first), key=int)
recs = taxonomy(first) + taxonomy(second)
json.dump(recs, open(os.path.join(OUT, "taxonomy.json"), "w"))
print(len(recs), "taxonomy records")
