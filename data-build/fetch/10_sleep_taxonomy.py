"""10_sleep_taxonomy.py OUT SPECIES_TABLE – the NCBI taxonomy of each sleep species: its taxon, order, family and
common name, which the student table gives every species (msleep's orders are older ones, and some species
have no common name there).

SPECIES_TABLE is a tab-separated table with a column 'species' (binomials, "Genus species"): the analysis
makes it, outputs/tables/sleep_species.tsv of 07_sleep_data.py; make_bundle.sh runs this script after that step
if the output is missing.
  OUT/sleep/sleep_species_taxonomy.json   {species: {taxid, name, order, family, common}}, in the table's order
One esearch per species (by scientific name, then by any name), slowly, as NCBI asks of users without an API
key, then one efetch for all.
"""
import csv
import json
import os
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def get(url):
    for k in range(6):
        try:
            return urllib.request.urlopen(url, timeout=60).read()
        except Exception:
            time.sleep(2 + 3 * k)
    raise SystemExit("failed " + url)


species = [r["species"] for r in csv.DictReader(open(sys.argv[2]), delimiter="\t")]
taxid = {}
for sp in species:
    ids = []
    for term in (f'"{sp}"[Scientific Name]', sp):
        d = json.loads(get(f"{E}/esearch.fcgi?db=taxonomy&term={urllib.parse.quote(term)}&retmode=json"))
        ids = d["esearchresult"]["idlist"]
        if ids:
            break
    taxid[sp] = ids[0] if ids else ""
    time.sleep(0.5)
x = get(f"{E}/efetch.fcgi?db=taxonomy&id=" + ",".join(v for v in taxid.values() if v) + "&retmode=xml")
info = {}
for t in ET.fromstring(x).findall("Taxon"):
    ranks = {}
    lex = t.find("LineageEx")
    for r in (lex if lex is not None else []):
        if r.findtext("Rank") in ("order", "family"):
            ranks[r.findtext("Rank")] = r.findtext("ScientificName")
    other = t.find("OtherNames")
    common = ""
    if other is not None:
        common = other.findtext("GenbankCommonName") or other.findtext("CommonName") or ""
    info[t.findtext("TaxId")] = dict(name=t.findtext("ScientificName"), order=ranks.get("order", ""),
                                     family=ranks.get("family", ""), common=common)
res = {sp: dict(taxid=tid, **info.get(tid, {})) for sp, tid in taxid.items()}
os.makedirs(os.path.join(sys.argv[1], "sleep"), exist_ok=True)
json.dump(res, open(os.path.join(sys.argv[1], "sleep", "sleep_species_taxonomy.json"), "w"), indent=1)
print(len(res), "species;", "without taxonomy:", [sp for sp, v in res.items() if not v.get("order")])
