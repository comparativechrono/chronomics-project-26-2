"""08_sleep.py OUT – the two compilations of mammal sleep.

  OUT/sleep/Complete_Sleep_Dataset_FINAL.xml  the Phylogeny of Sleep database (McNamara et al. 2008), an Excel 2003
                                              workbook (XML) from the project's web site, Boston University
  OUT/sleep/phylogeny_of_sleep_db.tsv         the same, as a table: one record (a study of a species) per row
  OUT/sleep/msleep.tsv                        msleep (Savage & West 2007), as shipped with the R package ggplot2
                                              (needs R with ggplot2; written by 08_msleep.R)
"""
import csv
import os
import subprocess
import sys
import urllib.request
import xml.etree.ElementTree as ET

OUT = os.path.join(sys.argv[1], "sleep")
os.makedirs(OUT, exist_ok=True)
URL = "https://www.bu.edu/phylogeny/datasets/Complete_Sleep_Dataset_FINAL.xml"
xml = os.path.join(OUT, "Complete_Sleep_Dataset_FINAL.xml")
req = urllib.request.Request(URL, headers={"User-Agent": "Mozilla/5.0"})
open(xml, "wb").write(urllib.request.urlopen(req, timeout=300).read())

# the worksheet's rows; a cell with ss:Index starts at that column (the cells before it are empty)
SS = "{urn:schemas-microsoft-com:office:spreadsheet}"
rows = []
for row in ET.parse(xml).getroot().iter(SS + "Row"):
    cells = []
    for cell in row.findall(SS + "Cell"):
        i = cell.get(SS + "Index")
        if i:
            cells += [""] * (int(i) - 1 - len(cells))
        d = cell.find(SS + "Data")
        # line breaks and tabs in a cell become spaces
        cells.append((d.text or "").replace("\n", " ").replace("\t", " ") if d is not None else "")
    rows.append(cells)
n = len(rows[0])
with open(os.path.join(OUT, "phylogeny_of_sleep_db.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    for r in rows:
        w.writerow((r + [""] * n)[:n])
print(len(rows) - 1, "records of the Phylogeny of Sleep database")

here = os.path.dirname(os.path.abspath(__file__))
subprocess.run(["Rscript", os.path.join(here, "08_msleep.R"), os.path.join(OUT, "msleep.tsv")], check=True)
