"""11_scan_sra.py OUT [FIRST LAST] – the reads of the Brazilian porcupine (Coendou prehensilis; DNA Zoo whole-genome
sequencing, SRA run SRR10331604, 463,919,562 read pairs) that contain the 20 bases on either side of NPSR1 codon 206,
on either strand.

  OUT/reads/coendou_SRR10331604_NPSR1_reads.txt   one read pair per line (both reads, 151 + 151 bases, joined),
                                                  as vdb-dump writes them

The whole run is read once from NCBI's servers with vdb-dump (SRA Toolkit; set VDB_DUMP to its path if it is not
on the PATH), in chunks of 20 million pairs, six at a time, and each chunk is passed through grep. vdb-dump keeps
a cache of what it has read, so the chunks are kept small enough for the disk. A chunk is done when
OUT/reads/scan/FIRST-LAST.done exists: run the script again to finish an interrupted scan. It takes many hours.
FIRST and LAST (rows of the run) limit the scan, to try it out.
"""
import os
import subprocess
import sys
import time

RUN = "SRR10331604"
SPOTS = 463919562
FLANKS = ["CTCTGTGGCCTGATGACTCC", "TGGACTCCATACATGACCGT"]  # the 20 bases before and after codon 206 (coding strand)
CHUNK = 20_000_000
NWORK = 6
VDB = os.environ.get("VDB_DUMP", "vdb-dump")

rc = lambda s: s[::-1].translate(str.maketrans("ACGT", "TGCA"))
pats = []
for f in FLANKS:
    pats += ["-e", f, "-e", rc(f)]
out = os.path.join(sys.argv[1], "reads")
work = os.path.join(out, "scan")
os.makedirs(work, exist_ok=True)
first, last = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (1, SPOTS)
chunks = [(s, min(last, s + CHUNK - 1)) for s in range(first, last + 1, CHUNK)]
todo = [c for c in chunks if not os.path.exists(os.path.join(work, f"{c[0]}-{c[1]}.done"))]
running = {}
while todo or running:
    while todo and len(running) < NWORK:
        s, e = todo.pop(0)
        base = os.path.join(work, f"{s}-{e}")
        fo, fe = open(base + ".txt", "w"), open(base + ".err", "w")
        p1 = subprocess.Popen([VDB, "-C", "READ", "-f", "tab", "-R", f"{s}-{e}", RUN], stdout=subprocess.PIPE, stderr=fe)
        p2 = subprocess.Popen(["grep", "-F"] + pats, stdin=p1.stdout, stdout=fo)
        p1.stdout.close()
        running[(s, e)] = (p1, p2, fo, fe, time.time())
    time.sleep(10)
    for key in list(running):
        p1, p2, fo, fe, t0 = running[key]
        if p1.poll() is not None and p2.poll() is not None:
            fo.close()
            fe.close()
            if p1.returncode == 0:
                open(os.path.join(work, f"{key[0]}-{key[1]}.done"), "w").write(f"{time.time() - t0:.0f}s\n")
            else:
                todo.append(key)  # a dropped connection: read the chunk again
            del running[key]
            print("chunk", key, "exit", p1.returncode, f"{time.time() - t0:.0f}s;", len(todo), "to go", flush=True)
with open(os.path.join(out, f"coendou_{RUN}_NPSR1_reads.txt"), "w") as fh:
    for s, e in chunks:
        fh.write(open(os.path.join(work, f"{s}-{e}.txt")).read())
print("reads:", sum(1 for _ in open(os.path.join(out, f"coendou_{RUN}_NPSR1_reads.txt"))))
