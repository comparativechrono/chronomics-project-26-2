"""Shared code for the FNSS comparative analysis.

Reading sequences, translating codons, and calling the residue that each species carries at a human
residue position by pairwise alignment to the human reference protein (with quality measures for the
call).
"""
import csv
import gzip
import os
from collections import OrderedDict

from Bio import Align
from Bio.Align import substitution_matrices

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "outputs")
TAB = os.path.join(OUT, "tables")
WORK = os.path.join(OUT, "work")
for _d in (OUT, TAB, WORK):
    os.makedirs(_d, exist_ok=True)

# window used for the quality of a call: the human residues pos-WIN .. pos+WIN
WIN = 10
# a call is 'confident' if the residue is aligned (not a gap, not masked) and the window
# is well aligned: at least MIN_WIN_ID identical residues (fraction of 2*WIN+1) and at most
# MAX_WIN_GAPS window positions without an aligned residue
MIN_WIN_ID = 0.5
MAX_WIN_GAPS = 4

_bases = "TCAG"
_aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {}
_i = 0
for _a in _bases:
    for _b in _bases:
        for _c in _bases:
            CODE[_a + _b + _c] = _aas[_i]
            _i += 1


def translate_codon(c):
    c = c.upper()
    if c == "---":
        return "-"
    return CODE.get(c, "X")


def readfa(path):
    """FASTA to an ordered dict {first word of header: sequence}; headers kept in .headers."""
    seqs = OrderedDict()
    heads = {}
    op = gzip.open if path.endswith(".gz") else open
    name = None
    with op(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0] if " " in line or "\t" not in line else line[1:]
                heads[name] = line[1:]
                seqs[name] = []
            else:
                seqs[name].append(line.strip())
    out = OrderedDict((k, "".join(v)) for k, v in seqs.items())
    out.headers = heads
    return out


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, rows, fields=None):
    if fields is None:
        fields = list(rows[0].keys()) if rows else []
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def load_variants(include_comparison=True):
    rows = read_tsv(os.path.join(DATA, "human", "variants.tsv"))
    for r in rows:
        r["pos"] = int(r["pos"])
    if not include_comparison:
        rows = [r for r in rows if r["category"] == "FNSS"]
    return rows


def human_proteins():
    s = readfa(os.path.join(DATA, "human", "human_refseq.fasta"))
    return dict(s)


def make_aligner():
    al = Align.PairwiseAligner()
    al.mode = "global"
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score = -10.0
    al.extend_gap_score = -0.5
    # free end gaps on both sequences (fragments, different isoform ends)
    al.end_insertion_score = 0.0
    al.end_deletion_score = 0.0
    return al


_ALIGNER = None


def map_positions(human, other):
    """Align `other` to `human`; return a list m with m[i] = index in `other` aligned to human
    residue i (0-based), or None if human residue i faces a gap; and the number of residues of
    `other` inserted next to each human residue (insertions after human residue i)."""
    global _ALIGNER
    if _ALIGNER is None:
        _ALIGNER = make_aligner()
    # X (masked codons) are scored like any mismatch by BLOSUM62's X row
    aln = _ALIGNER.align(human, other)[0]
    m = [None] * len(human)
    ins_after = [0] * len(human)
    hb, ob = aln.aligned
    for (h0, h1), (o0, o1) in zip(hb, ob):
        for k in range(h1 - h0):
            m[h0 + k] = o0 + k
    # insertions: residues of other between consecutive aligned blocks
    for j in range(1, len(hb)):
        h_prev_end = hb[j - 1][1]
        o_prev_end = ob[j - 1][1]
        o_start = ob[j][0]
        if o_start > o_prev_end and h_prev_end - 1 >= 0:
            ins_after[h_prev_end - 1] += o_start - o_prev_end
    return m, ins_after, aln.score


def call_site(human, other, m, ins_after, pos, win=WIN):
    """Residue of `other` at human position `pos` (1-based) with quality measures."""
    i = pos - 1
    j = m[i]
    aa = other[j] if j is not None else "-"
    lo, hi = max(0, i - win), min(len(human) - 1, i + win)
    n = hi - lo + 1
    ident = 0
    gaps = 0
    for k in range(lo, hi + 1):
        jj = m[k]
        if jj is None or other[jj] == "X":
            gaps += 1
        elif other[jj] == human[k]:
            ident += 1
    ins = sum(ins_after[lo:hi])  # insertions strictly inside the window
    ctx_h = human[max(0, i - 7): i] + "[" + human[i] + "]" + human[i + 1: i + 8]
    ctx_o = []
    for k in range(max(0, i - 7), min(len(human), i + 8)):
        jj = m[k]
        c = other[jj] if jj is not None else "-"
        ctx_o.append("[" + c + "]" if k == i else c)
    # global identity / coverage over all human residues
    aligned = sum(1 for k in range(len(human)) if m[k] is not None)
    same = sum(1 for k in range(len(human)) if m[k] is not None and other[m[k]] == human[k])
    return {
        "aa": aa,
        "other_index": (j + 1) if j is not None else "",
        "win_id": round(ident / n, 3),
        "win_gaps": gaps,
        "win_ins": ins,
        "context_human": ctx_h,
        "context_other": "".join(ctx_o),
        "global_id": round(same / max(aligned, 1), 3),
        "coverage": round(aligned / len(human), 3),
    }


def confident(c):
    return (c["aa"] not in ("-", "X", "*")) and c["win_id"] >= MIN_WIN_ID and c["win_gaps"] <= MAX_WIN_GAPS \
        and c["win_ins"] <= MAX_WIN_GAPS


# ---------------------------------------------------------------- TOGA alignments
TOGA_FILES = {
    "BHLHE41": ("ENST00000242728.BHLHE41", "NP_110389.1"),
    "ADRB1": ("ENST00000369295.ADRB1", "NP_000675.1"),
    "NPSR1": ("ENST00000381539.NPSR1", "NP_001287864.1"),  # longest CDS = isoform G; same numbering to 340
    "GRM1": ("ENST00000361719.GRM1", "NP_001264993.1"),  # isoform alpha; R889 of beta not included
    "SIK3": ("ENST00000445177.SIK3", "NP_001353615.1"),
}


def toga_overview():
    rows = read_tsv(os.path.join(DATA, "toga", "overview.table.tsv"))
    return {r["Assembly name"]: r for r in rows}


def toga_sequences(gene):
    """{assembly: (protein, [codons])} from TOGA's multiple codon alignment; the reference too."""
    f = os.path.join(DATA, "toga", TOGA_FILES[gene][0] + ".fasta.gz")
    aln = readfa(f)
    out = OrderedDict()
    for name, s in aln.items():
        head = aln.headers[name]
        asm = "REFERENCE" if head.startswith("REFERENCE") else head.split("\t")[0].split()[0].replace("vs_", "")
        proj = head.split("\t")[1] if "\t" in head else ""
        prot, cods = [], []
        for k in range(0, len(s) - 2, 3):
            c = s[k:k + 3]
            if c == "---":
                continue
            a = translate_codon(c)
            prot.append(a)
            cods.append(c.upper())
        out.setdefault(asm, []).append(("".join(prot), cods, proj))
    return out
