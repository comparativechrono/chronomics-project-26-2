"""Step 02: the residue at each FNSS position in every mammal assembly of TOGA's multiple codon
alignments (Kirilenko et al. 2023; human hg38 reference).

Each assembly's sequence is re-aligned to the human reference protein (the alignment of a single
pair is more robust than the column of a 400-sequence alignment in low-complexity regions), and the
residue facing each human FNSS position is called with quality measures (fnsslib.call_site).

Output: outputs/tables/toga_site_calls.tsv (one row per assembly x projection x site)
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fnsslib as L

variants = L.load_variants()
hp = L.human_proteins()
ov = L.toga_overview()


def isoform_map(var, target_acc):
    """Map a variant position on its own isoform to the isoform used by TOGA; None if the residue
    and its +-10 neighbourhood are not shared by the two isoforms."""
    src = hp[var["refseq_protein"]]
    tgt = hp[target_acc]
    if src == tgt:
        return var["pos"]
    m, ins, _ = L.map_positions(src, tgt)
    c = L.call_site(src, tgt, m, ins, var["pos"])
    if c["aa"] == var["ref_aa"] and c["win_id"] == 1.0:
        return m[var["pos"] - 1] + 1
    return None


rows = []
for gene, (fname, ref_acc) in L.TOGA_FILES.items():
    human = hp[ref_acc]
    seqs = L.toga_sequences(gene)
    # the TOGA reference must translate to the human RefSeq protein over the FNSS positions
    tref = seqs["REFERENCE"][0][0]
    gvars = [v for v in variants if v["gene"] == gene]
    sites = []
    for v in gvars:
        p = isoform_map(v, ref_acc)
        if p is None:
            print(f"{v['variant_id']}: not represented in the TOGA alignment ({fname})")
            continue
        assert human[p - 1] == v["ref_aa"], (v, p)
        sites.append((v, p))
    mref, insref, _ = L.map_positions(human, tref)
    for v, p in sites:
        c = L.call_site(human, tref, mref, insref, p)
        assert c["aa"] == v["ref_aa"] and c["win_id"] == 1.0, (gene, v["variant_id"], c)
    n = 0
    for asm, projs in seqs.items():
        if asm == "REFERENCE":
            continue
        info = ov.get(asm, {})
        for prot, cods, proj in projs:
            m, ins, score = L.map_positions(human, prot)
            for v, p in sites:
                c = L.call_site(human, prot, m, ins, p)
                j = m[p - 1]
                codon = cods[j] if j is not None else ""
                rows.append({
                    "variant_id": v["variant_id"], "gene": gene, "source": "TOGA",
                    "assembly": asm, "projection": proj, "n_projections": len(projs),
                    "species": info.get("Species", ""), "common_name": info.get("Common name", ""),
                    "taxid": info.get("Species Taxonomy ID", ""),
                    "lineage": info.get("Taxonomic Lineage", ""),
                    "assembly_accession": info.get("NCBI accession / source", ""),
                    "contig_n50": info.get("contig N50 (bp)", ""),
                    "human_pos_toga_isoform": p, "ref_aa": v["ref_aa"], "alt_aa": v["alt_aa"],
                    "aa": c["aa"], "codon": codon, "win_id": c["win_id"], "win_gaps": c["win_gaps"],
                    "win_ins": c["win_ins"], "global_id": c["global_id"], "coverage": c["coverage"],
                    "context_human": c["context_human"], "context_other": c["context_other"],
                    "confident": int(L.confident(c)),
                })
            n += 1
    print(f"{gene}: {n} projections in {len(seqs) - 1} assemblies; sites {[v['variant_id'] for v, _ in sites]}")

fields = list(rows[0].keys())
L.write_tsv(os.path.join(L.TAB, "toga_site_calls.tsv"), rows, fields)
print("rows", len(rows))
