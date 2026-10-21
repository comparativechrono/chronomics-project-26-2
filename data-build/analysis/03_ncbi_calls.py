"""Step 03: the residue at each FNSS position in the annotated orthologs of NCBI (NCBI Gene orthologs,
RefSeq proteins of all isoforms; vertebrates).

For each species and site, every protein isoform is aligned to the human isoform that carries the
variant, and the isoform with the best-aligned window around the site is kept (for GRM1 R889W only
isoforms with the mGluR1b C-terminus align there). The residues of all isoforms are listed too.

Output: outputs/tables/ncbi_site_calls.tsv (one row per species x site)
"""
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fnsslib as L

variants = L.load_variants()
hp = L.human_proteins()
tax = {r["tax_id"]: r for r in json.load(open(os.path.join(L.DATA, "ncbi", "taxonomy.json")))}

rows = []
for gid in sorted({v["gene_id"] for v in variants}, key=int):
    gvars = [v for v in variants if v["gene_id"] == gid]
    gene = gvars[0]["gene"]
    orth = json.load(open(os.path.join(L.DATA, "ncbi", f"orthologs_{gid}.json")))["reports"]
    ginfo = {r["gene"]["gene_id"]: r["gene"] for r in orth}
    fa = L.readfa(os.path.join(L.DATA, "ncbi", f"proteins_{gid}.faa"))
    by_gene = {}
    for acc, seq in fa.items():
        head = fa.headers[acc]
        m = re.search(r"\[GeneID=(\d+)\]", head)
        iso = re.search(r"\[isoform=([^\]]+)\]", head)
        by_gene.setdefault(m.group(1), []).append((acc, seq, iso.group(1) if iso else ""))
    # one alignment per (protein, human reference isoform)
    refs = sorted({v["refseq_protein"] for v in gvars})
    for g, prots in by_gene.items():
        gi = ginfo.get(g)
        if gi is None:
            continue
        t = tax.get(gi["tax_id"], {})
        cache = {}
        for ref in refs:
            human = hp[ref]
            for acc, seq, iso in prots:
                cache[(ref, acc)] = L.map_positions(human, seq)
        for v in gvars:
            human = hp[v["refseq_protein"]]
            cands = []
            for acc, seq, iso in prots:
                m, ins, score = cache[(v["refseq_protein"], acc)]
                c = L.call_site(human, seq, m, ins, v["pos"])
                c.update(accession=acc, isoform=iso, length=len(seq))
                cands.append(c)
            # best isoform: confident first, then window identity, then global identity, NP_ before XP_
            cands.sort(key=lambda c: (L.confident(c), c["win_id"], -c["win_gaps"], c["global_id"],
                                      c["accession"].startswith("NP_"), c["length"]), reverse=True)
            b = cands[0]
            all_res = ";".join(f"{c['accession']}:{c['aa']}" for c in cands)
            conf_res = sorted({c["aa"] for c in cands if L.confident(c)})
            rows.append({
                "variant_id": v["variant_id"], "gene": gene, "source": "NCBI",
                "ncbi_gene_id": g, "species": gi["taxname"], "common_name": gi.get("common_name", ""),
                "taxid": gi["tax_id"], "class": t.get("cls", ""), "order": t.get("order", ""),
                "family": t.get("family", ""), "lineage": t.get("lineage", ""),
                "n_isoforms": len(prots), "accession": b["accession"], "isoform": b["isoform"],
                "ref_aa": v["ref_aa"], "alt_aa": v["alt_aa"], "aa": b["aa"], "position_in_protein": b["other_index"],
                "win_id": b["win_id"], "win_gaps": b["win_gaps"], "win_ins": b["win_ins"],
                "global_id": b["global_id"], "coverage": b["coverage"],
                "context_human": b["context_human"], "context_other": b["context_other"],
                "confident": int(L.confident(b)), "confident_residues_all_isoforms": ",".join(conf_res),
                "all_isoforms": all_res,
            })
    print(gene, "species", len(by_gene))

L.write_tsv(os.path.join(L.TAB, "ncbi_site_calls.tsv"), rows)
print("rows", len(rows))
