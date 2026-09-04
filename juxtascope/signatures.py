"""
Two modes:
  derive_signatures(group_key)         -> each group vs the REST,
                                          keep genes with logFC > min_lfc,
                                          top-N by logFC (specificity filter).
  derive_sibling_signatures(comp_key,  -> within  each compartment,
                            group_key)     rank each group vs its SIBLINGS only,
                                           so shared pan-lineage genes (Ptprc,
                                           Vim, ...) are stripped out.

The logFC filter is what keeps broad genes (Vim, Cd74, Col4a1) OUT of the
signatures — without it every cell scores high on everything.
"""
import numpy as np
import pandas as pd
import scanpy as sc


def derive_signatures(adata, group_key, n_genes=25, min_lfc=1.0,
                      method="wilcoxon", exclude=None, verbose=True):
    """each group vs rest, logFC>min_lfc, top-N by logFC."""
    exclude = exclude or {}
    if verbose:
        print(f"  [sig] {group_key}: {adata.obs[group_key].nunique()} groups, "
              f"logFC>{min_lfc}")
    tmp = adata
    tmp.obs[group_key] = tmp.obs[group_key].astype("category")
    vc = tmp.obs[group_key].value_counts()
    keep = vc[vc >= 2].index.tolist()
    tmp = tmp[tmp.obs[group_key].isin(keep)].copy()
    sc.tl.rank_genes_groups(tmp, group_key, method=method,
                            key_added="_rk", n_genes=80)
    sigs = {}
    for g in tmp.obs[group_key].cat.categories:
        d = sc.get.rank_genes_groups_df(tmp, group=g, key="_rk")
        d = d[d["logfoldchanges"] > min_lfc].sort_values(
            "logfoldchanges", ascending=False)
        ex = exclude.get(g, set())
        genes = [x for x in d["names"].tolist() if x not in ex][:n_genes]
        sigs[g] = genes
    return sigs


def derive_sibling_signatures(adata, comp_key, group_key, n_genes=20,
                              min_lfc=0.5, min_per_group=100, exclude=None,
                              method="wilcoxon", verbose=True):
    """within each compartment, rank each group vs its SIBLINGS.
    Returns {compartment: {group: [genes]}}. Strips pan-lineage shared genes."""
    exclude = exclude or {}
    out = {}
    for comp in adata.obs[comp_key].astype("category").cat.categories:
        sub = adata[adata.obs[comp_key] == comp].copy()
        sub.obs[group_key] = sub.obs[group_key].astype("category")
        vc = sub.obs[group_key].value_counts()
        keep = vc[vc >= min_per_group].index.tolist()
        if len(keep) < 2:
            continue  # single-group compartment: no within-doublets possible
        sub = sub[sub.obs[group_key].isin(keep)].copy()
        sub.obs[group_key] = sub.obs[group_key].cat.remove_unused_categories()
        sc.tl.rank_genes_groups(sub, group_key, method=method,
                                key_added="_rk", n_genes=80)
        comp_sigs = {}
        for g in sub.obs[group_key].cat.categories:
            d = sc.get.rank_genes_groups_df(sub, group=g, key="_rk")
            d = d[d["logfoldchanges"] > min_lfc].sort_values(
                "logfoldchanges", ascending=False)
            ex = exclude.get(g, set())
            comp_sigs[g] = [x for x in d["names"].tolist() if x not in ex][:n_genes]
        out[comp] = comp_sigs
        if verbose:
            print(f"  [sib] {comp}: {list(comp_sigs.keys())}")
    return out


def score_signatures(adata, signatures, prefix="sig_", verbose=False):
    cols = []
    for g, genes in signatures.items():
        present = [x for x in genes if x in adata.var_names]
        if len(present) < 3:
            continue
        name = f"{prefix}{g}"
        sc.tl.score_genes(adata, present, score_name=name, use_raw=False)
        cols.append(name)
    return cols
