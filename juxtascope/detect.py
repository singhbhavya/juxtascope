"""
JuxtaScope detect() — a faithful port of the ImmGen colon 16a-16d doublet
pipeline, generalized + with juxtaposition-pair labeling.

Scoring (ported):
  * signatures: rank genes per group, keep logFC>threshold (16a: cross vs rest;
    16c: within-compartment vs siblings, stripping shared pan-lineage genes).
  * contamination = 2nd-highest foreign signature score (16a).
  * threshold = DATA-DRIVEN PERCENTILE of the contamination distribution
    (16b cross p95; 16d within p99) -- NOT a fixed z-score.
  * doublet = high contamination AND oversized (area OR counts, same percentile).
  * ambiguous_embedded = high contamination, normal size (kept separate: IELs).
  * misclassified = strong foreign identity, low contamination, normal size.
"""
import numpy as np
import pandas as pd

from .signatures import (derive_signatures, derive_sibling_signatures,
                         score_signatures)
from .morphology import morphology_scores, classify_morph


def _total_counts(adata, counts_key):
    if counts_key and counts_key in adata.obs.columns:
        return adata.obs[counts_key].values.astype(float)
    X = adata.X
    return np.asarray(X.sum(axis=1)).ravel() if hasattr(X, "sum") else X.sum(1)


def _second_highest(S):
    """per row: index & value of top and 2nd-highest score."""
    order = np.argsort(-S, axis=1)
    ti, si = order[:, 0], order[:, 1]
    top = S[np.arange(len(S)), ti]; sec = S[np.arange(len(S)), si]
    return ti, si, top, sec


def _oversized(adata, pctile, area_key, counts_vals):
    over = np.zeros(adata.n_obs, bool)
    if area_key and area_key in adata.obs.columns:
        a = adata.obs[area_key].values.astype(float)
        over |= a > np.nanquantile(a, pctile)
    if counts_vals is not None:
        over |= counts_vals > np.nanquantile(counts_vals, pctile)
    return over


def detect(adata, compartment_key, celltype_key=None,
           custom_broad=None, custom_fine=None,
           size_gate=True, embedded_as_doublet=False,
           cross_pctile=0.95, within_pctile=0.99, size_pctile=None,
           min_lfc_cross=1.0, min_lfc_within=0.5,
           sig_exclude=None, area_key="cell_area", counts_key="transcript_counts",
           n_genes=25, morphology=False, morph_extreme=3.5,
           inplace=True, save_to=None, verbose=True):
    # NOTE: results are added to the returned object; disk .h5ad is only written if save_to is set.
    """
    cross_pctile / within_pctile : percentile of the contamination distribution
        used as the cutoff (0.95 cross like 16b, 0.99 within like 16d). This is
        the KEY threshold -- data-driven, not a fixed z.
    size_pctile : percentile for 'oversized' (defaults to cross/within pctile).
    min_lfc_cross / min_lfc_within : logFC specificity filter for signatures
        (1.0 / 0.5, as 16a / 16c) -- keeps broad genes (Vim, Cd74) OUT.
    sig_exclude : {group: {genes}} to force-drop bleed genes (like 16c SIG_EXCLUDE).
    """
    ad = adata if inplace else adata.copy()
    counts_vals = _total_counts(ad, counts_key)
    uns = {}

    # ============ CROSS-COMPARTMENT (16a/16b) ============
    if verbose: print("[juxtascope] cross-compartment signatures (16a) ...")
    sig_broad = custom_broad or derive_signatures(
        ad, compartment_key, n_genes=n_genes, min_lfc=min_lfc_cross,
        exclude=sig_exclude, verbose=verbose)
    cols = score_signatures(ad, sig_broad, prefix="jsb_")
    groups_c = [c.replace("jsb_", "") for c in cols]
    Sc = ad.obs[cols].values
    ti, si, top_c, sec_c = _second_highest(Sc)
    own_c = ad.obs[compartment_key].astype(str).values
    topsig_c = np.array([groups_c[i] for i in ti])
    secsig_c = np.array([groups_c[i] for i in si])
    # contamination = 2nd-highest score; cutoff = its data-driven percentile
    contam_cut_c = np.nanquantile(sec_c, cross_pctile)
    hi_c = sec_c > contam_cut_c
    sp_c = size_pctile if size_pctile is not None else cross_pctile
    over_c = _oversized(ad, sp_c, area_key, counts_vals) if size_gate else np.ones(ad.n_obs, bool)
    # foreign identity for a cell = its top signature if that disagrees with its
    # own compartment, else the 2nd (the contaminating one)
    foreign_c = np.where(topsig_c != own_c, topsig_c, secsig_c)

    ad.obs["js_own"] = pd.Categorical(own_c)
    ad.obs["js_foreign"] = pd.Categorical(foreign_c)
    ad.obs["js_contamination"] = sec_c
    if verbose: print(f"    contamination p{int(cross_pctile*100)} cutoff = {contam_cut_c:.3f}")

    # ============ WITHIN-COMPARTMENT (16c/16d) ============
    hi_w = np.zeros(ad.n_obs, bool); over_w = np.zeros(ad.n_obs, bool)
    own_w = own_c.copy(); foreign_w = foreign_c.copy()
    sib_sigs = None
    if celltype_key and (custom_fine is not None or celltype_key in ad.obs):
        if verbose: print("[juxtascope] within-compartment sibling signatures (16c) ...")
        sib_sigs = custom_fine or derive_sibling_signatures(
            ad, compartment_key, celltype_key, n_genes=max(20, n_genes-5),
            min_lfc=min_lfc_within, exclude=sig_exclude, verbose=verbose)
        own_w = ad.obs[celltype_key].astype(str).values
        sp_w = size_pctile if size_pctile is not None else within_pctile
        over_w_all = _oversized(ad, sp_w, area_key, counts_vals) if size_gate else np.ones(ad.n_obs, bool)
        # score within each compartment against its siblings
        for comp, csigs in sib_sigs.items():
            m = (ad.obs[compartment_key].astype(str) == comp).values
            if m.sum() == 0 or len(csigs) < 2: continue
            sub = ad[m].copy()
            wcols = score_signatures(sub, csigs, prefix="jsw_")
            wgroups = [c.replace("jsw_", "") for c in wcols]
            Sw = sub.obs[wcols].values
            _, wsi, _, wsec = _second_highest(Sw)
            cut = np.nanquantile(wsec, within_pctile)  # per-compartment cutoff
            idx = np.where(m)[0]
            hi_w[idx] = wsec > cut
            over_w[idx] = over_w_all[idx]
            foreign_w_sub = np.array([wgroups[j] for j in wsi])
            foreign_w[idx] = foreign_w_sub
            if verbose: print(f"    {comp}: within-contam p{int(within_pctile*100)} = {cut:.3f}")

    # ============ COMBINE (priority, like 16d) ============
    status = np.array(["clean"]*ad.n_obs, dtype=object)
    level = np.array([""]*ad.n_obs, dtype=object)

    # misclassified: strong foreign single identity, low contam, normal size
    strong_top = top_c > np.nanquantile(top_c, 0.50)
    misclass = (topsig_c != own_c) & ~hi_c & ~over_c & strong_top
    status[misclass] = "misclassified"

    # embedded: high contam, normal size (real cell, neighbor bleed)
    embedded = (hi_c & ~over_c) | (hi_w & ~over_w)
    status[embedded] = "ambiguous_embedded"

    # doublet: high contam AND oversized
    dbl_w = hi_w & over_w
    dbl_c = hi_c & over_c
    status[dbl_w] = "doublet"; level[dbl_w] = "within"
    status[dbl_c] = "doublet"; level[dbl_c] = "cross"

    if embedded_as_doublet:
        status[embedded] = "doublet"

    # pair labels (own+foreign at whichever level fired)
    own_final = own_c.copy(); for_final = foreign_c.copy()
    wmask = (level == "within")
    own_final[wmask] = own_w[wmask]; for_final[wmask] = foreign_w[wmask]
    flagged = np.isin(status, ["doublet", "ambiguous_embedded"])
    pair = np.array([""]*ad.n_obs, dtype=object)
    for i in np.where(flagged)[0]:
        pair[i] = f"Juxtaposed {own_final[i]} - {for_final[i]}"

    ad.obs["js_status"] = pd.Categorical(status)
    ad.obs["js_level"] = pd.Categorical(level)
    ad.obs["js_pair"] = pd.Categorical(pair)

    # co-expressed markers per flagged cell (from foreign broad program)
    from .transcriptomic import coexpressed_markers
    coexpr = np.array([""]*ad.n_obs, dtype=object)
    for fg in np.unique(for_final[flagged]) if flagged.any() else []:
        mask = flagged & (for_final == fg)
        mk = coexpressed_markers(ad, mask, fg, sig_broad, top=5)
        for i in np.where(mask)[0]: coexpr[i] = ",".join(mk)
    ad.obs["js_coexpr_markers"] = pd.Categorical(coexpr)

    if morphology:
        mo = morphology_scores(ad, celltype_key or compartment_key, verbose=verbose)
        ad.obs["js_morph_flag"] = pd.Categorical(
            classify_morph(mo["morph_score"].values, morph_extreme, morph_extreme))

    uns["signatures_broad"] = {k: list(v) for k, v in sig_broad.items()}
    if sib_sigs:
        uns["signatures_within"] = {c: {g: list(gs) for g, gs in cs.items()}
                                    for c, cs in sib_sigs.items()}
    uns["params"] = dict(size_gate=size_gate, embedded_as_doublet=embedded_as_doublet,
                         cross_pctile=cross_pctile, within_pctile=within_pctile,
                         min_lfc_cross=min_lfc_cross, min_lfc_within=min_lfc_within)
    ad.uns["juxtascope"] = uns

    if verbose:
        vc = pd.Series(status).value_counts()
        print("[juxtascope] status:")
        for k, v in vc.items(): print(f"    {k}: {v:,} ({v/ad.n_obs:.1%})")

    if save_to:
        # scrub index/columns to plain object dtype so pyarrow-backed strings
        # don't break the .h5ad write
        for _df in (ad.obs, ad.var):
            _df.index = pd.Index(np.asarray(_df.index.tolist(), dtype=object),
                                 name=_df.index.name)
        ad.write_h5ad(save_to)
        if verbose: print(f"[juxtascope] wrote annotated object -> {save_to}")
    return ad
