"""
Morphology-based doublet evidence, judged PER CELL TYPE.

Every feature (cell area, nucleus count, nuclear:cytoplasmic ratio) is converted
to a within-cell-type robust z-score (median / MAD). A cell is morphologically
anomalous only relative to OTHER CELLS OF ITS OWN PREDICTED TYPE -- so
multinucleate types (myofibers, osteoclasts) or large types (macrophages) set
their own baseline and are never flagged merely for being large or multinucleate.

Directionality: for doublet detection we care about cells that are anomalously
LARGE / EXTRA-nucleated / mis-proportioned for their type (segmentation merges),
not cells that are smaller than typical. Scores are therefore one-sided (high).
"""
import numpy as np
import pandas as pd


# candidate morphology columns (auto-detected; all optional)
_AREA_KEYS    = ["cell_area","area","cell_size","Area","segmentation_area"]
_NUC_AREA_KEYS= ["nucleus_area","nuclear_area","nuc_area"]
_NUC_CNT_KEYS = ["nucleus_count","n_nuclei","nuclei_count","num_nuclei"]


def _first_present(obs, keys):
    for k in keys:
        if k in obs.columns:
            return k
    return None


def _robust_z(values, labels, one_sided=True):
    """within-group robust z: (x - group_median)/group_MAD, per label group.
    one_sided keeps only positive deviations (anomalously HIGH), clips neg to 0."""
    v = np.asarray(values, float)
    z = np.full(len(v), np.nan)
    s = pd.Series(v)
    for g, idx in pd.Series(range(len(v))).groupby(np.asarray(labels)).groups.items():
        gi = np.asarray(idx)
        vg = v[gi]
        ok = np.isfinite(vg)
        if ok.sum() < 5:
            continue
        med = np.median(vg[ok])
        mad = np.median(np.abs(vg[ok] - med)) * 1.4826  # ~ std for normal
        if mad == 0:
            mad = np.nanstd(vg[ok]) or 1.0
        zz = (vg - med) / mad
        if one_sided:
            zz = np.clip(zz, 0, None)
        z[gi] = zz
    return z


def morphology_scores(adata, celltype_key, area_key=None, nuc_area_key=None,
                      nuc_count_key=None, verbose=True):
    """Return a DataFrame of per-cell within-type morphology anomaly z-scores
    plus a combined `morph_score`. Missing features are skipped gracefully.

    Columns (any present): area_z, ncratio_z, nuccount_z, morph_score, morph_flag
    morph_flag: 'na' | 'clean' | 'anomalous' | 'extreme'
    """
    obs = adata.obs
    labels = obs[celltype_key].astype(str).values
    area_key = area_key or _first_present(obs, _AREA_KEYS)
    nuc_area_key = nuc_area_key or _first_present(obs, _NUC_AREA_KEYS)
    nuc_count_key = nuc_count_key or _first_present(obs, _NUC_CNT_KEYS)

    out = pd.DataFrame(index=obs.index)
    feats = []

    if area_key:
        out["area_z"] = _robust_z(obs[area_key].values, labels, one_sided=True)
        feats.append("area_z")
        if verbose: print(f"  [morph] area from '{area_key}'")

    # nuclear:cytoplasmic ratio (needs both nucleus_area and cell area)
    if nuc_area_key and area_key:
        with np.errstate(divide="ignore", invalid="ignore"):
            ncr = obs[nuc_area_key].values.astype(float) / obs[area_key].values.astype(float)
        # anomaly in EITHER direction matters for N:C (wrong cell-type proportion),
        # so use two-sided magnitude
        z = _robust_z(ncr, labels, one_sided=False)
        out["ncratio_z"] = np.abs(z)
        feats.append("ncratio_z")
        if verbose: print(f"  [morph] N:C ratio from '{nuc_area_key}'/'{area_key}'")

    if nuc_count_key:
        # judged within type: multinucleate types (myofibers) have high baseline,
        # only WITHIN-type excess counts. one-sided (extra nuclei = merge).
        out["nuccount_z"] = _robust_z(obs[nuc_count_key].values, labels, one_sided=True)
        feats.append("nuccount_z")
        if verbose: print(f"  [morph] nucleus count from '{nuc_count_key}'")

    if not feats:
        if verbose: print("  [morph] no morphology columns found; skipping")
        out["morph_score"] = np.nan
        out["morph_flag"] = "na"
        return out

    # combined score = max anomaly across available features (worst offender),
    # which is robust to having 1, 2 or 3 features present.
    out["morph_score"] = out[feats].max(axis=1, skipna=True)
    return out


def classify_morph(morph_score, anomalous=2.0, extreme=3.5):
    """flag thresholds on the within-type z (not absolute units)."""
    f = np.full(len(morph_score), "clean", dtype=object)
    s = np.asarray(morph_score, float)
    f[~np.isfinite(s)] = "na"
    f[s >= anomalous] = "anomalous"
    f[s >= extreme] = "extreme"
    return f


def check_morphology_discriminates(adata, celltype_key, reference_flag_key,
                                   suspect_values=None, verbose=True):
    """Sanity check BEFORE trusting morphology: do cells already flagged as
    suspect (by transcriptomics or a prior pipeline) actually have higher
    within-type morphology anomaly than clean cells?

    If the two distributions overlap (AUC ~0.5), morphology adds no signal in
    this dataset and should stay disabled -- this is what we observed in Xenium
    colon. Returns a dict of AUC per feature; AUC >~0.6 means some signal.

    reference_flag_key : an obs column marking suspect cells (e.g. a prior
                         doublet_status). suspect_values : which of its values
                         count as suspect (default: anything != 'clean').
    """
    import numpy as np
    obs = adata.obs
    m = morphology_scores(adata, celltype_key, verbose=False)
    ref = obs[reference_flag_key].astype(str)
    if suspect_values is None:
        suspect = (ref != "clean").values
    else:
        suspect = ref.isin(suspect_values).values

    def auc(score, pos):
        s = np.asarray(score, float); ok = np.isfinite(s)
        s, pos = s[ok], pos[ok]
        if pos.sum()==0 or (~pos).sum()==0: return np.nan
        # rank-based AUC
        order = np.argsort(s)
        ranks = np.empty(len(s)); ranks[order] = np.arange(1, len(s)+1)
        n1 = pos.sum(); n0 = (~pos).sum()
        return (ranks[pos].sum() - n1*(n1+1)/2) / (n1*n0)

    results={}
    for col in [c for c in ["area_z","ncratio_z","nuccount_z","morph_score"] if c in m.columns]:
        a = auc(m[col].values, suspect)
        results[col]=a
        if verbose:
            verdict = ("no signal" if not np.isfinite(a) or a<0.55
                       else "weak" if a<0.62 else "some signal")
            print(f"  {col:12s} AUC={a:.3f}  -> {verdict}")
    if verbose:
        best=max([v for v in results.values() if np.isfinite(v)] or [0])
        print(f"\n  {'ENABLE morphology' if best>=0.62 else 'KEEP morphology OFF'} "
              f"(best AUC {best:.3f}); AUC~0.5 means area/nuclei do not separate "
              f"doublets in this data.")
    return results
