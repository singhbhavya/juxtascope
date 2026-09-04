"""
Transcriptomic juxtaposition/doublet evidence: does a cell co-express a FOREIGN
compartment (or cell-type) program in addition to its own? Records the top
foreign identity and the specific co-expressed markers, so downstream we can
label each flagged cell as "Juxtaposed <own> - <foreign>".
"""
import numpy as np
import pandas as pd
from .signatures import derive_signatures, score_signatures


def _z(obs, cols):
    return obs[cols].apply(lambda c: (c - c.mean())/(c.std() or 1))


def transcriptomic_scores(adata, group_key, signatures=None, n_genes=25,
                          prefix="js_", verbose=True):
    """Score each cell against ALL group signatures; return its own-group z, its
    top FOREIGN group + that group's z, and the markers of the foreign program
    the cell actually expresses (for labeling & reporting)."""
    obs = adata.obs
    if signatures is None:
        signatures = derive_signatures(adata, group_key, n_genes=n_genes, verbose=verbose)
    cols = score_signatures(adata, signatures, prefix=f"{prefix}sig_", verbose=verbose)
    grp_of = {f"{prefix}sig_{g}": g for g in signatures if f"{prefix}sig_{g}" in cols}
    Z = _z(obs, cols); zvals = Z.values
    zgroups = [grp_of[c] for c in cols]

    own = obs[group_key].astype(str).values
    own_z = np.full(len(obs), np.nan)
    ftop = np.array(["none"]*len(obs), dtype=object)
    fz = np.full(len(obs), np.nan)
    for i in range(len(obs)):
        row = zvals[i]
        try: oj = zgroups.index(own[i])
        except ValueError: oj = -1
        if oj>=0: own_z[i]=row[oj]
        bj,bv=-1,-np.inf
        for j,v in enumerate(row):
            if j==oj: continue
            if v>bv: bv,bj=v,j
        if bj>=0: ftop[i]=zgroups[bj]; fz[i]=bv

    out = pd.DataFrame(index=obs.index)
    out["own_group"]=own
    out["own_z"]=own_z
    out["foreign_group"]=ftop
    out["foreign_z"]=fz
    out["cross_score"]=np.clip(fz,0,None)
    return out, signatures


def coexpressed_markers(adata, cell_mask, foreign_group, signatures, top=5):
    """For a set of cells assigned to a foreign_group, return which of that
    group's signature genes they actually express most (the co-expression
    evidence). Used to annotate each juxtaposition pair."""
    genes = [g for g in signatures.get(foreign_group,[]) if g in adata.var_names]
    if not genes or cell_mask.sum()==0:
        return []
    sub = adata[cell_mask, genes]
    X = sub.X
    mean = np.asarray(X.mean(axis=0)).ravel() if hasattr(X,"mean") else X.mean(0)
    order = np.argsort(mean)[::-1][:top]
    return [genes[i] for i in order]
