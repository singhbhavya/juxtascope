"""
JuxtaScope reporting: UMAPs, juxtaposition-pair summary, chosen markers, and
per-cell-type metric tables + bar charts. Call report() after detect().
"""
import os
import numpy as np
import pandas as pd


# OPTIONAL preferred hues for commonly-seen compartment names, so figures stay
# consistent when these names appear. Any compartment NOT in here (i.e. another
# tissue's naming) is assigned a distinct hue automatically -- nothing is
# hardcoded as required.
PREFERRED_HUE = {"Epithelial":"#1f77b4","Immune":"#d62728","Stroma":"#ff7f0e",
                 "Vascular":"#2ca02c","Neural":"#9467bd"}


def _assign_compartment_hues(compartments):
    """Build {compartment: hex hue} for WHATEVER compartments the data has.
    Known names get their preferred hue; unknown names get evenly-spaced hues
    from HSV space that avoid collisions with the preferred ones."""
    import colorsys, matplotlib.colors as mcolors
    comps = list(compartments)
    used_hues = set()
    hue_map = {}
    # 1) assign preferred hues to any recognized names
    for c in comps:
        if c in PREFERRED_HUE:
            hue_map[c] = PREFERRED_HUE[c]
            h,_,_ = colorsys.rgb_to_hsv(*mcolors.to_rgb(PREFERRED_HUE[c]))
            used_hues.add(round(h,2))
    # 2) auto-assign the rest to evenly spaced hues not near the used ones
    unknown = [c for c in comps if c not in hue_map]
    if unknown:
        # candidate hues spread around the wheel
        cands = [i/ max(len(comps),1) for i in range(len(comps)*2)]
        cands = [h for h in cands if all(abs(h-u) > 0.06 for u in used_hues)]
        k = 0
        for c in unknown:
            h = cands[k % len(cands)] if cands else (k/ max(len(unknown),1))
            k += 1
            hue_map[c] = mcolors.to_hex(colorsys.hsv_to_rgb(h, 0.65, 0.85))
    return hue_map


# module-level cache filled per-call by the palette functions below
COMP_HUE = dict(PREFERRED_HUE)   # back-compat default; overwritten at runtime


def _shades(hue, n):
    """n perceptually-spaced shades of a base hue, dark -> light."""
    import matplotlib.colors as mcolors
    import numpy as _np
    if n <= 0: return []
    base = _np.array(mcolors.to_rgb(hue))
    if n == 1: return [mcolors.to_hex(base)]
    out = []
    for f in _np.linspace(-0.4, 0.55, n):
        c = _np.clip(base*(1+f),0,1) if f<=0 else _np.clip(base+(1-base)*f,0,1)
        out.append(mcolors.to_hex(c))
    return out


def _compartment_palette(obs, comp_key):
    cats = [c for c in pd.Categorical(obs[comp_key].astype(str)).categories]
    return _assign_compartment_hues(cats)


def _hier_palette(obs, comp_key, child_key):
    """each child (Middle/MiddleSub) -> a shade of its Compartment's hue,
    with a per-group hue nudge so siblings stay distinguishable. Compartment
    hues are assigned from WHATEVER compartments the data has (tissue-agnostic)."""
    import colorsys, matplotlib.colors as mcolors
    palette = {}
    if comp_key not in obs.columns or child_key not in obs.columns:
        return _compartment_palette(obs, child_key)
    comps = [c for c in pd.Categorical(obs[comp_key].astype(str)).categories]
    hue_map = _assign_compartment_hues(comps)
    comp_of = (obs.groupby(child_key, observed=True)[comp_key]
               .agg(lambda s: s.astype(str).mode().iloc[0]).to_dict())
    by_comp = {}
    for ch in pd.Categorical(obs[child_key].astype(str)).categories:
        by_comp.setdefault(comp_of.get(ch,"other"), []).append(ch)
    for comp, kids in by_comp.items():
        base = hue_map.get(comp, "#888888")
        r,g,b = mcolors.to_rgb(base); h,s,v = colorsys.rgb_to_hsv(r,g,b)
        for j, ch in enumerate(sorted(kids)):
            hj = (h + (j-len(kids)/2)*0.045) % 1.0
            palette[ch] = mcolors.to_hex(colorsys.hsv_to_rgb(hj, s, v))
    return palette


def _counts_per_cell(adata):
    X = adata.X
    if hasattr(X, "sum"):
        tc = np.asarray(X.sum(axis=1)).ravel()
    else:
        tc = X.sum(1)
    return tc


def pair_summary(adata, signatures_key="signatures_broad", top_markers=5):
    """(2) which cell-type PAIRS are most often juxtaposed, with counts and the
    co-expressed markers driving each. Returns a DataFrame."""
    obs = adata.obs
    fl = obs["js_status"].astype(str).isin(["doublet","ambiguous_embedded"])
    if fl.sum() == 0:
        return pd.DataFrame(columns=["pair","n","own","foreign","coexpr_markers"])
    sub = obs[fl].copy()
    # canonical unordered pair so A-B and B-A aggregate
    a = sub["js_own"].astype(str); b = sub["js_foreign"].astype(str)
    canon = [" | ".join(sorted([x,y])) for x,y in zip(a,b)]
    sub["_canon"] = canon
    rows = []
    for cp, g in sub.groupby("_canon", observed=True):
        # most common co-expressed markers across these cells
        mk = (g["js_coexpr_markers"].astype(str).str.split(",").explode())
        mk = mk[mk.str.len()>0].value_counts().head(top_markers).index.tolist()
        own_ex = g["js_own"].astype(str).mode().iloc[0]
        for_ex = g["js_foreign"].astype(str).mode().iloc[0]
        rows.append({"pair":cp, "n":len(g), "example_own":own_ex,
                     "example_foreign":for_ex, "coexpr_markers":",".join(mk)})
    return pd.DataFrame(rows).sort_values("n", ascending=False)


def markers_table(adata):
    """(3) which markers were chosen per broad (and fine) type."""
    uns = adata.uns.get("juxtascope", {})
    rows=[]
    for grp, genes in uns.get("signatures_broad", {}).items():
        rows.append({"level":"broad","group":grp,"markers":",".join(genes)})
    for comp, csigs in uns.get("signatures_within", {}).items():
        for grp, genes in csigs.items():
            rows.append({"level":f"within:{comp}","group":grp,
                         "markers":",".join(genes)})
    return pd.DataFrame(rows)


def celltype_metrics(adata, celltype_key, extra_group="js_pair"):
    """(5) per-cell-type metrics: avg cell area, nucleus count, transcripts, etc.
    Also computed per juxtaposition-pair label."""
    obs = adata.obs.copy()
    obs["_n_counts"] = _counts_per_cell(adata)
    feats = {}
    for name, keys in {"cell_area":["cell_area","area"],
                       "nucleus_area":["nucleus_area"],
                       "nucleus_count":["nucleus_count","n_nuclei"]}.items():
        for k in keys:
            if k in obs.columns: feats[name]=k; break
    def agg(gkey):
        g = obs.groupby(gkey, observed=True)
        d = pd.DataFrame({"n_cells": g.size()})
        d["mean_transcripts"] = g["_n_counts"].mean().round(1)
        for name,col in feats.items():
            d[f"mean_{name}"] = g[col].mean().round(1)
        return d.reset_index().rename(columns={gkey:"group"})
    per_type = agg(celltype_key); per_type["kind"]="cell_type"
    pairs = obs[obs["js_status"].astype(str).isin(["doublet","ambiguous_embedded"])]
    per_pair = pd.DataFrame()
    if len(pairs):
        gp = pairs.groupby("js_pair", observed=True)
        per_pair = pd.DataFrame({"group":gp.size().index,"n_cells":gp.size().values})
        per_pair["mean_transcripts"]=gp["_n_counts"].mean().round(1).values
        for name,col in feats.items():
            per_pair[f"mean_{name}"]=gp[col].mean().round(1).values
        per_pair["kind"]="juxtaposition"
    return pd.concat([per_type, per_pair], ignore_index=True)


def report(adata, celltype_key, compartment_key, outdir="juxtascope_report",
           make_umaps=True, verbose=True):
    """Write the full report: CSVs + UMAPs + metric bar charts."""
    os.makedirs(outdir, exist_ok=True)
    import matplotlib as mpl; mpl.use("Agg"); mpl.rcParams["pdf.fonttype"]=42
    import matplotlib.pyplot as plt

    # --- CSVs ---
    ps = pair_summary(adata); ps.to_csv(f"{outdir}/juxtaposition_pairs.csv", index=False)
    mt = markers_table(adata); mt.to_csv(f"{outdir}/chosen_markers.csv", index=False)
    cm = celltype_metrics(adata, celltype_key)
    cm.to_csv(f"{outdir}/celltype_metrics.csv", index=False)
    if verbose:
        print(f"[report] wrote juxtaposition_pairs.csv ({len(ps)} pairs), "
              f"chosen_markers.csv ({len(mt)} groups), celltype_metrics.csv")
        if len(ps): print("[report] top juxtaposition pairs:\n",
                          ps.head(8).to_string(index=False))

    # --- (1) UMAPs ---
    if make_umaps and "X_umap" in adata.obsm:
        xy = adata.obsm["X_umap"]
        obs = adata.obs

        # ----- hierarchical color scheme (compartment hue -> shades) -----
        STATUS_COLORS = {"clean":"#dcdcdc", "ambiguous_embedded":"#ff7f0e",
                         "misclassified":"#9467bd", "doublet":"#d62728"}
        # draw order for status: clean underneath ... doublets on top  (#2)
        STATUS_ORDER = ["clean","ambiguous_embedded","misclassified","doublet"]
        FLAG_VALUES = {"doublet","ambiguous_embedded","misclassified"}

        def _plain_umap(color_key, title, fname, palette):
            labels = obs[color_key].astype(str).values
            colors = np.array([palette.get(l,"#cccccc") for l in labels],dtype=object)
            fig,ax=plt.subplots(figsize=(10,9))
            ax.scatter(xy[:,0],xy[:,1],s=2,c=list(colors),linewidths=0,rasterized=True)
            ax.set_title(title,fontsize=14); ax.axis("off"); ax.set_aspect("equal")
            uniq=[l for l in dict.fromkeys(labels) if l in palette]
            if len(uniq)<=30:
                h=[plt.Line2D([0],[0],marker='o',color='w',markerfacecolor=palette[l],
                   markersize=7,label=l) for l in sorted(uniq)]
                ax.legend(handles=h,loc="center left",bbox_to_anchor=(1,0.5),
                          fontsize=7,frameon=False,ncol=1 if len(uniq)<=25 else 2)
            plt.tight_layout()
            fig.savefig(f"{outdir}/{fname}.png",dpi=300,bbox_inches="tight")
            fig.savefig(f"{outdir}/{fname}.pdf",bbox_inches="tight"); plt.close()
            if verbose: print(f"[report] wrote {fname}")

        # broad = compartment base hues; granular = shades within compartment (#1)
        pal_broad = _compartment_palette(obs, compartment_key)
        _plain_umap(compartment_key, "Broad cell types", "umap_broad", pal_broad)
        pal_gran = _hier_palette(obs, compartment_key, celltype_key)
        _plain_umap(celltype_key, "Granular subtypes", "umap_granular", pal_gran)

        # (#2) status UMAP with explicit layering: clean, ambiguous, misclassified, doublet on top
        if "js_status" in obs:
            labels = obs["js_status"].astype(str).values
            fig,ax=plt.subplots(figsize=(10,9))
            for lvl in STATUS_ORDER:
                m = labels==lvl
                if not m.any(): continue
                ax.scatter(xy[m,0],xy[m,1],s=2,c=STATUS_COLORS.get(lvl,"#999999"),
                           linewidths=0,rasterized=True)
            ax.set_title("Juxtaposition / doublet detection",fontsize=14)
            ax.axis("off"); ax.set_aspect("equal")
            seen=[l for l in STATUS_ORDER if (labels==l).any()]
            h=[plt.Line2D([0],[0],marker='o',color='w',markerfacecolor=STATUS_COLORS[l],
               markersize=8,label=l) for l in seen]
            ax.legend(handles=h,loc="center left",bbox_to_anchor=(1,0.5),
                      fontsize=8,frameon=False)
            plt.tight_layout()
            fig.savefig(f"{outdir}/umap_juxtaposition.png",dpi=300,bbox_inches="tight")
            fig.savefig(f"{outdir}/umap_juxtaposition.pdf",bbox_inches="tight"); plt.close()
            if verbose: print("[report] wrote umap_juxtaposition")

            # (#3) compartment + granular colored, with ALL flagged cells greyed on top
            for tkey, pal, nm in [(compartment_key,pal_broad,"umap_broad_flags_grey"),
                                  (celltype_key,pal_gran,"umap_granular_flags_grey")]:
                lab = obs[tkey].astype(str).values
                flagged = obs["js_status"].astype(str).isin(FLAG_VALUES).values
                colors = np.array([pal.get(l,"#cccccc") for l in lab],dtype=object)
                colors[flagged] = "#4d4d4d"
                fig,ax=plt.subplots(figsize=(10,9))
                for m in (~flagged, flagged):        # clean first, grey flags on top
                    ax.scatter(xy[m,0],xy[m,1],s=2,c=list(colors[m]),
                               linewidths=0,rasterized=True)
                ax.set_title(f"{tkey} (flagged cells greyed)",fontsize=14)
                ax.axis("off"); ax.set_aspect("equal")
                plt.tight_layout()
                fig.savefig(f"{outdir}/{nm}.png",dpi=300,bbox_inches="tight")
                fig.savefig(f"{outdir}/{nm}.pdf",bbox_inches="tight"); plt.close()
                if verbose: print(f"[report] wrote {nm}")

    # --- (5) metric bar charts ---
    numcols=[c for c in cm.columns if c.startswith("mean_")]
    if numcols:
        ct = cm[cm["kind"]=="cell_type"].sort_values("group")
        for col in numcols:
            fig,ax=plt.subplots(figsize=(max(10,len(ct)*0.22),5))
            ax.bar(range(len(ct)),ct[col].values,color="#4a7ba6")
            ax.set_xticks(range(len(ct)));ax.set_xticklabels(ct["group"],rotation=90,fontsize=6)
            ax.set_ylabel(col); ax.set_title(f"{col} per cell type")
            plt.tight_layout()
            fig.savefig(f"{outdir}/metric_{col}.png",dpi=200,bbox_inches="tight");plt.close()
        if verbose: print(f"[report] wrote {len(numcols)} metric bar charts")
    print(f"[report] done -> {outdir}/")
