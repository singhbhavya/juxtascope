# JuxtaScope

> **Status: v0.2 beta.** 

**Juxtaposition & doublet detection** for segmentation-based spatial transcriptomics (Xenium, MERFISH, CosMx). This is meant to be run it *after* annotation to find cells whose transcriptome mixes two identities.The pipeline **labels each by the mixing/juxtaposed pair** (e.g. *Juxtaposed T cell - Macrophage*).

This will work with your **own broad + granular annotations** (as this is a POST-ANNOTATION pipeline), and accepts **custom marker signatures**. Marker genes are identified for broad cellular identities (epithelial, vascular, immune, etc.) as well as within-compartment cellular identities (B cells, T cells, etc) to identify cross-compartment doublets (e.g., Immmune - Vascular) and within-compartment doublets (e.g., Dendritic cell - T cell). 

## Methodology:

A cell is only called a **confident doublet** when it is **both**:

1. **high in a foreign identity's markers** (co-expresses another cell type's program), **AND**
2. **oversized** — larger cell area or more transcripts than the pXX (default: p95 for cross-compartment and p99 for within-compartment) cutoff.

This AND-gate is deliberate, and it separates two physically different things that co-expression *alone* cannot tell apart:

- **`doublet`** — high co-expression **and oversized**. Two cells merged into one segment: a segmentation error producing one chimeric "cell." In most cases, these should be excluded, unless biologically interesting (.
- **`ambiguous_embedded`** — high co-expression but **normal size**. One real, correctly-sized cell whose transcriptome picks up a *physically adjacent neighbor's* signal because it sits embedded in another tissue — the classic  case being an **intraepithelial lymphocyte (IEL)**: a genuine single T cell living inside the epithelial layer. **This is not a segmentation error — it is often real biology you want to KEEP.** If you labeled it a doublet and filtered doublets out, you would be deleting real IELs.
- **`misclassified`** — a strong *foreign* single identity, low co-expression, normal size: one cell that simply carries the wrong label.

Because embedded cells are frequently real, JuxtaScope keeps them **separate from doublets by default**. Set `embedded_as_doublet=True` if you want the stricter filter that treats them as doublets too. Preserve the information by default; collapse only when you choose to.

Size gating can be turned off entirely (`size_gate=False`) — then high co-expression alone flags a cell as `ambiguous_embedded`.

## Two levels & cutoffs (mirroring the original pipeline)

- **Cross-compartment** (broad key) — default contamination+size cutoff **p95**.
- **Within-compartment** (granular key) — default **p99**.

Both, plus the co-expression z thresholds, are user-settable:
`cross_pctile`, `within_pctile`, `cross_z`, `within_z`, `size_pctile`.

## Quick start

```python
import scanpy as sc, juxtascope as js
ad = sc.read_h5ad("tissue.h5ad")

ad = js.detect(ad,
    compartment_key="Compartment",   # broad identity  -> cross-compartment (p95)
    celltype_key="cell_type",        # granular identity -> within (p99)
    size_gate=True,                  # require oversized for a doublet (default)
    embedded_as_doublet=False,       # keep IEL-type embedded cells separate (default)
    # custom_broad={...}, custom_fine={...},   # optional marker overrides
)

js.report(ad, celltype_key="cell_type", compartment_key="Compartment",
          outdir="juxtascope_report")
```

## Outputs

**obs columns** (`js_`): `js_status` (clean / doublet / ambiguous_embedded /
misclassified), `js_pair` (the labeled mixing pair), `js_level` (cross/within),
`js_own`, `js_foreign`, `js_cross_score`, `js_coexpr_markers`.

**report/ folder:**
- `umap_broad`, `umap_granular`, `umap_juxtaposition` (png+pdf)
- `juxtaposition_pairs.csv` — which type-pairs mix most, counts, co-expressed markers
- `chosen_markers.csv` — markers used to delineate each broad & granular type
- `celltype_metrics.csv` — per cell type AND per pair: mean cell area, nucleus count, transcripts
- `metric_*.png` — bar charts of those metrics

## Custom markers

```python
ad = js.detect(ad, compartment_key="Compartment", celltype_key="cell_type",
    custom_broad={"Immune":["Ptprc","Cd3e"], "Epithelial":["Epcam","Krt8"], ...},
    custom_fine={"CD8 T":["Cd8a","Gzmb"], "Macrophage":["Cd68","Adgre1"], ...})
```

## Morphology (off by default)

Area/nucleus *anomaly* scoring is separate from the size gate and is OFF by
default — in Xenium colon it did not add doublet signal beyond co-expression +
the size gate. Check before enabling:
```python
js.check_morphology_discriminates(ad, celltype_key="cell_type",
                                  reference_flag_key="doublet_status")
```

## Install
```bash
pip install -e .        # from this directory (contains setup.py)
```
