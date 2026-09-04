"""
JuxtaScope — juxtaposition & doublet detection for segmentation-based spatial
transcriptomics.

Flags cells whose transcriptome mixes two identities (adjacent/juxtaposed cells
or segmentation merges), LABELS each by the pair ("Juxtaposed T cell -
Macrophage"), and reports which markers drove every call. Works off the user's
own broad + granular annotations OR reference/Pointillist predictions, with
optional user-supplied marker signatures.

    import scanpy as sc, juxtascope as js
    ad = sc.read_h5ad("tissue.h5ad")
    ad = js.detect(ad, compartment_key="Compartment", celltype_key="cell_type")
    js.report(ad, celltype_key="cell_type", compartment_key="Compartment",
              outdir="juxtascope_report")
"""
from .detect import detect
from .report import (report, pair_summary, markers_table, celltype_metrics)
from .signatures import derive_signatures, derive_sibling_signatures
from .morphology import morphology_scores, check_morphology_discriminates

__version__ = "0.2.0"
__all__ = ["detect","report","pair_summary","markers_table","celltype_metrics",
           "derive_signatures","morphology_scores","check_morphology_discriminates"]
