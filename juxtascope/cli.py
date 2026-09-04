"""Command-line entry: run juxtascope on an .h5ad and write results."""
import argparse, sys


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="juxtascope",
        description="Post-annotation doublet / mixed-population detection for "
                    "segmentation-based spatial transcriptomics.")
    p.add_argument("h5ad", help="input .h5ad (annotated)")
    p.add_argument("-o","--out", required=True, help="output .h5ad")
    p.add_argument("--compartment-key", required=True,
                   help="obs column with broad compartment labels")
    p.add_argument("--celltype-key", default=None,
                   help="obs column with fine cell-type labels "
                        "(defaults to compartment-key)")
    p.add_argument("--no-morphology", action="store_true",
                   help="disable morphology evidence")
    p.add_argument("--cross-z", type=float, default=2.0)
    p.add_argument("--morph-anomalous", type=float, default=2.0)
    p.add_argument("--morph-extreme", type=float, default=3.5)
    p.add_argument("--report-dir", default=None,
                   help="write full report (UMAPs, pairs, metrics) to this dir")
    p.add_argument("--report", default=None,
                   help="optional CSV summary of flags per cell type")
    a = p.parse_args(argv)

    import scanpy as sc
    from .detect import detect
    from .report import report as _report

    ad = sc.read_h5ad(a.h5ad)
    ad = detect(ad, compartment_key=a.compartment_key,
                celltype_key=a.celltype_key,
                morphology=not a.no_morphology,
                cross_z=a.cross_z, morph_anomalous=a.morph_anomalous,
                morph_extreme=a.morph_extreme)
    ad.write_h5ad(a.out)
    if a.report_dir:
        _report(ad, celltype_key=a.celltype_key or a.compartment_key,
                compartment_key=a.compartment_key, outdir=a.report_dir)
    print(f"[juxtascope] wrote {a.out}")

    if a.report:
        import pandas as pd
        ctk = a.celltype_key or a.compartment_key
        rep = (ad.obs.groupby(ctk, observed=True)["ds_status"]
               .value_counts(normalize=True).unstack(fill_value=0).round(3))
        rep.to_csv(a.report)
        print(f"[juxtascope] wrote report {a.report}")


if __name__ == "__main__":
    main()
