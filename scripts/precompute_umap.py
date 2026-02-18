#!/usr/bin/env python
"""
Precompute UMAP for all *_final.h5ad files that are missing X_umap,
and write X_umap back into each file in-place.

Run once from the repo root:
    python scripts/precompute_umap.py

Uses the same parameters as the old background computation:
    n_hvgs=1000, n_pcs=20, n_neighbors=15
"""

import glob
import os
import sys
import warnings

warnings.filterwarnings("ignore")

import anndata as ad
import scanpy as sc

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
N_HVGS = 1000
N_PCS = 20
N_NEIGHBORS = 15


def compute_umap(path: str) -> None:
    name = os.path.basename(path)
    print(f"\n{'='*60}")
    print(f"Processing: {name}")

    adata = ad.read_h5ad(path)

    if "X_umap" in adata.obsm:
        print(f"  → X_umap already present, skipping.")
        return

    if adata.X is None or adata.shape[1] < 2:
        print(f"  → No expression matrix, skipping.")
        return

    print(f"  shape: {adata.shape}")

    # Work on a copy so the original adata stays intact until we write
    _ad = adata.copy()

    sc.pp.normalize_total(_ad, target_sum=1e4)
    sc.pp.log1p(_ad)

    n_hvgs = min(N_HVGS, _ad.shape[1])
    print(f"  HVGs: {n_hvgs}")
    sc.pp.highly_variable_genes(_ad, n_top_genes=n_hvgs, flavor="seurat",
                                 subset=True, inplace=True)

    sc.pp.scale(_ad, max_value=10)

    n_pcs = min(N_PCS, _ad.shape[1] - 1)
    print(f"  PCA: {n_pcs} components")
    sc.tl.pca(_ad, n_comps=n_pcs)

    n_pcs_used = min(n_pcs, _ad.obsm["X_pca"].shape[1])
    sc.pp.neighbors(_ad, n_pcs=n_pcs_used, n_neighbors=N_NEIGHBORS)

    print(f"  Computing UMAP…")
    sc.tl.umap(_ad)

    # Copy result into original adata and write back to disk
    adata.obsm["X_umap"] = _ad.obsm["X_umap"]
    print(f"  Writing X_umap ({adata.obsm['X_umap'].shape}) back to {name}…")
    adata.write_h5ad(path)
    print(f"  ✓ Done.")


def main():
    files = sorted(glob.glob(os.path.join(DATA_DIR, "*_final.h5ad")))
    if not files:
        print(f"No *_final.h5ad files found in {DATA_DIR}")
        sys.exit(1)

    print(f"Found {len(files)} files in {os.path.abspath(DATA_DIR)}")
    for path in files:
        compute_umap(path)

    print(f"\n{'='*60}")
    print("All done.")


if __name__ == "__main__":
    main()
