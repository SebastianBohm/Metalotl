"""Benchmark individual steps in the load pipeline."""
import time

import anndata as ad
import pandas as pd
from metalotl._constants import DATA_DIR, GENE_ANNOTATION

files = sorted(DATA_DIR.glob('Meta_*_final.h5ad'))

for p in files[:3]:
    print(f"\n{p.name}")

    # 1. Raw h5ad read
    t0 = time.perf_counter()
    adata = ad.read_h5ad(str(p))
    t1 = time.perf_counter()
    print(f"  read_h5ad:          {t1-t0:.3f}s  shape={adata.shape}")

    # 2. _rename_amex_var_names (list comprehension over GENE_ANNOTATION dict)
    t0 = time.perf_counter()
    new_names = [GENE_ANNOTATION.get(str(g), str(g)) for g in adata.var_names]
    idx = pd.Index(new_names)
    t1 = time.perf_counter()
    print(f"  rename list+Index:  {t1-t0:.3f}s  n_genes={len(new_names)}")

    # 3. var_names_make_unique
    t0 = time.perf_counter()
    adata.var_names = idx
    adata.var_names_make_unique()
    t1 = time.perf_counter()
    print(f"  make_unique:        {t1-t0:.3f}s")

    # 4. _build_gene_choices
    t0 = time.perf_counter()
    choices = {'': ''}
    for name in map(str, adata.var_names):
        choices[name] = name
    t1 = time.perf_counter()
    print(f"  build_choices:      {t1-t0:.3f}s  n={len(choices)}")

    # 5. Cache hit simulation (just dict lookup)
    t0 = time.perf_counter()
    _ = adata.var_names
    t1 = time.perf_counter()
    print(f"  cache hit (lookup): {t1-t0:.6f}s")

print("\nDone.")
