import os
import threading
import pandas as pd
import numpy as np
import scanpy as sc
from metalotl._constants import DATA_DIR, GENE_ANNOTATION


# In-memory cache: filepath -> (adata, gene_choices, mtime)
_data_cache: dict = {}
_cache_lock = threading.Lock()

# Pre-built numpy array for vectorised rename (built once at import)
_amex_keys = np.array(list(GENE_ANNOTATION.keys()), dtype=object)
_amex_vals = np.array(list(GENE_ANNOTATION.values()), dtype=object)
_amex_map = dict(zip(_amex_keys, _amex_vals))


def _rename_amex_var_names(adata) -> None:
    """Rename raw AMEX gene IDs in adata.var_names to annotated names in-place.

    Uses a pre-built dict for O(1) per-gene lookup; only rewrites var_names
    when at least one name actually changed (avoids unnecessary copies).
    """
    current = adata.var_names.tolist()
    new_names = [_amex_map.get(g, g) for g in current]
    if new_names != current:
        adata.var_names = pd.Index(new_names)
        adata.var_names_make_unique()


def _build_gene_choices(adata) -> dict:
    """Build the selectize choices dict from adata.var_names — called once per file."""
    names = adata.var_names.tolist()
    choices = {'': ''}
    choices.update(zip(names, names))
    return choices


def _load_h5ad(filepath: str):
    """Load an h5ad fully into memory with mtime-based caching.

    Thread-safe: a per-file lock prevents two reactive effects from
    reading the same file simultaneously on first load.
    Returns (adata, gene_choices).
    """
    try:
        mtime = os.path.getmtime(filepath)
    except OSError:
        mtime = None

    with _cache_lock:
        entry = _data_cache.get(filepath)
        if entry is not None and (mtime is None or entry[2] == mtime):
            return entry[0], entry[1]

    # Read outside the lock so other files can be loaded concurrently
    adata = sc.read_h5ad(filepath)
    _rename_amex_var_names(adata)
    choices = _build_gene_choices(adata)

    with _cache_lock:
        _data_cache[filepath] = (adata, choices, mtime)
    return adata, choices


def _resolve_filepath(name: str):
    """Return the resolved Path for a dataset name, or None if not found."""
    stem = name.removesuffix('_final')
    for suffix in ('_final.h5ad', '.h5ad'):
        p = DATA_DIR / (stem + suffix)
        if p.exists():
            return p
    return None


def _prewarm_cache() -> None:
    """Load all *_final.h5ad files into the cache in a background thread."""
    for p in sorted(DATA_DIR.glob('*_final.h5ad')):
        try:
            _load_h5ad(str(p))
        except Exception as e:
            print(f"[metalotl] prewarm failed for {p.name}: {e}")


# Kick off background prewarming immediately at import time
threading.Thread(target=_prewarm_cache, daemon=True, name="metalotl-prewarm").start()


def clear_data_cache():
    """Clear the in-memory AnnData cache."""
    with _cache_lock:
        _data_cache.clear()


def get_data(input_or_name):
    """Return an AnnData for the given dataset name or Shiny inputs object."""
    try:
        name = input_or_name.select_dataset() if hasattr(input_or_name, 'select_dataset') else input_or_name
    except Exception:
        name = input_or_name

    if not name:
        return None

    p = _resolve_filepath(str(name))
    if p is None:
        print(f"[metalotl] dataset not found: {name}")
        return None

    try:
        adata, _ = _load_h5ad(str(p))
        return adata
    except Exception as e:
        print(f"[metalotl] load error for {name}: {e}")
        return None


def get_gene_choices(input_or_name) -> dict:
    """Return the pre-built gene selectize choices dict for the given dataset."""
    try:
        name = input_or_name.select_dataset() if hasattr(input_or_name, 'select_dataset') else input_or_name
    except Exception:
        name = input_or_name

    if not name:
        return {}

    p = _resolve_filepath(str(name))
    if p is None:
        return {}

    try:
        _, choices = _load_h5ad(str(p))
        return choices
    except Exception:
        return {}
