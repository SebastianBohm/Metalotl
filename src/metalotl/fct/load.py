import os
import pandas as pd
import scanpy as sc
from metalotl._constants import DATA_DIR, GENE_ANNOTATION


# In-memory cache: filepath -> (adata, gene_choices, mtime)
_data_cache: dict = {}


def _rename_amex_var_names(adata) -> None:
    """Rename any raw AMEX gene IDs in adata.var_names to annotated names in-place.

    The _final.h5ad files were saved before the notebook's gene-renaming step ran,
    so var_names may contain a mix of AMEX IDs and already-annotated names.
    We apply GENE_ANNOTATION to every var_name, falling back to the original name
    when no annotation exists (preserving already-annotated names unchanged).
    """
    new_names = [GENE_ANNOTATION.get(str(g), str(g)) for g in adata.var_names]
    adata.var_names = pd.Index(new_names)
    adata.var_names_make_unique()


def _build_gene_choices(adata) -> dict:
    """Build the selectize choices dict from adata.var_names — called once per file.

    At this point var_names are already annotated display names, so the label
    is just the name itself (no AMEX ID suffix needed).
    """
    choices = {'': ''}
    for name in map(str, adata.var_names):
        choices[name] = name
    return choices


def _load_h5ad(filepath):
    """Load an h5ad fully into memory with mtime-based caching.

    Full in-memory load is used (not backed='r') because gene-expression slicing
    on a backed dataset is ~7x slower, and these files load in <0.1s regardless.
    Returns (adata, gene_choices).
    """
    try:
        mtime = os.path.getmtime(filepath)
    except Exception:
        adata = sc.read_h5ad(filepath)
        _rename_amex_var_names(adata)
        return adata, _build_gene_choices(adata)

    entry = _data_cache.get(filepath)
    if entry is not None:
        cached_adata, cached_choices, cached_mtime = entry
        if cached_mtime == mtime:
            return cached_adata, cached_choices

    adata = sc.read_h5ad(filepath)
    _rename_amex_var_names(adata)
    choices = _build_gene_choices(adata)
    _data_cache[filepath] = (adata, choices, mtime)
    return adata, choices


def clear_data_cache():
    """Clear the in-memory AnnData cache."""
    _data_cache.clear()


def get_data(input_or_name):
    """Return an AnnData for the given dataset name or Shiny inputs object."""
    name = None
    try:
        if hasattr(input_or_name, 'select_dataset') and callable(input_or_name.select_dataset):
            name = input_or_name.select_dataset()
        else:
            name = input_or_name
    except Exception:
        name = input_or_name

    if not name:
        return None

    try:
        stem = name[:-len('_final')] if name.endswith('_final') else name
        final_path = DATA_DIR / (stem + '_final.h5ad')
        fallback_path = DATA_DIR / (stem + '.h5ad')
        filepath = final_path if final_path.exists() else fallback_path

        if not filepath.exists():
            raise FileNotFoundError(f"Neither {final_path} nor {fallback_path} exist")

        adata, _ = _load_h5ad(str(filepath))
        return adata

    except (FileNotFoundError, IOError) as e:
        print("File not found:", e)
        return None


def get_gene_choices(input_or_name) -> dict:
    """Return the pre-built gene selectize choices dict for the given dataset."""
    name = None
    try:
        if hasattr(input_or_name, 'select_dataset') and callable(input_or_name.select_dataset):
            name = input_or_name.select_dataset()
        else:
            name = input_or_name
    except Exception:
        name = input_or_name

    if not name:
        return {}

    try:
        stem = name[:-len('_final')] if name.endswith('_final') else name
        final_path = DATA_DIR / (stem + '_final.h5ad')
        fallback_path = DATA_DIR / (stem + '.h5ad')
        filepath = final_path if final_path.exists() else fallback_path

        if not filepath.exists():
            return {}

        _, choices = _load_h5ad(str(filepath))
        return choices

    except Exception:
        return {}
