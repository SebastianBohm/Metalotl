import os
import scanpy as sc
from functools import lru_cache

from metalotl._constants import *

# Cache for loaded datasets - stores up to 5 datasets in memory
_data_cache = {}

def _load_h5ad(filepath):
    """Load h5ad file with caching."""
    if filepath not in _data_cache:
        _data_cache[filepath] = sc.read_h5ad(filepath)
    return _data_cache[filepath]

def get_data(input):
        
    if not (name := input.select_dataset()):
        return None
    
    try:
        filepath = os.path.join(DATA_DIR, name + '.h5ad')
        adata = _load_h5ad(filepath)
        return adata
    
    except (FileNotFoundError, IOError):
        print("File not found")
        return None
