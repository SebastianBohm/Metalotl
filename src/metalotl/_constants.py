from pathlib import Path
import numpy as np
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent.parent
DATA_DIR = BASE_DIR / 'data'

# Clustering columns present in the datasets, in preference order
CLUSTERING_OPTIONS = {
    'spatial_leiden_e30_s8': 'Leiden clustering',
    'structure': 'Structure annotation',
    'seurat_clusters': 'Seurat clusters',
}

# Discover datasets from *_final.h5ad files in data/
try:
    DATA = sorted([p.stem.replace('_final', '') for p in DATA_DIR.glob('*_final.h5ad')])
    if not DATA:
        DATA = np.load(DATA_DIR / 'samples.npy', allow_pickle=True).tolist()
except Exception:
    DATA = []

# Gene annotation lookup: AMEX gene_id -> annotated gene name
ANNOTATION_FILE = DATA_DIR / 'Adult_meta_DGE_markers.csv'
try:
    _ann_df = pd.read_csv(ANNOTATION_FILE, index_col=0, low_memory=False)
    _ann_df = _ann_df[~_ann_df.index.duplicated(keep='first')]
    GENE_ANNOTATION = {}
    for gene_id, row in _ann_df.iterrows():
        name = row.get('Axolotl_tanaka_annotated_gene', None)
        if pd.isna(name) or name in ('-', 'N/A') or name is None:
            name = gene_id
        GENE_ANNOTATION[str(gene_id)] = str(name)
    del _ann_df
except Exception:
    GENE_ANNOTATION = {}

# Set of annotated display names — used to identify expression traces for removal
GENES_LABEL = set(GENE_ANNOTATION.values())

# G2M marker genes used to compute a proliferation score
G2M_GENES = [
    'MKI67', 'TOP2A', 'CCNB1', 'CCNB2', 'CCNA2', 'CDC20', 'PLK1', 'AURKA', 'AURKB',
    'BUB1', 'BUB1B', 'CENPA', 'CENPE', 'CENPF', 'KIF11', 'KIF20A', 'KIF23', 'NUSAP1',
    'TPX2', 'UBE2C', 'BIRC5', 'NUF2', 'PRC1', 'SMC1A', 'NEK2', 'DLGAP5', 'H2AZ1',
    'PTTG1', 'CDK1', 'HMGB2'
]
