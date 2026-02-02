from pathlib import Path
import numpy as np
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent.parent
DATA_DIR = BASE_DIR / 'data'

# Possible clustering options - will be filtered based on what's available in each dataset
CLUSTERING_OPTIONS = {
    'spatial_leiden_e30_s8': 'Leiden clustering',
    'structure': 'Structure annotation',
    'seurat_clusters': 'Seurat clusters'
}

DATA = np.load(DATA_DIR / 'samples.npy', allow_pickle=True).tolist()

# Load annotation metadata for differential gene expression markers
ANNOTATION_FILE = BASE_DIR.parent / 'Result' / 'Adult_meta_DGE_markers.csv'
try:
    ANNOTATION_DF = pd.read_csv(ANNOTATION_FILE, index_col=0, low_memory=False)
    
    # Create gene annotation lookup (AMEX gene_id -> annotated gene name)
    GENE_ANNOTATION = {}
    for gene_id in ANNOTATION_DF.index.unique():
        row = ANNOTATION_DF.loc[gene_id]
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        # Get Axolotl_tanaka_annotated_gene
        gene_name = row.get('Axolotl_tanaka_annotated_gene', None)
        if pd.isna(gene_name) or gene_name == '-' or gene_name == 'N/A' or gene_name is None:
            gene_name = gene_id  # fallback to AMEX ID
        GENE_ANNOTATION[gene_id] = gene_name

    # Get unique cell types from the DGE markers
    CELLTYPES = sorted(ANNOTATION_DF['Celltype'].unique().tolist())
except FileNotFoundError:
    ANNOTATION_DF = None
    GENE_ANNOTATION = {}
    CELLTYPES = []

# Load genes from npy and create display labels using annotations
GENES_AMEX = np.load(DATA_DIR / 'genes.npy', allow_pickle=True).tolist()

# Create gene choices for dropdown: display annotated name but use AMEX ID as value
# Format: {amex_id: annotated_name} for genes that have annotations
GENES_DISPLAY = {}
for amex_id in GENES_AMEX:
    if amex_id in GENE_ANNOTATION:
        annotated = GENE_ANNOTATION[amex_id]
        # Show "AnnotatedName (AMEX_ID)" for clarity
        GENES_DISPLAY[amex_id] = f"{annotated} ({amex_id})" if annotated != amex_id else amex_id
    else:
        GENES_DISPLAY[amex_id] = amex_id

# For backwards compatibility
GENES = GENES_AMEX
GENES_LABEL = [GENE_ANNOTATION.get(g, g) for g in GENES_AMEX]
