# Metalotl 🦎

[![Install](https://github.com/quadbio/metalotl/actions/workflows/install.yml/badge.svg)](https://github.com/quadbio/metalotl/actions/workflows/install.yml)
[![Python 3.12+](https://img.shields.io/badge/python-3.12+-blue.svg)](https://www.python.org/downloads/)
[![Shiny](https://img.shields.io/badge/shiny-python-green.svg)](https://shiny.posit.co/py/)

An interactive Shiny application for exploring spatial transcriptomics data from the metamorphosed axolotl brain.

## 📋 Overview

Metalotl provides an interactive web interface to explore spatial gene expression data from metamorphosed axolotl brain regions, including:

| Brain Region | Replicates |
|--------------|------------|
| 🧠 Metencephalon (hindbrain) | 2 |
| 👃 Olfactory bulb | 2 |
| 🔴 Pituitary | 1 |
| 🧩 Telencephalon (forebrain) | 1 |
| 🔷 Thalamencephalon (diencephalon) | 3 |

## ✨ Features

- **Spatial visualization**: View gene expression patterns in their spatial context
- **PCA projection**: Explore cell clusters in 3D reduced dimensional space
- **Gene expression plots**: Visualize expression levels across cell clusters
- **Differential expression**: Identify marker genes for each cluster
- **Annotated gene names**: Gene IDs mapped to Axolotl Tanaka annotations (~8,200 genes)
- **Fast loading**: Cached data for responsive interactions after initial load

---

## 🚀 Installation

### Requirements
- Python 3.12 or higher
- A package manager: [conda](https://docs.conda.io/) or [mamba](https://mamba.readthedocs.io/) (recommended)

### Setup Instructions

**a. On a remote server:** Connect via SSH  
**b. On a local machine (MacOS):** Open Terminal  
**c. On a local machine (Windows):** Press `Windows key` + `X`, select Windows Terminal

```bash
# Create a new conda environment
mamba create -n metalotl python=3.12

# Activate the environment
mamba activate metalotl

# Clone the repository
git clone --branch main https://github.com/quadbio/metalotl.git

# Navigate to the directory
cd metalotl

# Install the package
pip install -e .
```

---

## 📊 Data Setup

### Required Files

Place the `.h5ad` data files in the `data/` directory:

```
data/
├── Meta_metencephalon_rep1_DP8400015649BR_C1-2_region_ann.h5ad
├── Meta_metencephalon_rep2_DP8400015649BR_C1-1_region_ann.h5ad
├── Meta_olfactory_bulb_rep1_DP8400015234BL_B2-1_region_ann.h5ad
├── Meta_olfactory_bulb_rep2_DP8400015234BL_B3-1_region_ann.h5ad
├── Meta_pituitary_rep5_DP8400015234BL_B3-2_region_ann.h5ad
├── Meta_telencephalon_rep3_DP8400015234BL_B5-1_region_ann.h5ad
├── Meta_thalamencephalon_rep1_DP8400015234BL_B4-1_region_ann.h5ad
├── Meta_thalamencephalon_rep3_DP8400015234BL_B5-2_region_ann.h5ad
├── Meta_thalamencephalon_rep5_DP8400015234BL_B3-2_region_ann.h5ad
├── genes.npy
└── samples.npy
```

### Gene Annotations

Gene annotations are automatically loaded from `../Result/Adult_meta_DGE_markers.csv` (relative to the Metalotl directory).

---

## 🖥️ Running the App

### On a Local Machine

```bash
# Activate the environment
mamba activate metalotl

# Navigate to the project directory
cd metalotl

# Run the Shiny app
shiny run src/metalotl/app.py
```

Open your browser: **http://localhost:8000**

### On a Remote Server

1. Connect with port forwarding:
```bash
ssh -L 12345:localhost:8000 username@server
```

2. On the server, run:
```bash
mamba activate metalotl
cd metalotl
shiny run src/metalotl/app.py --port 8000
```

3. Access locally at: **http://localhost:12345**

---

## 🎮 Usage Guide

1. **Select a dataset** from the dropdown menu
2. **Choose clustering** (Leiden clustering, Structure annotation, or Seurat clusters)
3. **Toggle cluster visualization** with the "Show clusters" switch
4. **Search for a gene** using the annotated gene names (e.g., "GLUL", "GAD1")
5. **Enable expression plotting** with the "Plot gene expression" switch
6. **Adjust visualization** using the dot size sliders
7. **Explore markers** in the differential expression accordion panel

---

## 📁 Project Structure

```
Metalotl/
├── data/                       # H5AD data files
├── src/metalotl/
│   ├── __init__.py
│   ├── _constants.py           # Configuration & gene annotations
│   ├── app.py                  # Main Shiny app entry point
│   ├── fct/
│   │   ├── expression.py       # Gene expression plotting
│   │   ├── load.py             # Data loading with caching
│   │   ├── spatial_widget.py   # Spatial plot functions
│   │   └── umap_widget.py      # PCA/UMAP plot functions
│   ├── js/
│   │   └── _format.py          # Dropdown formatting
│   └── mod/
│       ├── server.py           # Shiny server logic
│       └── ui.py               # Shiny UI definition
├── scripts/
│   └── create_tarball.sh       # Data packaging script
├── setup.py
├── pyproject.toml
└── README.md
```

---

## 🔧 Dependencies

| Package | Purpose |
|---------|---------|
| [Shiny for Python](https://shiny.posit.co/py/) | Web application framework |
| [Scanpy](https://scanpy.readthedocs.io/) | Single-cell analysis |
| [Plotly](https://plotly.com/python/) | Interactive visualizations |
| [Glasbey](https://github.com/lmcinnes/glasbey) | Color palette generation |
| [Pandas](https://pandas.pydata.org/) | Data manipulation |
| [NumPy](https://numpy.org/) | Numerical computing |

---

## 🙏 Acknowledgments

- **Adnan** for the template

---
