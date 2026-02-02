from metalotl.fct.load import get_data
from metalotl._constants import GENES_LABEL, GENE_ANNOTATION

import numpy as np
import plotly.graph_objects as go
import glasbey

def plot_umap(input):

    adata = get_data(input)
    if not adata:
        return None

    # Check if UMAP or PCA coordinates exist
    has_umap = 'X_umap' in adata.obsm
    has_pca = 'X_pca' in adata.obsm
    
    if not has_umap and not has_pca:
        return None

    fig = go.Figure()
    fig.update_layout(
        template="plotly_dark",
        showlegend=False,
        autosize=True,
        scene=dict(
            xaxis=dict(showgrid=False, showticklabels=False, title='', zeroline=False),
            yaxis=dict(showgrid=False, showticklabels=False, title='', zeroline=False),
            zaxis=dict(showgrid=False, showticklabels=False, title='', zeroline=False))
    )

    return fig

def add_umap_clusters(input, widget):
    
    if not (adata := get_data(input)):
        return None
    
    resolution = input.select_resolution()
    if not resolution or resolution not in adata.obs.columns:
        return None

    # Check if we have UMAP or PCA
    if 'X_umap' in adata.obsm:
        coords = adata.obsm['X_umap']
    elif 'X_pca' in adata.obsm:
        coords = adata.obsm['X_pca']
    else:
        return None

    widget.data = [trace for trace in widget.data if not trace.name.isdigit()]

    if input.switch_clusters():
        
        cluster_ids = np.unique(adata.obs[resolution].dropna())
        colors = glasbey.create_palette(palette_size=len(cluster_ids), lightness_bounds=(50, 100))

        for color, cluster_id in zip(colors, cluster_ids):
            mask = adata.obs[resolution] == cluster_id

            x_coords = coords[mask, 0]
            y_coords = coords[mask, 1]
            z_coords = coords[mask, 2] if coords.shape[1] > 2 else np.zeros(mask.sum())

            widget.add_trace(
                go.Scatter3d(
                name = str(cluster_id),
                x = x_coords,
                y = y_coords,
                z = z_coords,
                mode='markers',
                marker=dict(
                    color = color,
                    size = input.slider_dotsize_umap()
                )
                )
            )

def add_umap_expression(input, widget):

    if not (adata := get_data(input)):
        return None

    # Check if we have UMAP or PCA
    if 'X_umap' in adata.obsm:
        coords = adata.obsm['X_umap']
    elif 'X_pca' in adata.obsm:
        coords = adata.obsm['X_pca']
    else:
        return None

    # Remove existing gene expression traces
    widget.data = [trace for trace in widget.data if trace.name not in GENES_LABEL and trace.name not in GENE_ANNOTATION.values()]

    if input.switch_expression() and input.select_gene():

        gene_id = input.select_gene()  # This is the AMEX ID
        if gene_id in adata.var_names:
            # Handle sparse matrices by converting to dense array
            gene_data = adata[:,gene_id].X
            if hasattr(gene_data, 'toarray'):
                gene_expression = np.array(gene_data.toarray()).flatten()
            else:
                gene_expression = np.array(gene_data).flatten()
            x_coords = coords[:, 0]
            y_coords = coords[:, 1]
            z_coords = coords[:, 2] if coords.shape[1] > 2 else np.zeros(len(x_coords))
            
            # Use annotated name for display
            display_name = GENE_ANNOTATION.get(gene_id, gene_id)
            
            widget.add_trace(
                go.Scatter3d(
                    name = display_name,
                    x = x_coords,
                    y = y_coords,
                    z = z_coords,
                    mode='markers',
                    marker=dict(
                        color = gene_expression,
                        colorscale = 'Viridis',
                        showscale = True,
                        size = input.slider_dotsize_umap(),
                        colorbar=dict(
                            orientation = 'h',
                            lenmode='fraction',
                            len=0.25,
                            thickness=10,
                            y = 0.05,
                            x = 0.15)
                    )
                )
            )
