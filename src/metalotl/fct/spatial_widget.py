import glasbey
import numpy as np

import plotly.graph_objects as go
from metalotl.fct.load import get_data
from metalotl._constants import GENES_LABEL, GENE_ANNOTATION

def plot_space(input):

    if not get_data(input):
        return None

    fig = go.Figure()
    fig.update_layout(
        template="plotly_dark",
        showlegend=True,
        autosize=True,
        yaxis_scaleanchor="x",
        xaxis=dict(title="x"),
        yaxis=dict(title="y")
    )
    
    return fig

def add_space_clusters(input, widget):

    if not (adata := get_data(input)):
        return None
    
    resolution = input.select_resolution()
    if not resolution or resolution not in adata.obs.columns:
        return None

    widget.data = [trace for trace in widget.data if not trace.name.isdigit() and trace.name != 'cluster']

    if input.switch_clusters():

        cluster_ids = np.unique(adata.obs[resolution].dropna())
        colors = glasbey.create_palette(palette_size=len(cluster_ids), lightness_bounds=(50, 100))
        
        # Get spatial coordinates
        x_coords = adata.obsm['spatial'][:, 0]
        y_coords = adata.obsm['spatial'][:, 1]

        for color, cluster_id in zip(colors, cluster_ids):
            mask = adata.obs[resolution] == cluster_id
            widget.add_trace(
                go.Scatter(
                    x=x_coords[mask],
                    y=y_coords[mask],
                    name=str(cluster_id),
                    mode='markers',
                    marker=dict(
                        color=color,
                        size=input.slider_dotsize_space()
                    )
                )
            )

def add_space_expression(input, widget):

    if not (adata := get_data(input)):
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
            x_coords = adata.obsm['spatial'][:, 0]
            y_coords = adata.obsm['spatial'][:, 1]
            
            # Use annotated name for display
            display_name = GENE_ANNOTATION.get(gene_id, gene_id)
            
            widget.add_trace(
                go.Scatter(
                    name = display_name,
                    x = x_coords,
                    y = y_coords,
                    mode='markers',
                    marker=dict(
                        color = gene_expression,
                        colorscale = 'Viridis',
                        showscale = True,
                        size = input.slider_dotsize_space(),
                        colorbar=dict(
                            orientation = 'h',
                            lenmode='fraction',
                            len=0.25,
                            thickness=10,
                            y = 0.05,
                            x = 0.15)
                    ),
                )
            )
