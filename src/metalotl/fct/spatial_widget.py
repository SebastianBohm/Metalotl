import glasbey
import numpy as np

import plotly.graph_objects as go
from metalotl.fct.load import get_data
from metalotl._constants import GENES_LABEL, GENE_ANNOTATION

# Cache glasbey palettes by size — creation takes ~100 ms, reuse is free
_palette_cache: dict = {}

def _get_palette(n: int) -> list:
    if n not in _palette_cache:
        _palette_cache[n] = glasbey.create_palette(palette_size=n, lightness_bounds=(50, 100))
    return _palette_cache[n]

def plot_space(input):

    if not get_data(input):
        return None

    fig = go.Figure()
    fig.update_layout(
        template="plotly_dark",
        showlegend=True,
        autosize=True,
        yaxis_scaleanchor="x",
        xaxis=dict(title="x", automargin=True),
        yaxis=dict(title="y", autorange="reversed", automargin=True),
        margin=dict(l=40, r=40, t=20, b=40),
    )
    return fig

def add_space_clusters(input, widget):

    if not (adata := get_data(input)):
        return None

    resolution = input.select_resolution()
    if not resolution or resolution not in adata.obs.columns:
        return None

    widget.data = [trace for trace in widget.data if getattr(trace, 'legendgroup', None) != 'cluster']

    if 'spatial' not in adata.obsm:
        return None

    # Get spatial coordinates
    x_coords = adata.obsm['spatial'][:, 0]
    y_coords = adata.obsm['spatial'][:, 1]

    # Set explicit axis ranges from data so the view always fits correctly
    pad = 0.02
    x_min, x_max = float(x_coords.min()), float(x_coords.max())
    y_min, y_max = float(y_coords.min()), float(y_coords.max())
    x_pad = (x_max - x_min) * pad
    y_pad = (y_max - y_min) * pad
    widget.update_layout(
        xaxis=dict(title="x", range=[x_min - x_pad, x_max + x_pad],
                   automargin=True),
        yaxis=dict(title="y", range=[y_max + y_pad, y_min - y_pad],  # reversed
                   automargin=True, autorange=False),
    )

    if input.switch_clusters():
        cluster_ids = np.unique(adata.obs[resolution].dropna())
        colors = _get_palette(len(cluster_ids))

        for color, cluster_id in zip(colors, cluster_ids):
            mask = adata.obs[resolution] == cluster_id
            widget.add_trace(
                go.Scatter(
                    x=x_coords[mask],
                    y=y_coords[mask],
                    name=str(cluster_id),
                    legendgroup='cluster',
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

    # Remove existing expression / G2M traces (GENES_LABEL is a set of all annotated names)
    widget.data = [trace for trace in widget.data
                   if trace.name not in GENES_LABEL
                   and trace.name != 'G2M Score']

    if not (input.switch_expression() or input.switch_g2m()):
        return

    if 'spatial' in adata.obsm:
        coords = adata.obsm['spatial']
    elif 'X_spatial' in adata.obsm:
        coords = adata.obsm['X_spatial']
    else:
        return

    if input.switch_g2m():
        try:
            from metalotl.fct.umap_widget import compute_g2m_on_demand
            compute_g2m_on_demand(adata)
        except Exception:
            pass

        if 'g2m_score' in adata.obs.columns:
            try:
                score = np.asarray(adata.obs['g2m_score'].values).flatten()
                widget.add_trace(
                    go.Scatter(
                        name='G2M Score',
                        x=coords[:, 0],
                        y=coords[:, 1],
                        mode='markers',
                        marker=dict(
                            color=score,
                            colorscale='RdBu_r',
                            cmin=-np.percentile(np.abs(score), 99),
                            cmax=np.percentile(np.abs(score), 99),
                            showscale=True,
                            size=input.slider_dotsize_space(),
                            colorbar=dict(orientation='v', thickness=15,
                                          xanchor='left', x=1.02)
                        )
                    )
                )
            except Exception:
                pass

    if not input.switch_expression():
        return

    try:
        gene_id = input.select_gene() or None
    except Exception:
        gene_id = None

    if not gene_id or gene_id not in adata.var_names:
        return

    try:
        gene_data = adata[:, gene_id].X
        if hasattr(gene_data, 'toarray'):
            gene_expression = np.array(gene_data.toarray()).flatten()
        else:
            gene_expression = np.array(gene_data).flatten()

        display_name = GENE_ANNOTATION.get(gene_id, gene_id)
        widget.add_trace(
            go.Scatter(
                name=display_name,
                x=coords[:, 0],
                y=coords[:, 1],
                mode='markers',
                marker=dict(
                    color=gene_expression,
                    colorscale='Viridis',
                    showscale=True,
                    size=input.slider_dotsize_space(),
                    colorbar=dict(orientation='v', thickness=15,
                                  xanchor='left', x=1.02),
                ),
            )
        )
    except Exception:
        return
