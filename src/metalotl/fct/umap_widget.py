from metalotl.fct.load import get_data
from metalotl._constants import GENES_LABEL, GENE_ANNOTATION, G2M_GENES

import numpy as np
import plotly.graph_objects as go
import glasbey

# Cache glasbey palettes by size — creation takes ~100 ms, reuse is free
_palette_cache: dict = {}

def _get_palette(n: int) -> list:
    if n not in _palette_cache:
        _palette_cache[n] = glasbey.create_palette(palette_size=n, lightness_bounds=(50, 100))
    return _palette_cache[n]


def _get_umap_or_pca(adata):
    """Return (coords, label) — UMAP if available, PCA as fallback."""
    if 'X_umap' in adata.obsm:
        return adata.obsm['X_umap'], 'UMAP'
    if 'X_pca' in adata.obsm:
        return adata.obsm['X_pca'], 'PCA'
    return None, None


def _layout_2d(axis_title_1='UMAP 1', axis_title_2='UMAP 2'):
    return dict(
        template="plotly_dark",
        showlegend=True,
        legend=dict(itemsizing='constant', tracegroupgap=0),
        autosize=True,
        xaxis=dict(showgrid=False, showticklabels=False, title=axis_title_1,
                   zeroline=False, automargin=True),
        yaxis=dict(showgrid=False, showticklabels=False, title=axis_title_2,
                   zeroline=False, automargin=True),
        margin=dict(l=40, r=40, t=20, b=40),
    )


def compute_g2m_on_demand(adata):
    """Compute G2M score lazily and store in adata.obs['g2m_score'] if not present."""
    if 'g2m_score' in adata.obs.columns:
        return True
    try:
        genes_in_data = [g for g in G2M_GENES if g in adata.var_names]
        if len(genes_in_data) == 0:
            return False
        gene_indices = [adata.var_names.get_loc(g) for g in genes_in_data]
        if hasattr(adata.X, 'toarray'):
            expr = adata.X[:, gene_indices].toarray()
        else:
            expr = adata.X[:, gene_indices]
        expr = np.asarray(expr, dtype=float)
        means = np.nanmean(expr, axis=0)
        stds = np.nanstd(expr, axis=0)
        stds[stds == 0] = 1.0
        z = (expr - means) / stds
        z = np.nan_to_num(z)
        adata.obs['g2m_score'] = z.mean(axis=1)
        return True
    except Exception:
        return False


def plot_umap(input):

    adata = get_data(input)
    if not adata:
        return None

    if 'X_umap' not in adata.obsm and 'X_pca' not in adata.obsm:
        return None

    _, label = _get_umap_or_pca(adata)
    fig = go.Figure()
    fig.update_layout(**_layout_2d(f"{label} 1", f"{label} 2"))
    return fig


def add_umap_clusters(input, widget):

    if not (adata := get_data(input)):
        return None

    resolution = input.select_resolution()
    if not resolution or resolution not in adata.obs.columns:
        return None

    coords, label = _get_umap_or_pca(adata)
    if coords is None:
        return None

    widget.data = [trace for trace in widget.data if getattr(trace, 'legendgroup', None) != 'cluster']

    # Set axis titles and explicit ranges from data so the view always fits
    pad = 0.05
    x_min, x_max = float(coords[:, 0].min()), float(coords[:, 0].max())
    y_min, y_max = float(coords[:, 1].min()), float(coords[:, 1].max())
    x_pad = (x_max - x_min) * pad
    y_pad = (y_max - y_min) * pad
    widget.update_layout(
        xaxis=dict(title=f"{label} 1", range=[x_min - x_pad, x_max + x_pad],
                   showgrid=False, showticklabels=False, zeroline=False, automargin=True),
        yaxis=dict(title=f"{label} 2", range=[y_min - y_pad, y_max + y_pad],
                   showgrid=False, showticklabels=False, zeroline=False, automargin=True),
    )

    if input.switch_clusters():

        cluster_ids = np.unique(adata.obs[resolution].dropna())
        colors = _get_palette(len(cluster_ids))

        for color, cluster_id in zip(colors, cluster_ids):
            mask = adata.obs[resolution] == cluster_id
            widget.add_trace(
                go.Scatter(
                    name=str(cluster_id),
                    legendgroup='cluster',
                    x=coords[mask, 0],
                    y=coords[mask, 1],
                    mode='markers',
                    marker=dict(color=color, size=input.slider_dotsize_umap())
                )
            )


def add_umap_expression(input, widget):

    if not (adata := get_data(input)):
        return None

    coords, label = _get_umap_or_pca(adata)
    if coords is None:
        return None

    # Remove existing expression / G2M traces (GENES_LABEL is a set of all annotated names)
    widget.data = [trace for trace in widget.data
                   if trace.name not in GENES_LABEL
                   and trace.name != 'G2M Score']

    if input.switch_g2m():
        try:
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
                            size=input.slider_dotsize_umap(),
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
                    size=input.slider_dotsize_umap(),
                    colorbar=dict(orientation='v', thickness=15,
                                  xanchor='left', x=1.02),
                )
            )
        )
    except Exception:
        return
