from metalotl.fct.load import get_data
from metalotl._constants import GENE_ANNOTATION
import matplotlib.pyplot as plt
import scanpy as sc

def plot_gene_expression(input):  

    adata = get_data(input)
    plot_genes = list(input.select_gene_expression())
    resolution = input.select_resolution()

    if adata is None or not plot_genes:
        return
    
    if not resolution or resolution not in adata.obs.columns:
        return

    plot_ex = sc.pl.dotplot(adata,
                            var_names = plot_genes,
                            swap_axes = True,
                            groupby = resolution,
                            return_fig = True
                            )

    main_ax_ex = plot_ex.get_axes()['mainplot_ax']
    # Convert AMEX IDs to annotated gene names
    y_axis_labels = [tick.get_text() for tick in main_ax_ex.get_yticklabels()]
    annotated_labels = [GENE_ANNOTATION.get(label, label) for label in y_axis_labels]
    main_ax_ex.set_yticklabels(annotated_labels)

    plt.subplots_adjust(top=0.95, bottom=0.1, left=0.1, right=0.95)

def plot_de(input):  

    if not (adata := get_data(input)):
        return None
    
    resolution = input.select_resolution()
    if not resolution or resolution not in adata.obs.columns:
        return None
    
    key = 'rank_' + resolution
    if key not in adata.uns:
        # Run differential expression if not already computed
        sc.tl.rank_genes_groups(adata, groupby=resolution, key_added=key)
    
    plot_de = sc.pl.rank_genes_groups_dotplot(
        adata,
        n_genes = input.slider_n_genes(),
        key = key,
        min_logfoldchange=input.slider_lfc(),
        return_fig = True)

    main_ax_de = plot_de.get_axes()['mainplot_ax']
    # Convert AMEX IDs to annotated gene names
    x_axis_labels = [tick.get_text() for tick in main_ax_de.get_xticklabels()]
    annotated_labels = [GENE_ANNOTATION.get(label, label) for label in x_axis_labels]
    main_ax_de.set_xticklabels(annotated_labels)
    
    plt.subplots_adjust(top=1, bottom=0.2, left=0.05, right=0.95)
