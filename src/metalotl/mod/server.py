from shiny import reactive, render, ui
from shinywidgets import render_widget

from metalotl.fct.load import get_data as _get_data, get_gene_choices as _get_gene_choices

from metalotl.fct.umap_widget import plot_umap as _plot_umap
from metalotl.fct.umap_widget import add_umap_clusters as _add_umap_clusters
from metalotl.fct.umap_widget import add_umap_expression as _add_umap_expression

from metalotl.fct.spatial_widget import plot_space as _plot_space
from metalotl.fct.spatial_widget import add_space_clusters as _add_space_clusters
from metalotl.fct.spatial_widget import add_space_expression as _add_space_expression


def server(input, output, session):

    # Load the data
    @reactive.Calc
    def get_data():
        return _get_data(input)

    # Update gene selector to only show genes present in the selected dataset
    @reactive.effect
    def update_gene_choices():
        choices = _get_gene_choices(input)
        ui.update_selectize("select_gene", choices=choices, selected=None, session=session)

    # Set up the UMAP/PCA widget (primary) - show clusters only
    @render_widget
    def plot_umap():
        return _plot_umap(input)

    @reactive.effect
    def add_umap_clusters():
        return _add_umap_clusters(input, plot_umap.widget)

    # Extra UMAP widget for expression/G2M visualizations
    @render_widget
    def plot_umap_extra():
        return _plot_umap(input)

    @reactive.effect
    def update_umap_extra():
        input.select_dataset()
        input.select_gene()
        if input.switch_expression() or input.switch_g2m():
            return _add_umap_expression(input, plot_umap_extra.widget)
        plot_umap_extra.widget.data = []

    # Set up the spatial widget (primary) - show clusters only
    @render_widget
    def plot_space():
        return _plot_space(input)

    @reactive.effect
    def add_clusters():
        return _add_space_clusters(input, plot_space.widget)

    # Extra spatial widget for expression plotting below primary spatial
    @render_widget
    def plot_space_extra():
        return _plot_space(input)

    @reactive.effect
    def update_space_extra():
        input.select_dataset()
        input.select_gene()
        if input.switch_expression() or input.switch_g2m():
            return _add_space_expression(input, plot_space_extra.widget)
        plot_space_extra.widget.data = []
