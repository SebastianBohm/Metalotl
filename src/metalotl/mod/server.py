from shiny import reactive, render
from shinywidgets import render_widget

from metalotl.fct.load import get_data as _get_data

from metalotl.fct.umap_widget import plot_umap as _plot_umap
from metalotl.fct.umap_widget import add_umap_clusters as _add_umap_clusters
from metalotl.fct.umap_widget import add_umap_expression as _add_umap_expression

from metalotl.fct.spatial_widget import plot_space as _plot_space
from metalotl.fct.spatial_widget import add_space_clusters as _add_space_clusters
from metalotl.fct.spatial_widget import add_space_expression as _add_space_expression

from metalotl.fct.expression import plot_gene_expression as _plot_gene_expression
from metalotl.fct.expression import plot_de as _plot_de

import matplotlib.pyplot as plt
plt.style.use('dark_background')
plt.rcParams.update({'font.size': 8})

def server(input, output, session):

    # Load the data
    @reactive.Calc
    def get_data():
        return _get_data(input)

    # Set up the UMAP/PCA widget
    @render_widget
    def plot_umap():
        return _plot_umap(input)
    
    @reactive.effect
    def add_umap_clusters():
        return _add_umap_clusters(input, plot_umap.widget)
    
    @reactive.effect
    def add_umap_expression():
        return _add_umap_expression(input, plot_umap.widget)

    # Set up the spatial widget
    @render_widget
    def plot_space():
        return _plot_space(input)

    @reactive.effect
    def add_clusters():
        return _add_space_clusters(input, plot_space.widget)
    
    @reactive.effect
    def add_expression():
        return _add_space_expression(input, plot_space.widget)

    # Gene expression plot
    @render.plot()
    def plot_gene_expression():
        return _plot_gene_expression(input)
    
    # Differential gene expression plot
    @render.plot()  
    def plot_de():
        return _plot_de(input)
