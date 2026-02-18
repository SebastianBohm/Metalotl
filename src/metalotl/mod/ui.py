from shiny import ui
from shinywidgets import output_widget

from metalotl._constants import DATA, CLUSTERING_OPTIONS
from metalotl.js._format import DROPDOWN_CONFIG
from metalotl import __version__

app_ui = ui.page_navbar(

    # Spatial Data
    ui.nav_panel(ui.HTML("Spatial<br><span style='font-size: smaller;'>data v1.0.0</span>"),

        # Sidebar
        ui.page_sidebar(
            ui.sidebar(
                ui.input_selectize(
                    "select_dataset",
                    "Select dataset",
                    ['', *DATA],
                    selected=None,
                    options={
                        "render": DROPDOWN_CONFIG
                        }
                    ),
                ui.input_selectize(
                    "select_resolution",
                    "Cluster resolution",
                    choices={},
                    selected=None,
                    ),
                ui.input_switch("switch_clusters", "Show clusters", True),
                ui.input_selectize(
                    "select_gene",
                    "Select gene",
                    choices={},
                    selected=None,
                    options={
                        "render": DROPDOWN_CONFIG
                        }
                    ),
                ui.input_switch("switch_expression", "Plot gene expression", False),
                ui.input_switch("switch_g2m", "Show G2M score", False),
                ui.input_slider("slider_dotsize_umap", "Dot size UMAP", 1, 20, value=2),
                ui.input_slider("slider_dotsize_space", "Dot size Space", 1, 20, value=2),
            ),
            ui.layout_columns(
                ui.card(
                    ui.card_header("UMAP Projection"),
                    output_widget('plot_umap', height="400px"),
                    output_widget('plot_umap_extra', height="400px"),
                    full_screen=True),
                ui.card(
                    ui.card_header("Spatial plot"),
                    output_widget("plot_space", height="400px"),
                    output_widget("plot_space_extra", height="400px"),
                    full_screen=True),
                col_widths={"sm": (5, 7, 12)}
                )

            )
    ),

    # Other components of the header
    ui.nav_spacer(),
    ui.nav_control(ui.input_dark_mode(id="mode", mode='dark')),

    title=ui.HTML(f"Metalotl - Metamorphosed Axolotl Spatial Data<br><span style='font-size: 12px; display: block; '>app v{__version__}</span>"),
    id="page"

)
