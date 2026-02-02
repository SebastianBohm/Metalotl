from shiny import App
from metalotl.mod.ui import app_ui as _ui
from metalotl.mod.server import server as _server

app = App(_ui, _server)
