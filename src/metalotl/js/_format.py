from shiny import ui

DROPDOWN_CONFIG = ui.js_eval(
    '{option: function(item, escape) {return "<div style=\\"font-size: 12px; padding: 5px; display: flex; justify-content: space-between; align-items: center;\\">" + escape(item.label) + "<span style=\\"flex-grow: 1;\\"></span></div>";},'
    + 'item: function(item, escape) {return "<div style=\\"font-size: 12px; display: flex; justify-content: space-between; align-items: center;\\">" + escape(item.label) + "<span style=\\"flex-grow: 1;\\"></span></div>";}}'
)
