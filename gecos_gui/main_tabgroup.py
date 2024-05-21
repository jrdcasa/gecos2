import PySimpleGUI as Sg
from gecos_gui.general_inputs_tab import general_inputs_layout
from gecos_gui.prop_conf_tab import prop_conf_layout
from gecos_gui.common_elements import theme

Sg.ChangeLookAndFeel(theme)

tabgroup_layout = Sg.TabGroup([
                   [Sg.Tab('Conformer optimization', general_inputs_layout)],
                   [Sg.Tab('Conformer properties', prop_conf_layout)]
                 ], key="-TABGROUP-")

tabgroup_layout_col = Sg.Column([
                            [tabgroup_layout]
                        ])
