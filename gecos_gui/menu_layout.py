import PySimpleGUI as sg
from gecos_gui.common_elements import theme

sg.ChangeLookAndFeel(theme)

# ------ Menu Definition ------ #
# menu_def = [['File', ['Import Json...', 'Export Json...', 'Import Python Script...', 'Write Python Script...',
# 'Exit']],
#             ['Edit', ['Delete database and done file', 'Clean Form', ], ],
#             ['Help', 'About...'], ]
menu_def = [['File', ['Import Python Script...', 'Write Python Script...', 'Exit']],
            ['Edit', ['Delete database and done file', 'Clean Form', ], ],
            ['Help', 'About...'], ]

menu_layout = sg.Menu(menu_def, tearoff=False, pad=(200, 1),
                      font=("Helvetica", 14), key="-MENU-")
