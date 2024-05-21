import PySimpleGUI as Sg
from gecos_gui.common_elements import theme


Sg.ChangeLookAndFeel(theme)

pad1 = ((5, 5), (5, 5))

molecule_input = Sg.Frame(title="Molecule Input",
                          layout=[
                              [Sg.Text('Molecule file:', size=(12, 1)),
                               Sg.Input(key='-MOLECULE_INPUT_EXTRACT-', size=(120, 1),
                                        tooltip='Allowed files: pdb, mol2, sdf.', enable_events=True),
                               Sg.FileBrowse(button_text="Browse", key="-MOLECULE_INPUT_EXTRACT_BROWSER-",
                                             file_types=(('ALL Files', '*'),),)]
                          ], title_color='blue', pad=((20, 20), (10, 10)))

method_extract_input = Sg.Frame(title="Method to extract molecule:",
                                layout=[
                                       [Sg.Text('Method:', size=(12, 1)),
                                        Sg.Combo(['Sphere', 'Voronoi'], disabled=False, key='-METHOD_EXTRACT-',
                                        default_value='Sphere',
                                        tooltip="Method to calculate the negihbors",
                                        size=(15, 1), pad=pad1, enable_events=True),
                                        Sg.Text('Radius (angstroms):', size=(18, 1)),
                                        Sg.Input(key='-RADIUS_SPHERE-', size=(20, 1), enable_events=True,
                                                 tooltip='Sphere method needs the radius sphere.',
                                                 pad=pad1, default_text=7.0),
                                ]
                                ], title_color='blue', pad=((20, 20), (10, 10)))

method_buttons = Sg.Frame(title= "Run extract methods:",
                          layout=[
                                   [Sg.Button('Run Extract', key='-BUTTONRUN_EXTRACT-', disabled=True),
                                    Sg.Button('Import Python Extract', key='-BUTTON_IMPORTPYTHON_EXTRACT-',
                                              disabled=False),
                                    Sg.Button('Export Python Extract', key='-BUTTON_EXPORTPYTHON_EXTRACT-',
                                              disabled=True),
                                    Sg.Input('No python to extract is loaded.', key='-EXTRACT_HIDEINPUTSCRIPT-',
                                             visible=True, readonly=True, justification="left"),
                                   ],
                          ])

extract_mol_layout = [[molecule_input], [method_extract_input], [method_buttons]]
