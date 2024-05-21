import PySimpleGUI as Sg
from gecos_gui.common_elements import theme

Sg.ChangeLookAndFeel(theme)

listmol2 = []

pad1 = ((5, 5), (5, 5))
pad2 = ((20, 20), (10, 10))

mol_grid_conf = Sg.Frame(title="Load QM Conformers",
                         layout=[
                             [Sg.Text('Opt Results Folder:', size=(22, 1), pad=pad1),
                              Sg.Input(key='-QM_PROP_MOL2LOCAL_DIR-', size=(120, 1), enable_events=True,
                                       tooltip='Folder containing the mol2 files of the alligned QM molecules.\n'
                                               'This folder is created from a previous optimization job from GeCos.\n'
                                               'The directory must exist in the local server', pad=pad1),
                              Sg.FolderBrowse(button_text="Browse", key='-QMMOL_MOL2DIR_BROWSER-'),
                              ],
                             [Sg.Text('Prop Local Folder:', size=(22, 1), pad=pad1),
                              Sg.Input(key='-QM_PROP_LOCAL_DIR-', size=(120, 1), enable_events=True,
                                       tooltip='Local Folder to store property QM calculations.'
                                               'The directory must exist in the local server', pad=pad1),
                              Sg.FolderBrowse(button_text="Browse", key='-QMMOL_DIR_BROWSER-'),
                              ],
                             [Sg.Text('Prop Remote Folder:', size=(22, 1), pad=pad1),
                              Sg.Input(key='-QM_PROP_REMOTE_DIR-', size=(120, 1), enable_events=True,
                                       tooltip='Folder where QM calculations will be stored.'
                                               'The directory must exist in the remote server', pad=pad1),
                              ],
                             [Sg.Text('Cut-off energy (kcal/mol):', size=(22, 1), pad=pad1),
                              Sg.Input(key='-QM_PROP_CUTOFF_ENERGY-', size=(10, 1), enable_events=True,
                                       tooltip='Cut-off energy to send the optimized structure\n'
                                               'to property calculation.\n'
                                               'Example: 10 kcal/mol, send all structures\n'
                                               'within 10kcal/mol from the minimimun', default_text="25.0",
                                       pad=pad1, justification='right'),
                              ],
                         ], title_color='blue', pad=pad2, size=(1350, 180))

browse_analysis = Sg.Column([
    [Sg.Frame('QM optimized conformers',
              layout=[
                  [Sg.Text('Number of files:', size=(14, 1)), Sg.Text('0', key='-NUM_MOL2FILES-')],
                  [Sg.Listbox(listmol2, key='-LISTBOX_MOL2-', size=(40, 20),
                              tooltip='Keywords for Gaussian16.', pad=pad1,
                              background_color="white")
                   ],
              ], title_color='blue', pad=pad2)
     ],

])

pad5 = ((5, 0), (0, 0))
pad3 = ((5, 20), (20, 0))
pad4 = ((15, 0), (20, 0))

browse_analysis2 = Sg.Column([
    [Sg.Frame('Keyword line',
              layout=[
                  [Sg.Multiline('#p M062X/6-311G**', size=(600, 2), pad=pad5, key='-KEYWORD_LINE-', disabled=True)],
              ], title_color='blue', pad=pad1, size=(900, 80))],
    [Sg.Frame('Properties',
              layout=[
                  [Sg.Text('Method:', size=(8, 1), pad=pad4),
                   Sg.Input(key='-INPUT_METHOD_PROP-', size=(16, 1),
                            tooltip='Theory method valid for gaussian16', enable_events=True,
                            default_text="M062X", pad=pad3, disabled=False, text_color="black"),
                   Sg.Text('Basis Set:', size=(8, 1), pad=pad3),
                   Sg.Input(key='-INPUT_BASISSET_PROP-', size=(16, 1), enable_events=True,
                            tooltip='Basis set valid for gaussian16',
                            default_text="6-311G**", pad=pad3, disabled=False, text_color="black"),
                   Sg.Checkbox('Optimize?', key='-CHECKBOX_PROP_OPT-', size=(22, 1),
                               tooltip='Reoptimize the molecule', enable_events=True,
                               default=False, pad=pad3, disabled=False, checkbox_color="white"),
                   Sg.Checkbox('Run Gaussian?', key='-CHECKBOX_RUN_GAUSSIAN_OPT_PROP-',
                               pad=pad3, enable_events=True, default=False)
                   ],
                  [
                   Sg.Text('Solvent model:', size=(12, 1), pad=pad3),
                   Sg.Combo(values=['None', 'PCM', 'IEFPCM', 'SMD', 'IPCM', 'SCIPCM', 'CPCM'],
                            key='-COMBO_MODELSOLVENT_PROP-',
                            enable_events=True, pad=pad3, size=(16, 9), default_value='None'),

                   Sg.Text('Solvent:', size=(8, 1), pad=pad4),
                   Sg.Combo(values=['None', 'Water', 'DMSO', 'Nitromethane', 'Acetonitrile', 'Methanol',
                                    'Ethanol', 'Acetone', 'Dichloromethane', 'Dichloroethane', 'THF',
                                    'Aniline', 'Chlorobenzene', 'Chloroform', 'Diethyl eter', 'Toluene',
                                    'Benzene', 'CCl4', 'Cyclohexane', 'Heptane', 'Other...'],
                            key='-COMBO_SOLVENT_PROP-',
                            enable_events=True, pad=pad3, size=(16, 21), default_value='None'),
                   Sg.Input(key='-INPUT_OTHERSOLV_PROP-', size=(16, 1),
                            tooltip='More solvents in gaussian web page',
                            pad=pad3, text_color="black", visible=False, enable_events=True),
                   Sg.Text('https://gaussian.com/scrf/', key='-TEXT_OTHERSOLVWEB_PROP-',
                           size=(24, 1), visible=False, enable_events=True),
                   ],
                  [
                   Sg.Text('Database filename:', size=(18, 1), pad=pad4),
                   Sg.Input(key='-INPUT_DATABASE_PROP-', size=(16, 1),
                            tooltip='Name of the database to follow the property calculations',
                            pad=pad3, text_color="black", enable_events=True),
                   Sg.Text('Output (log) filename:', size=(20, 1), pad=pad4),
                   Sg.Input(key='-INPUT_LOG_PROP-', size=(16, 1),
                            tooltip='Log for properties',
                            pad=pad3, text_color="black", enable_events=True),
                  ],
                  [Sg.Radio('Single point', key='-RADIO_SP-', group_id=1,
                            enable_events=True, default=True, pad=pad3)],
                  [Sg.Radio('Frequency calculations', key='-RADIO_FREQ-', group_id=1,
                            enable_events=True, default=False, pad=pad3),
                   Sg.Text('Scale value:', size=(10, 1), pad=pad3),
                   Sg.Input(key='-INPUT_SCALEVAL_FREQ-', size=(16, 1),
                            tooltip='Scale value for IR calculations', enable_events=True,
                            default_text="0.987", pad=pad3, disabled=True, text_color="gray"),
                   Sg.Checkbox('Anharmonic corrections', key='-CHECKBOX_ANHAR_FREQ-', size=(22, 1),
                               tooltip='Applied anharmocity corrections', enable_events=True,
                               default=False, pad=pad3, disabled=True, text_color="gray", checkbox_color="white")
                   ],
                  [Sg.Radio('Electrostatic potential & Wavefunction analysis', key='-RADIO_ESPWFN-',
                            enable_events=True, default=False, pad=pad3, group_id=1)
                   ],
                  [Sg.Radio('NMR', key='-RADIO_NMR-',
                            enable_events=True, default=False, pad=pad3, group_id=1)
                   ],
              ], title_color='blue', pad=pad1, size=(900, 400))
    ],
    [Sg.Button('Check/Run Properties', key='-BUTTONRUN_PROP-', disabled=True),
     Sg.Button('Import Python Properties', key='-IMPORTPYTHON_PROP-', disabled=False),
     Sg.Button('Export Python Properties', key='-EXPORTPYTHON_PROP-', disabled=True),
     ],
    [
     Sg.Input('No python for properties is loaded.', key='-PROP_HIDEINPUTSCRIPT-', visible=True,
              readonly=True, justification="right", tooltip='Name of the loaded python script for properties.',
              size=(110, 1)),
    ],
])

prop_conf_layout = [[mol_grid_conf], [browse_analysis, browse_analysis2]]
