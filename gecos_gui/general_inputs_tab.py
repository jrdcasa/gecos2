import PySimpleGUI as Sg
from gecos_gui.common_elements import theme, row_buttons


Sg.ChangeLookAndFeel(theme)

molecule_input = Sg.Frame(title="Molecule Input",
                          layout=[
                              [Sg.Text('Molecule file:', size=(12, 1)),
                               Sg.Input(key='-MOLECULE_INPUT-', size=(142, 1),
                                        tooltip='Allowed files: pdb, mol2, sdf.', enable_events=True),
                               Sg.FileBrowse(button_text="Browse", key="-MOLECULE_INPUT_BROWSER-",
                                             file_types=(('ALL Files', '*'),),)]
                          ], title_color='blue', pad=((10, 0), (10, 0)))

pad1 = ((5, 5), (5, 0))
pad2 = ((5, 5), (5, 5))
pad3 = ((0, 10), (5, 5))

server_input = Sg.Frame(title="Server Options",
                        layout=[
                            [Sg.Text('Name server:', size=(18, 1), pad=pad1),
                             Sg.Input(key='-NAME_SERVER-', size=(20, 1), enable_events=True,
                                      tooltip='Name of the remote server to perform QM calculations.', pad=pad1),
                             Sg.Text('Username:', size=(9, 1), pad=pad1),
                             Sg.Input(key='-USER_NAME-', size=(20, 1), enable_events=True,
                                      tooltip='Username to connect in the remote server.', pad=pad1),
                             Sg.Text('Key SSH file:', size=(14, 1), pad=pad2),
                             Sg.Input(key='-KEY_SSH_FILE-', size=(52, 1), enable_events=True,
                                      tooltip='Private key for SSH connection.', pad=pad1),
                             Sg.FileBrowse(button_text="Browse", key='-KEY_SSH_BROWSER-', pad=pad1,
                                           file_types=(('ALL Files', '*'),)),
                             ],
                            [
                              Sg.Text('Encrypted passwd file:', size=(19, 1), pad=pad1),
                              Sg.Input(key='-ENCRYPT_PASS_FILE-', size=(52, 1), enable_events=True,
                                       tooltip='Encrypted file containing the pass for a SSH connection.',
                                       pad=pad1, disabled=True),

                              Sg.Text('SLURM partition:', size=(14, 1), pad=pad1),
                              Sg.Input(key='-SLURM_PART-', size=(16, 1), enable_events=True,
                                       tooltip='SLURM partition to perform QM calculations.', pad=pad1),
                              Sg.Text('SLURM partition master:', size=(21, 1), pad=pad1),
                              Sg.Input(key='-SLURM_PART_MASTER-', size=(22, 1), enable_events=True,
                                       tooltip='SLURM partition to run the deamon script.', pad=pad1),
                             ],
                            [
                             Sg.Text('Exclude Nodes:', size=(14, 1), pad=pad1),
                             Sg.Input(key='-EXCLUDE_NODES-', size=(30, 1), enable_events=True,
                                      tooltip='Nodes excluded in the partition.', pad=pad3),
                             Sg.Text('Time SLURM:', size=(12, 1), pad=pad1),
                             Sg.Input(key='-TIME_NODES-', size=(15, 1), enable_events=True,
                                      tooltip='Set a limit on the total run time of the job allocation.', pad=pad3),
                             Sg.Text('Max. number of slurm jobs:', size=(25, 1), pad=pad1),
                             Sg.Input(key='-MAX_JOBS_SLURM-', size=(10, 1), enable_events=True,
                                      tooltip='Maximum number of jobs to run at the same time in the remote host. '
                                              'A value of 50 is recommend',
                                      pad=pad3),

                             Sg.Text('Node Master:', size=(12, 1), pad=pad1),
                             Sg.Input(key='-NODE_MASTER-', size=(23, 1), enable_events=True,
                                      tooltip='Node to run the daemon script.', pad=pad1),

                             ],
                          ], title_color='blue', pad=((10, 0), (10, 0)))

directories = Sg.Column([
                 [Sg.Frame(title="Directories",
                       layout=[
                           [Sg.Text('Local dir:', size=(10, 1), pad=pad1),
                            Sg.Input(key='-LOCAL_DIR-', size=(142, 1), enable_events=True,
                                     tooltip='Local directory to put the outputs. '
                                             'The directory must exist in the local server', pad=pad1),
                            Sg.FolderBrowse(button_text="Browse", key='-LOCAL_DIR_BROWSER-'),
                            ],
                           [Sg.Text('Remote dir:', size=(10, 1), pad=pad1),
                            Sg.Input(key='-REMOTE_DIR-', size=(142, 1), enable_events=True,
                                     tooltip='Remote directory where QM calculation carry out. '
                                             'The directory must exist in the remote server', pad=pad1),
                            ],
                           [Sg.Text('Pattern:', size=(10, 1), pad=pad1),
                            Sg.Input(key='-PATTERN-', size=(35, 1), enable_events=True,
                                     tooltip='Pattern used to give the name to the generated files.', pad=pad1),
                            Sg.Text('Database Name:', size=(15, 1), pad=pad1),
                            Sg.Input(key='-DATABASE_NAME-', size=(35, 1), enable_events=True,
                                     tooltip='Name of the database.', pad=pad1),
                            Sg.Text('Filename Log:', size=(15, 1), pad=pad1),
                            Sg.Input(key='-FILENAME_LOG-', size=(35, 1),
                                     tooltip='Name of the database.', pad=pad1, enable_events=True),
                            ],
                          ], title_color='blue', pad=((10, 0), (10, 0)))
                        ]])

gaussian_col = Sg.Column([
    [Sg.Frame('Gaussian Keywords',
              layout=[
                        [Sg.Text('G16 keywords:', size=(14, 1)),
                         Sg.Input(key='-G16_KEYWORDS-', size=(44, 1),
                                  tooltip='Keywords for Gaussian16.', pad=pad1, enable_events=True)
                         ],
                        [Sg.Text('G16 nproc:', size=(14, 1)),
                         Sg.Input(key='-G16_NPROC-', size=(5, 1),
                                  tooltip='Number of processors for Gaussian16.', pad=pad1, enable_events=True),
                         Sg.Text('G16 mem (Mb):', size=(14, 1)),
                         Sg.Input(key='-G16_MEM-', size=(5, 1),
                                  tooltip='Amount of memory in Mb for Gaussian16.', pad=pad1, enable_events=True),
                         Sg.Checkbox('Write Gaussian?', key='-WRITE_GAUSSIAN-',
                                     pad=pad1, enable_events=True, default=False)
                         ],
                        [Sg.Text('Charge:', size=(14, 1)),
                         Sg.Input(key='-CHARGE-', size=(5, 1), default_text="",
                                  tooltip='Total molecular charge.', pad=pad1, enable_events=True),
                         Sg.Text('Multiplicity:', size=(14, 1)),
                         Sg.Input(key='-MULTIPLICITY-', size=(5, 1), default_text="",
                                  tooltip='Molecular multiplicity for Gaussian16.', pad=pad1, enable_events=True),
                         Sg.Checkbox('Run Gaussian?', key='-RUN_GAUSSIAN-',
                                     pad=pad1, enable_events=True, default=False)
                         ],
                        [Sg.Text("Additional info following the molecule specification: ",
                                 key="-GAUSSIAN16_EXTRAINFO_URL-", text_color='blue', enable_events=True)],
                        [Sg.Multiline(size=(50, 4), key="-GAUSSIAN16_EXTRAINFO-", enable_events=True)],
                        [Sg.Text("Path to the Gaussian16 program.: ",
                                 text_color='blue', enable_events=True)],
                        [Sg.Input(key='-GAUSSIAN16PACK-', size=(60, 1),
                                  tooltip='Path to the Gaussian16 program.', pad=pad1, enable_events=True),],
                        [Sg.Text("Gaussian program: ", size=(16, 1),
                                 pad=pad1, text_color='red'),
                         Sg.Text('https://gaussian.com/keywords/', size=(40, 1),
                                 enable_events=True, key="-GAUSSIAN16_URL-", pad=pad1, text_color='blue')
                         ],
                ], title_color='blue', pad=pad1)
     ]
])

conformers_col = Sg.Column([

    [Sg.Frame(title="Modules and environment variables to send QM package. Bash script to send to SLURM",
              layout=[
                 [Sg.Multiline(size=(76, 4), key="-BASH_EXTRAINFO-", enable_events=True,
                               tooltip="Modules and enviorement variables to be loaded by the SLURM system\n"
                                       "to run the QM package on the server")],
              ],
              title_color="blue", pad=((10, 0), (0, 0))
              ),
     Sg.Combo(['None', 'Env Var Gaussian', 'Drago', 'Cesga'],
              disabled=False, key='-ENV_COMBO-',
              default_value='None',
              tooltip="Some default parameters",
              size=(15, 1), pad=pad1, enable_events=True),

     ],
    [Sg.Frame('Conformers Keywords',
              layout=[
                        [Sg.Text("Conformer program: ", size=(22, 1), pad=pad1),
                         Sg.Combo(['rdkit', 'openbabel', 'extract from md frame', 'systematic grid',
                                   'folder with pdb files'],
                                  disabled=False, key='-CONFPACK-',
                                  default_value='rdkit',
                                  size=(36, 1), pad=pad1, enable_events=True),
                         Sg.Text("URL", enable_events=True, key="-CONF_URL-", text_color='blue'),
                         Sg.Button('Advanced conformer options', key='-BUTTONADVANCE-', disabled=False, size=(25, 1))],
                        [Sg.Text("Clustering QM method: ", size=(22, 1), pad=pad1),
                         Sg.Combo(['None', 'RMSD'],
                                  disabled=False, key='-COMBO_CLUSTER-',
                                  default_value='None',
                                  size=(36, 1), pad=pad1, enable_events=True),
                         Sg.Button('Advanced clustering options', key='-BUTTONCLUSTERADVANCE-',
                                            pad=(46, 5), disabled=False, size=(25, 1))],
                        #cJ  [Sg.Text('Number of Conformers:', size=(22, 1), enable_events=True),
                        #cJ   Sg.Input(key='-NCONF-', size=(8, 1),
                        #cJ           tooltip='Initial number of conformers.', pad=pad1, enable_events=True
                        #cJ           , default_text=0),
                        #cJ  Sg.Text('Cutoff RMSD Cluster QM (A):', size=(29, 1), enable_events=True),
                        #cJ  Sg.Input(key='-CUTOFF_RMSD_QM-', size=(12, 1),
                        #cJ           tooltip='RMSD Cutoff to clusterize the conformers after QM calculations.', pad=pad1,
                        #cJ           enable_events=True, default_text=0.0)
                        #cJ  ],
                        #cJ [Sg.Text('Minimize Iterations MM:', size=(22, 1), enable_events=True),
                        #cJ Sg.Input(key='-MIN_ITER_MM-', size=(8, 1),
                        #cJ           tooltip='Number of iterations to minimize the conformers.', pad=pad1,
                        #cJ           enable_events=True, default_text=0),
                        #cJ  Sg.Text('Cutoff Energy Threshold(kcal/mol):', size=(30, 1), enable_events=True),
                        #cJ  Sg.Input(key='-CUTOFF_ENERGY_QM-', size=(12, 1),
                        #cJ           tooltip='Energy Cutoff to clusterize the conformers after QM calculations.', pad=pad1,
                        #cJ           enable_events=True, default_text=99999.0),
                          [Sg.Checkbox('Bond Perception', key='-BOND_PERCEPTION-', size=(15, 2),
                                      pad=pad1, enable_events=True, default=False)]
                        #cJ  ],
                ], title_color='blue', pad=((10, 0), (0, 0)))
     ],
    [Sg.Frame('Path to the DockRMSD program.',
              layout=[
                  [Sg.Input(key='-DOCKRMSDPACK-', size=(97, 1), justification='right',
                            tooltip='Path to the DockRMSD program.', pad=pad1),
                   ],
                  [Sg.Text("DockRMSD program: ", size=(20, 1),
                           pad=pad1, text_color='red'),
                   Sg.Text('https://zhanglab.ccmb.med.umich.edu/DockRMSD/', size=(56, 1),
                           enable_events=True, key="-DOCK_URL-", pad=pad1, text_color='blue')
                   ],
              ], title_color='blue', pad=((10, 0), (0, 0)))
     ]
])

hide_input = Sg.Input(key='-HIDEINPUTSCRIPT-', visible=False)

general_inputs_layout = [[molecule_input], [server_input], [directories], [gaussian_col, conformers_col],
                         [hide_input], [row_buttons]]
