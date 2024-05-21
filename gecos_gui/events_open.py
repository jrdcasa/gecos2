import glob

import PySimpleGUI as Sg
import os
import sqlite3
import pandas as pd
import matplotlib.pyplot as plt
import webbrowser
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


# =============================================================================
def popup_window(loc, window, title='Window', msg=''):
    # noinspection PyBroadException
    try:
        x, y = window.size
    except Exception:
        x, y = (400, 400)
    newloc = (loc[0] + (x / 2.) * 0.5, loc[1] + (y / 2.) * 0.5)
    Sg.popup(msg,
             title=title,
             grab_anywhere=True, location=newloc)


# =============================================================================
def open_advance_window_rdkit(loc, rdkit_dict_options):
    col1 = [
        [
            Sg.Text('Number of Conformers:', size=(22, 1), enable_events=True),
            Sg.Input(key='-RDKIT_NCONF-', size=(10, 1),
                     tooltip='Number of conformers.', enable_events=True
                     , default_text=rdkit_dict_options['-RDKIT_NCONF-'])
         ],
        [
            Sg.Text('Minimize Iterations MM:', size=(22, 1), enable_events=True),
            Sg.Input(key='-RDKIT_MIN_ITER_MM-', size=(8, 1),
                     tooltip='Number of iterations to minimize the conformers.',
                     enable_events=True, default_text=rdkit_dict_options['-RDKIT_MIN_ITER_MM-']),
         ],
        [
            Sg.Text('maxattempts:', size=(22, 1)),
            Sg.Input(key='-RDKIT_MAXATTEMPTS-', enable_events=True,
                     tooltip="RDKIT: The maximum number of attempts to try embedding",
                     default_text=rdkit_dict_options['-RDKIT_MAXATTEMPTS-'], size=(10, 1))
         ],
        [
            Sg.Text('prunermsthresh:', size=(22, 1)),
            Sg.Input(key='-RDKIT_PRUNERMSTHRESH-', enable_events=True,
                     tooltip="RDKIT: Retain only the conformations out of ‘numConfs’\n"
                             "after embedding that are at least\n"
                             "this far apart from each other. RMSD is computed on the heavy atoms.\n"
                             "Pruning is greedy; i.e. the first embedded conformation is retained\n"
                             "and from then on only those that are at least\n"
                             "pruneRmsThresh away from all retained conformations are kept.\n"
                             "The pruning is done after embedding\n"
                             "and bounds violation minimization. No pruning by default\n",
                     default_text=rdkit_dict_options['-RDKIT_PRUNERMSTHRESH-'], size=(10, 1))
        ],

        [
            Sg.Checkbox('useexptorsionangleprefs', enable_events=True,
                        default=rdkit_dict_options['-RDKIT_USEEXPTORSIONANGLEPREFS-'],
                        tooltip="Impose experimental torsion angle preferences",
                        key='-RDKIT_USEEXPTORSIONANGLEPREFS-', size=(26, 1)),

        ],
        [
            Sg.Checkbox('usebasicknowledge', enable_events=True,
                        default=rdkit_dict_options['-RDKIT_USEBASICKNOWLEDGE-'],
                        tooltip="Impose basic knowledge such as flat rings",
                        key='-RDKIT_USEBASICKNOWLEDGE-', size=(26, 1)),
        ],
        [
            Sg.Checkbox('enforcechirality', enable_events=True,
                        default=rdkit_dict_options['-RDKIT_ENFORCECHIRALITY-'],
                        tooltip="Enforce the correct chirality if chiral centers are present",
                        key='-RDKIT_ENFORCECHIRALITY-', size=(26, 1)),
        ],
    ]

    # col2 = [
    #     [
    #         Sg.Checkbox('useexptorsionangleprefs', enable_events=True,
    #                     default=rdkit_dict_options['-RDKIT_USEEXPTORSIONANGLEPREFS-'],
    #                     tooltip="Impose experimental torsion angle preferences",
    #                     key='-RDKIT_USEEXPTORSIONANGLEPREFS-', size=(26, 1)),
    #
    #     ],
    #     [
    #         Sg.Checkbox('usebasicknowledge', enable_events=True,
    #                     default=rdkit_dict_options['-RDKIT_USEBASICKNOWLEDGE-'],
    #                     tooltip="Impose basic knowledge such as flat rings",
    #                     key='-RDKIT_USEBASICKNOWLEDGE-', size=(26, 1)),
    #     ],
    #     [
    #         Sg.Checkbox('enforcechirality', enable_events=True,
    #                     default=rdkit_dict_options['-RDKIT_ENFORCECHIRALITY-'],
    #                     tooltip="Enforce the correct chirality if chiral centers are present",
    #                     key='-RDKIT_ENFORCECHIRALITY-', size=(26, 1)),
    #     ],
    # ]

    col3 = [

        [Sg.Checkbox('RMSD only heavy atoms', enable_events=True,
                     default=rdkit_dict_options['-RDKIT_RMSD_ONLY_HEAVY-'],
                     tooltip="RMSD computes using only heavy atoms (not hydrogens)",
                     key='-RDKIT_RMSD_ONLY_HEAVY-', size=(22, 1)),
         ],

        [Sg.Text('Energy_threshold (kcal/mol):', size=(22, 1)),
         Sg.Input(key='-RDKIT_ENERGY_THRES-', enable_events=True,
                  tooltip="Energy threshold to assign the conformer to a cluster.\n",
                  disabled=True, text_color="gray",
                  default_text=rdkit_dict_options['-RDKIT_ENERGY_THRES-'], size=(10, 1))
         ],
        [Sg.Text('RMSD_threshold (A):', size=(22, 1)),
         Sg.Input(key='-RDKIT_RMSD_THRES-', enable_events=True,
                  tooltip="RMSD threshold to assign the conformer to a cluster.\n"
                          "RMSD is computed on the heavy atoms.\n"
                          "This parameter only works with Cluster method = RMSD",
                  disabled=True, text_color="gray",
                  default_text=rdkit_dict_options['-RDKIT_RMSD_THRES-'], size=(10, 1))
         ],
        [Sg.Text('RotConst_threshold (cm-1):', size=(22, 1)),
         Sg.Input(key='-RDKIT_ROTCONST_THRES-', enable_events=True,
                  tooltip="Rotational constants threshold to assign the conformer to a cluster.\n",
                  disabled=True, text_color="gray",
                  default_text=rdkit_dict_options['-RDKIT_ROTCONST_THRES-'], size=(10, 1))
         ],
        [Sg.Text('Window energy (kcal/mol):', size=(22, 1)),
         Sg.Input(key='-RDKIT_WINDOW_ENERGY-', enable_events=True,
                  tooltip="Window energy. Conformations with energy greatest that this are discarded\n",
                  disabled=True, text_color="gray",
                  default_text=rdkit_dict_options['-RDKIT_WINDOW_ENERGY-'], size=(10, 1))
         ],
    ]

    row = [
        [Sg.Text('force field:', size=(15, 1)),
         Sg.Combo(['MMFF', 'UFF'], enable_events=True, disabled=False, key='-RDKIT_FFNAME-',
                  default_value=rdkit_dict_options['-RDKIT_FFNAME-'], size=(40, 1))],
        [Sg.Text('cluster method:', size=(15, 1)),
         Sg.Combo(['None', 'RMSD', 'Conformer Ensemble'], enable_events=True, disabled=False,
                  key='-RDKIT_CLUSTER_METHOD-',
                  tooltip="None method: No cluster method is performed. All generated conformers are selected.\n"
                          "RMSD method: Assign conformers to different clusters based on the RMSD (heavy atoms).\n"
                          "             of the conformer to the lowest energy conformer in its cluster.\n"
                          "             The result is a number of conformers, in which the lowest energy conformer\n"
                          "             is selected. The rest are discarded.\n",
                  default_value=rdkit_dict_options['-RDKIT_CLUSTER_METHOD-'], size=(40, 1))]
    ]

    row_buttons = Sg.Column([
        [Sg.Button('Default Values', key='-RDKIT_DEFAULTVALUES-', disabled=False, size=(40, 1), enable_events=True),
         Sg.Button('CLOSE', key='-RDKIT_CLOSE-', disabled=False, size=(40, 1), enable_events=True)]
    ], )

    layout = [[Sg.Text("RDKIT options", justification='c', size=(53, 1),
                       key='-LINK_RDKITOPTIONS-', enable_events=True, text_color='blue',
                       tooltip='Click to RDKIT Manual'),
               Sg.Text("More info", justification='c',
                       key='-LINK2_RDKITOPTIONS-', enable_events=True, text_color='blue',
                       tooltip='Click to article with information of the algorithm.', size=(20, 1))],
              [Sg.Column(col1), Sg.Column(col3)], [row], [row_buttons]]
    window2 = Sg.Window("RDKIT options docs", layout, modal=True, location=loc,
                        background_color='lightblue', size=(600, 400), finalize=True)

    window2['-LINK_RDKITOPTIONS-'].bind('<Enter>', '+MOUSE OVER+')
    window2['-LINK_RDKITOPTIONS-'].bind('<Leave>', '+MOUSE AWAY+')
    window2['-LINK2_RDKITOPTIONS-'].bind('<Enter>', '+MOUSE OVER+')
    window2['-LINK2_RDKITOPTIONS-'].bind('<Leave>', '+MOUSE AWAY+')

    while True:

        event2, values2 = window2.read()

        if event2 == "-LINK_RDKITOPTIONS-+MOUSE OVER+":
            window2['-LINK_RDKITOPTIONS-'].update(text_color="green")
        if event2 == "-LINK_RDKITOPTIONS-+MOUSE AWAY+":
            window2['-LINK_RDKITOPTIONS-'].update(text_color="blue")

        if event2 == "-LINK2_RDKITOPTIONS-+MOUSE OVER+":
            window2['-LINK2_RDKITOPTIONS-'].update(text_color="green")
        if event2 == "-CONF_URL-+MOUSE AWAY+":
            window2['-LINK2_RDKITOPTIONS-'].update(text_color="blue")

        if event2 == '-LINK_RDKITOPTIONS-''-RDKIT_RMSD_THRES-':
            url = "https://www.rdkit.org/docs/source/rdkit.Chem.rdDistGeom.html"
            webbrowser.open(url)

        if event2 == '-LINK2_RDKITOPTIONS-':
            url = "https://pubs.acs.org/doi/10.1021/acs.jcim.5b00654"
            webbrowser.open(url)

        if event2 == "Exit" or event2 == Sg.WIN_CLOSED:
            break

        if event2 == '-RDKIT_CLOSE-':
            break

        if event2 == '-RDKIT_CLUSTER_METHOD-':
            if window2['-RDKIT_CLUSTER_METHOD-'].get().upper() == 'NONE':
                window2['-RDKIT_RMSD_ONLY_HEAVY-'].update(disabled=True, text_color="gray")
                window2['-RDKIT_RMSD_THRES-'].update(disabled=True, text_color="gray")
                window2['-RDKIT_ENERGY_THRES-'].update(disabled=True, text_color="gray")
                window2['-RDKIT_ROTCONST_THRES-'].update(disabled=True, text_color="gray")
                window2['-RDKIT_WINDOW_ENERGY-'].update(disabled=True, text_color="gray")
            elif window2['-RDKIT_CLUSTER_METHOD-'].get().upper() == 'RMSD':
                window2['-RDKIT_RMSD_ONLY_HEAVY-'].update(disabled=True, text_color="gray")
                window2['-RDKIT_RMSD_THRES-'].update(disabled=False, text_color="black")
                window2['-RDKIT_ENERGY_THRES-'].update(disabled=True, text_color="gray")
                window2['-RDKIT_ROTCONST_THRES-'].update(disabled=True, text_color="gray")
                window2['-RDKIT_WINDOW_ENERGY-'].update(disabled=True, text_color="gray")
            elif window2['-RDKIT_CLUSTER_METHOD-'].get().upper() == 'CONFORMER ENSEMBLE':
                window2['-RDKIT_RMSD_ONLY_HEAVY-'].update(disabled=False, text_color="black")
                window2['-RDKIT_RMSD_THRES-'].update(disabled=False, text_color="black")
                window2['-RDKIT_ENERGY_THRES-'].update(disabled=False, text_color="black")
                window2['-RDKIT_ROTCONST_THRES-'].update(disabled=False, text_color="black")
                window2['-RDKIT_WINDOW_ENERGY-'].update(disabled=False, text_color="black")

        if event2 == '-RDKIT_NCONF-':
            rdkit_dict_options['-RDKIT_NCONF-'] = window2['-RDKIT_NCONF-'].get()
        if event2 == '-RDKIT_MIN_ITER_MM-':
            rdkit_dict_options['-RDKIT_MIN_ITER_MM-'] = window2['-RDKIT_MIN_ITER_MM-'].get()
        if event2 == '-RDKIT_MAXATTEMPTS-':
            rdkit_dict_options['-RDKIT_MAXATTEMPTS-'] = window2['-RDKIT_MAXATTEMPTS-'].get()
        if event2 == '-RDKIT_PRUNERMSTHRESH-':
            rdkit_dict_options['-RDKIT_PRUNERMSTHRESH-'] = window2['-RDKIT_PRUNERMSTHRESH-'].get()
        if event2 == '-RDKIT_USEEXPTORSIONANGLEPREFS-':
            rdkit_dict_options['-RDKIT_USEEXPTORSIONANGLEPREFS-'] = window2['-RDKIT_USEEXPTORSIONANGLEPREFS-'].get()
        if event2 == '-RDKIT_USEBASICKNOWLEDGE-':
            rdkit_dict_options['-RDKIT_USEBASICKNOWLEDGE-'] = window2['-RDKIT_USEBASICKNOWLEDGE-'].get()
        if event2 == '-RDKIT_ENFORCECHIRALITY-':
            rdkit_dict_options['-RDKIT_ENFORCECHIRALITY-'] = window2['-RDKIT_ENFORCECHIRALITY-'].get()
        if event2 == '-RDKIT_FFNAME-':
            rdkit_dict_options['-RDKIT_FFNAME-'] = values2['-RDKIT_FFNAME-']
        if event2 == '-RDKIT_CLUSTER_METHOD-':
            rdkit_dict_options['-RDKIT_CLUSTER_METHOD-'] = values2['-RDKIT_CLUSTER_METHOD-']
        if event2 == '-RDKIT_RMSD_THRES-':
            rdkit_dict_options['-RDKIT_RMSD_THRES-'] = window2['-RDKIT_RMSD_THRES-'].get()

        if event2 == '-RDKIT_DEFAULTVALUES-':
            # rdkit_dict_options = defaultdict()
            rdkit_dict_options['-RDKIT_NCONF-'] = 200
            rdkit_dict_options['-RDKIT_MIN_ITER_MM-'] = 0
            rdkit_dict_options['-RDKIT_MAXATTEMPTS-'] = 1000
            rdkit_dict_options['-RDKIT_PRUNERMSTHRESH-'] = 0.1
            rdkit_dict_options['-RDKIT_USEEXPTORSIONANGLEPREFS-'] = True
            rdkit_dict_options['-RDKIT_USEBASICKNOWLEDGE-'] = True
            rdkit_dict_options['-RDKIT_ENFORCECHIRALITY-'] = True
            rdkit_dict_options['-RDKIT_FFNAME-'] = "MMFF"
            rdkit_dict_options['-RDKIT_CLUSTER_METHOD-'] = "None"
            rdkit_dict_options['-RDKIT_RMSD_THRES-'] = 1.0  # A
            rdkit_dict_options['-RDKIT_ENERGY_THRES-'] = 1.0  # kcal/mol
            rdkit_dict_options['-RDKIT_ROTCONST_THRES-'] = 0.0005  # cm-1
            rdkit_dict_options['-RDKIT_WINDOW_ENERGY-'] = 100.0  # kcal/mol
            rdkit_dict_options['-RDKIT_RMSD_ONLY_HEAVY-'] = True
            window2['-RDKIT_MAXATTEMPTS-'].update(rdkit_dict_options['-RDKIT_MAXATTEMPTS-'])
            window2['-RDKIT_PRUNERMSTHRESH-'].update(rdkit_dict_options['-RDKIT_PRUNERMSTHRESH-'])
            window2['-RDKIT_USEEXPTORSIONANGLEPREFS-'].update(rdkit_dict_options['-RDKIT_USEEXPTORSIONANGLEPREFS-'])
            window2['-RDKIT_USEBASICKNOWLEDGE-'].update(rdkit_dict_options['-RDKIT_USEBASICKNOWLEDGE-'])
            window2['-RDKIT_ENFORCECHIRALITY-'].update(rdkit_dict_options['-RDKIT_ENFORCECHIRALITY-'])
            window2['-RDKIT_ENFORCECHIRALITY-'].update(rdkit_dict_options['-RDKIT_ENFORCECHIRALITY-'])
            window2['-RDKIT_FFNAME-'].update(rdkit_dict_options['-RDKIT_FFNAME-'])
            window2['-RDKIT_CLUSTER_METHOD-'].update(rdkit_dict_options['-RDKIT_CLUSTER_METHOD-'])
            window2['-RDKIT_RMSD_THRES-'].update(rdkit_dict_options['-RDKIT_RMSD_THRES-'])
            window2['-RDKIT_ENERGY_THRES-'].update(rdkit_dict_options['-RDKIT_ENERGY_THRES-'])
            window2['-RDKIT_ROTCONST_THRES-'].update(rdkit_dict_options['-RDKIT_ROTCONST_THRES-'])
            window2['-RDKIT_WINDOW_ENERGY-'].update(rdkit_dict_options['-RDKIT_WINDOW_ENERGY-'])
            window2['-RDKIT_RMSD_ONLY_HEAVY-'].update(rdkit_dict_options['-RDKIT_RMSD_ONLY_HEAVY-'])
            window2['-RDKIT_RMSD_ONLY_HEAVY-'].update(disabled=True, text_color="gray")
            window2['-RDKIT_RMSD_THRES-'].update(disabled=True, text_color="gray")
            window2['-RDKIT_ENERGY_THRES-'].update(disabled=True, text_color="gray")
            window2['-RDKIT_ROTCONST_THRES-'].update(disabled=True, text_color="gray")
            window2['-RDKIT_WINDOW_ENERGY-'].update(disabled=True, text_color="gray")

    window2.close()


# =============================================================================
def open_advance_window_openbabel(loc, openbabel_dict_options):
    col1 = [

        [
            Sg.Text('Number of Conformers:', size=(26, 1), enable_events=True),
            Sg.Input(key='-CONFAB_NCONF-', size=(10, 1),
                     tooltip='Number of conformers.', enable_events=True
                     , default_text=openbabel_dict_options['-CONFAB_NCONF-'])
        ],

        [
            Sg.Text('Minimize Iterations MM:', size=(26, 1), enable_events=True),
            Sg.Input(key='-CONFAB_MIN_ITER_MM-', size=(10, 1),
                     tooltip='Number of iterations to minimize the conformers.',
                     enable_events=True, default_text=openbabel_dict_options['-CONFAB_MIN_ITER_MM-']),
        ],

        [
            Sg.Text('rmsd for confab diversity:', size=(26, 1)),
            Sg.Input(key='-CONFAB_RMSD_CUTOFF-', enable_events=True,
                     default_text=openbabel_dict_options['-CONFAB_RMSD_CUTOFF-'], size=(10, 1))
        ],
        [
            Sg.Text('energy cutoff for confab:', size=(26, 1)),
            Sg.Input(key='-CONFAB_ENERGY_CUTOFF-', enable_events=True,
                     default_text=openbabel_dict_options['-CONFAB_ENERGY_CUTOFF-'], size=(10, 1))
        ],
        [
            Sg.Text('rmsd for confab rmsddock:', size=(26, 1)),
            Sg.Input(key='-CONFAB_RMSD_CUTOFF_RMSDDOCK-', enable_events=True,
                     default_text=openbabel_dict_options['-CONFAB_RMSD_CUTOFF_RMSDDOCK-'], size=(10, 1))
        ],
        [
            Sg.Text('Maximum number of clusters:', size=(26, 1)),
            Sg.Input(key='-CONFAB_MAX_ENERGY_CLUSTERS-', enable_events=True,
                     default_text=openbabel_dict_options['-CONFAB_MAX_ENERGY_CLUSTERS-'], size=(10, 1))
        ],
        [
            Sg.Text('Energy threshold for clusterize:', size=(26, 1)),
            Sg.Input(key='-CONFAB_ENERGY_THRESHOLD-', enable_events=True,
                     default_text=openbabel_dict_options['-CONFAB_ENERGY_THRESHOLD-'], size=(10, 1))
        ],
    ]

    col2 = [
        [
            Sg.Checkbox('Confab verbose', enable_events=True,
                        default=openbabel_dict_options['-CONFAB_VERBOSE-'],
                        key='-CONFAB_VERBOSE-', size=(26, 1)),

        ],
    ]

    row = [
        [Sg.Text('force field:', size=(15, 1)),
         Sg.Combo(['MMFF', 'UFF'], enable_events=True, disabled=False, key='-CONFAB_FFNAME-',
                  default_value=openbabel_dict_options['-CONFAB_FFNAME-'], size=(40, 1))],
    ]

    row_buttons = Sg.Column([
        [Sg.Button('Default Values', key='-OPENBABEL_DEFAULTVALUES-', disabled=False, size=(40, 1), enable_events=True),
         Sg.Button('CLOSE', key='-OPENBABEL_CLOSE-', disabled=False, size=(40, 1), enable_events=True)]
    ], )

    layout = [[Sg.Text("OPENBABEL options", justification='c', size=(500, 1),
                       key='-LINK_OPENBABELOPTIONS-', enable_events=True)],
              [Sg.Column(col1), Sg.Column(col2)], [row], [row_buttons]]
    window3 = Sg.Window("OPENBABEL options docs", layout, modal=True, location=loc,
                        background_color='lightblue', size=(550, 350))

    while True:

        event3, values3 = window3.read()

        if event3 == '-LINK_OPENBABELOPTIONS-':
            url = 'https://openbabel.github.io/api/3.0/index.shtml'
            webbrowser.open(url)

        if event3 == "Exit" or event3 == Sg.WIN_CLOSED:
            break

        if event3 == '-OPENBABEL_CLOSE-':
            break

        if event3 == '-CONFAB_NCONF-':
            openbabel_dict_options['-CONFAB_NCONF-'] = window3['-CONFAB_NCONF-'].get()
        if event3 == '-CONFAB_RMSD_CUTOFF-':
            openbabel_dict_options['-CONFAB_RMSD_CUTOFF-'] = window3['-CONFAB_RMSD_CUTOFF-'].get()
        if event3 == '-CONFAB_ENERGY_CUTOFF-':
            openbabel_dict_options['-CONFAB_ENERGY_CUTOFF-'] = window3['-CONFAB_ENERGY_CUTOFF-'].get()
        if event3 == '-CONFAB_RMSD_CUTOFF_RMSDDOCK-':
            openbabel_dict_options['-CONFAB_RMSD_CUTOFF_RMSDDOCK-'] = window3['-CONFAB_RMSD_CUTOFF_RMSDDOCK-'].get()
        if event3 == '-CONFAB_ENERGY_THRESHOLD-':
            openbabel_dict_options['-CONFAB_ENERGY_THRESHOLD-'] = window3['-CONFAB_ENERGY_THRESHOLD-'].get()
        if event3 == '-CONFAB_MAX_ENERGY_CLUSTERS-':
            openbabel_dict_options['-CONFAB_MAX_ENERGY_CLUSTERS-'] = window3['-CONFAB_MAX_ENERGY_CLUSTERS-'].get()
        if event3 == '-CONFAB_VERBOSE-':
            openbabel_dict_options['-CONFAB_VERBOSE-'] = window3['-CONFAB_VERBOSE-'].get()
        if event3 == '-CONFAB_FFNAME-':
            openbabel_dict_options['-CONFAB_FFNAME-'] = values3['-CONFAB_FFNAME-']
        if event3 == '-OPENBABEL_DEFAULTVALUES-':
            # openbabel_dict_options = defaultdict()
            openbabel_dict_options['-CONFAB_RMSD_CUTOFF-'] = 0.5  # Angstroms
            openbabel_dict_options['-CONFAB_ENERGY_CUTOFF-'] = 50.0  # kcal/mol
            openbabel_dict_options['-CONFAB_VERBOSE-'] = False  # Verbose
            openbabel_dict_options['-CONFAB_RMSD_CUTOFF_RMSDDOCK-'] = 0.5  # Angstroms
            openbabel_dict_options['-CONFAB_FFNAME-'] = "MMFF"
            openbabel_dict_options['-CONFAB_ENERGY_THRESHOLD-'] = 99999.0
            openbabel_dict_options['-CONFAB_MAX_ENERGY_CLUSTERS-'] = 100
            openbabel_dict_options['-CONFAB_NCONF-'] = 4000000
            window3['-CONFAB_RMSD_CUTOFF-'].update(openbabel_dict_options['-CONFAB_RMSD_CUTOFF-'])
            window3['-CONFAB_ENERGY_CUTOFF-'].update(openbabel_dict_options['-CONFAB_ENERGY_CUTOFF-'])
            window3['-CONFAB_RMSD_CUTOFF_RMSDDOCK-'].update(openbabel_dict_options['-CONFAB_RMSD_CUTOFF_RMSDDOCK-'])
            window3['-CONFAB_ENERGY_THRESHOLD-'].update(openbabel_dict_options['-CONFAB_ENERGY_THRESHOLD-'])
            window3['-CONFAB_MAX_ENERGY_CLUSTERS-'].update(openbabel_dict_options['-CONFAB_MAX_ENERGY_CLUSTERS-'])
            window3['-CONFAB_VERBOSE-'].update(openbabel_dict_options['-CONFAB_VERBOSE-'])
            window3['-CONFAB_FFNAME-'].update(openbabel_dict_options['-CONFAB_FFNAME-'])

    window3.close()


# =============================================================================
def open_advance_window_extractmd(loc, extractmd_dict_options):
    col1 = [
        [
            Sg.Text('Method to extract molecule:', size=(25, 1)),
            Sg.Combo(['Sphere_COM', 'Sphere_MINATOM', 'Voronoi'], enable_events=True, disabled=False,
                     key='-METHOD_EXTRACT-',
                     default_value=extractmd_dict_options['-METHOD_EXTRACT-'], size=(17, 1),
                     tooltip="Method to calculate the negihbors")
        ],
        [
            Sg.Text('Sphere radius (angstroms):', size=(25, 1)),
            Sg.Input(key='-RADIUS_SPHERE-', enable_events=True,
                     default_text=extractmd_dict_options['-RADIUS_SPHERE-'], size=(17, 1),
                     tooltip='Sphere method needs the radius sphere.', )
        ],
    ]

    col2 = [
        [
            Sg.Radio('Extract monomers', "EXTRACT_RADIO", enable_events=True,
                     default=extractmd_dict_options['-EXTRACT_MONOMER-'],
                     key='-EXTRACT_MONOMER-', size=(26, 1)),

        ],
        [
            Sg.Radio('Extract pairs', "EXTRACT_RADIO", enable_events=True,
                     default=extractmd_dict_options['-EXTRACT_PAIR-'],
                     key='-EXTRACT_PAIR-', size=(26, 1)),
        ],
        [
            Sg.Checkbox('Extract only different \nmolecules in pairs', key='-EXTRACT_ONLYDIFFMOL-',
                        enable_events=True, default=extractmd_dict_options['-EXTRACT_ONLYDIFFMOL-'])
        ],
    ]

    row_buttons = Sg.Column([
        [Sg.Button('Default Values', key='-EXTRACT_DEFAULTVALUES-', disabled=False, size=(40, 1), enable_events=True),
         Sg.Button('CLOSE', key='-EXTRACT_CLOSE-', disabled=False, size=(40, 1), enable_events=True)]
    ], )

    layout = [[Sg.Text("EXTRACT MD options", justification='c', size=(550, 1),
                       key='-LINK_EXTRACTMDOPTIONS-', enable_events=True)],
              [Sg.Column(col1), Sg.Column(col2)], [row_buttons]]
    window2 = Sg.Window("RDKIT options docs", layout, modal=True, location=loc,
                        background_color='lightblue', size=(600, 200))

    while True:

        event2, values2 = window2.read()

        if event2 == "Exit" or event2 == Sg.WIN_CLOSED:
            extractmd_dict_options['-EXTRACT_MONOMER-'] = window2['-EXTRACT_MONOMER-'].get()
            extractmd_dict_options['-EXTRACT_PAIR-'] = window2['-EXTRACT_PAIR-'].get()
            extractmd_dict_options['-EXTRACT_ONLYDIFFMOL-'] = window2['-EXTRACT_ONLYDIFFMOL-'].get()
            break

        if event2 == '-EXTRACT_CLOSE-':
            extractmd_dict_options['-EXTRACT_MONOMER-'] = window2['-EXTRACT_MONOMER-'].get()
            extractmd_dict_options['-EXTRACT_PAIR-'] = window2['-EXTRACT_PAIR-'].get()
            extractmd_dict_options['-EXTRACT_ONLYDIFFMOL-'] = window2['-EXTRACT_ONLYDIFFMOL-'].get()
            break

        if event2 == '-METHOD_EXTRACT-':
            extractmd_dict_options['-METHOD_EXTRACT-'] = values2['-METHOD_EXTRACT-']
        if event2 == '-RADIUS_SPHERE-':
            extractmd_dict_options['-RADIUS_SPHERE-'] = window2['-RADIUS_SPHERE-'].get()
        # if event2 == '-EXTRACT_MOL-':
        #     print(window2['-EXTRACT_MOL-'].values())
        #     print("H")
        if event2 == '-EXTRACT_DEFAULTVALUES-':
            # extractmd_dict_options = defaultdict()
            extractmd_dict_options['-METHOD_EXTRACT-'] = "Sphere_COM"
            extractmd_dict_options['-RADIUS_SPHERE-'] = 7.0
            extractmd_dict_options['-EXTRACT_MONOMER-'] = True
            extractmd_dict_options['-EXTRACT_PAIR-'] = False
            extractmd_dict_options['-EXTRACT_ONLYDIFFMOL-'] = False
            window2['-METHOD_EXTRACT-'].update(extractmd_dict_options['-METHOD_EXTRACT-'])
            window2['-RADIUS_SPHERE-'].update(extractmd_dict_options['-RADIUS_SPHERE-'])
            window2['-EXTRACT_MONOMER-'].update(True)
            window2['-EXTRACT_PAIR-'].update(False)
            window2['-EXTRACT_ONLYDIFFMOL-'].update(False)

    window2.close()


# =============================================================================
def open_advance_window_systematicgrid(loc, systematicgrid_dict_options):

    text = [Sg.Text('Atom indices start at one.',
                    tooltip="Generate conformers in a systematic way using "
                            "the dihedral angles below defined.")]
    maxdihedral = 6
    send_to_qm_freeze = False

    # Default values if the dictionary is not empty
    default_cell = []
    if systematicgrid_dict_options['-SG_NDIHEDRALS-'] != 0:
        dftext = systematicgrid_dict_options['-SG_NDIHEDRALS-']
        for i in range(0, maxdihedral):
            default_cell.append([])
            for j in range(0, 6):
                try:
                    vv = systematicgrid_dict_options['-SG_DIH_STEPS-'][i][j]
                except IndexError:
                    vv = None
                default_cell[i].append(vv)
    else:
        dftext = None
        for i in range(0, maxdihedral):
            default_cell.append([])
            for j in range(0, 6):
                default_cell[i].append(None)

    dfoptimize = systematicgrid_dict_options['-SG_MM_OPTIMIZATION-']
    dfiter = systematicgrid_dict_options['-SG_MM_MAX_ITER-']
    dfqmfreeze = systematicgrid_dict_options['-SG_ADD_QM-']

    layout = [[
        Sg.Text('Number of dihedrals (1 to {}):'.format(maxdihedral), size=(26, 1)),
        Sg.Input(key='-SG_NDIHEDRALS-', enable_events=True, background_color='white',
                 size=(8, 1), default_text=dftext, disabled=False,
                 tooltip='Number of dihedrals to perform systematic grid. A number between 1 to {}'
                 .format(maxdihedral))
    ],
    ]

    head = [Sg.Text('Atom 1', justification='center', size=(10, 1), pad=(1, 1)),
            Sg.Text('Atom 2', justification='center', size=(10, 1), pad=(1, 1)),
            Sg.Text('Atom 3', justification='center', size=(10, 1), pad=(1, 1)),
            Sg.Text('Atom 4', justification='center', size=(10, 1), pad=(1, 1)),
            Sg.Text('Step', justification='center', size=(10, 1), pad=(1, 1),
                    tooltip="Step angle to scan the dihedral angle between 0 to 360º")
            ]
    layout.append(head)

    for i in range(maxdihedral):
        row = []
        for j in range(0, 6):
            if j < 5:
                row.append(Sg.Input(size=(10, 1), key='-SG_CELL_{}_{}-'.format(i, j), justification='center',
                                    background_color='white', pad=(1, 1), default_text=default_cell[i][j],
                                    disabled=True))
            else:
                row.append(Sg.Checkbox('Only Trans(180), Gauche+() and Gauche-()', size=(40, 1),
                                       key='-SG_CELL_{}_{}-'.format(i, j),
                                       enable_events=True, default=default_cell[i][j]))

        layout.append(row)

    checkboxes = [Sg.Checkbox('Optimize conformer using MMFF?', key='-SG_MM_OPTIMIZATION-',
                              enable_events=True, default=dfoptimize),
                  Sg.Text("Max. number of MM iterations:"),
                  Sg.Input(size=(6, 1), key='-SG_MM_MAX_ITER-', default_text=dfiter, enable_events=True,
                           tooltip="Maximum number of iterations for the MM optimization.", disabled=False)]

    layout.append(checkboxes)

    checkboxes2 = [Sg.Checkbox('Add dihdedral freeze to QM program?', key='-SG_ADD_QM-',
                               enable_events=True, default=dfqmfreeze,
                               tooltip="Add the dihedral angles as freeze coordiantes in the QM package.")]

    layout.append(checkboxes2)

    row_buttons = \
        [Sg.Button('Clean Data', key='-SG_CLEAN_DATA-', disabled=False, size=(20, 1), enable_events=True),
         Sg.Button('CLOSE', key='-SG_CLOSE-', disabled=False, size=(20, 1), enable_events=True)]

    layout.append(row_buttons)

    window2 = Sg.Window("Systematic Grid options docs", [text, layout], modal=True, location=loc,
                        background_color='lightblue', size=(800, 180 + maxdihedral * 30 + 40))

    first = True

    window2.read(timeout=1)
    # Activate the cells
    try:
        ndihedrals = int(window2['-SG_NDIHEDRALS-'].get())
    except ValueError:
        ndihedrals = 0
    for idx in range(ndihedrals):
        window2['-SG_NDIHEDRALS-'].update(disabled=False)
        for j in range(5):
            window2['-SG_CELL_{}_{}-'.format(idx, j)].update(disabled=False)

    while True:

        event2, values2 = window2.read()

        window2['-SG_NDIHEDRALS-'].bind("<Return>", "_Enter")
        window2['-SG_NDIHEDRALS-'].bind("<Tab>", "_Enter")

        # Event close
        if event2 == "Exit" or event2 == Sg.WIN_CLOSED or event2 == '-SG_CLOSE-':
            try:
                ndihedrals = int(window2['-SG_NDIHEDRALS-'].get())
                systematicgrid_dict_options['-SG_NDIHEDRALS-'] = ndihedrals
                success = True
            except ValueError:
                ndihedrals = 0
                success = False

            for idihedral in range(ndihedrals):
                try:
                    systematicgrid_dict_options['-SG_DIH_STEPS-'][idihedral]
                except IndexError:
                    systematicgrid_dict_options['-SG_DIH_STEPS-'].append([])
                for icol in range(0, 6):
                    try:
                        val = int(window2['-SG_CELL_{}_{}-'.format(idihedral, icol)].get())
                        systematicgrid_dict_options['-SG_DIH_STEPS-'][idihedral][icol] = val
                    except IndexError:
                        systematicgrid_dict_options['-SG_DIH_STEPS-'][idihedral].append(val)
                    except ValueError:
                        if success and ndihedrals < 7:
                            success = False
                            popup_window(loc, window2,
                                         title='Check number and type of dihedrals',
                                         msg="ERROR!!!.Number of dihedrals and dihedral "
                                             "definitions are not consistent.\n"
                                             "Value in cell {},{} must be an integer".format(idihedral, icol))

            if success:
                # Calculate the number of conformers to be generated
                nconfs = 0
                for idx in range(ndihedrals):
                    idihedral = systematicgrid_dict_options['-SG_DIH_STEPS-'][idx]
                    if idihedral[5] == 1:
                        nsteps = 3
                    else:
                        nsteps = 360 / idihedral[4]
                    if nconfs == 0:
                        nconfs = nsteps
                    else:
                        nconfs *= nsteps
                send_to_qm_freeze = window2['-SG_ADD_QM-'].get()
                systematicgrid_dict_options['-SG_ADD_QM-'] = send_to_qm_freeze
                systematicgrid_dict_options['-SG_MM_OPTIMIZATION-'] = window2['-SG_MM_OPTIMIZATION-'].get()
                popup_window(loc, window2,
                             title='Number of conformers',
                             msg="The number of conformers to be generated are {}".format(int(nconfs)))

                break
            else:
                break

        elif event2 == '-SG_CLEAN_DATA-':
            window2['-SG_NDIHEDRALS-'].update(value="", disabled=False)
            for i in range(0, maxdihedral):
                for j in range(0, 5):
                    window2['-SG_CELL_{}_{}-'.format(i, j)].update(value="", disabled=True)
                    default_cell[i][j] = None

        # Event change ndihedrals
        elif event2 == "-SG_NDIHEDRALS-_Enter":
            try:
                ndihedrals = int(window2['-SG_NDIHEDRALS-'].get())
            except ValueError:
                ndihedrals = 0
                pass
            if ndihedrals < 1 or ndihedrals > maxdihedral:
                x, y = window2.size
                newloc = (loc[0] + (x / 2.) * 0.5, loc[1] + (y / 2.) * 0.5)
                Sg.popup('Allowed values are in the range of 1 to {}'.format(maxdihedral),
                         title='Check number of dihedrals',
                         grab_anywhere=True, location=newloc)
                window2['-SG_NDIHEDRALS-'].update(value="")
                continue
            # First disable all and after enabled only the ndihedral rows
            for i in range(0, maxdihedral):
                if i < ndihedrals:
                    for j in range(0, 5):
                        window2['-SG_CELL_{}_{}-'.format(i, j)].update(disabled=False)
                else:
                    for j in range(0, 5):
                        window2['-SG_CELL_{}_{}-'.format(i, j)].update(disabled=True)
        elif event2 == "-SG_ADD_QM-":
            send_to_qm_freeze = window2['-SG_ADD_QM-'].get()
            systematicgrid_dict_options['-SG_ADD_QM-'] = send_to_qm_freeze

    window2.close()
    return send_to_qm_freeze


# =============================================================================
def open_advance_window_folderwithpdb(loc, pdbfolder_dict_options):

    if pdbfolder_dict_options['-FOLDERPDB_INPUT-'] is None:
        default_value = ""
        npdb_files = 0
    else:
        default_value = pdbfolder_dict_options['-FOLDERPDB_INPUT-']
        list_pdb_files = sorted(glob.glob(os.path.join(pdbfolder_dict_options['-FOLDERPDB_INPUT-'], "*.pdb")))
        npdb_files = len(list_pdb_files)

    pdb_folder = Sg.Frame(title="PDB Folder",
                          layout=[
                                  [Sg.Text('Folder with PDB files:', size=(18, 1)),
                                   Sg.Input(key='-FOLDERPDB_INPUT-', size=(90, 1),
                                            tooltip='Folder containing PDB files', enable_events=True,
                                            default_text=default_value),
                                   Sg.FolderBrowse(button_text="Browse", key="-FOLDERPDB_INPUT-")],
                                  [Sg.Text('Number of PDBs in folder = {}'.format(npdb_files), size=(90, 1),
                                           key='-FOLDERPDB_NPDBS-')]
                              ], title_color='blue', pad=((10, 0), (10, 0)))

    row_buttons = Sg.Column([
        [Sg.Button('CLOSE', key='-FOLDERPDB_CLOSE-', disabled=False, size=(40, 1), enable_events=True)]
    ], )

    layout = [[pdb_folder], [row_buttons]]
    window2 = Sg.Window("Folder with PDB Files", layout, modal=True, location=loc,
                        background_color='lightblue', size=(1000, 150), finalize=True)

    while True:

        event2, values2 = window2.read()

        window2['-FOLDERPDB_INPUT-'].bind("<Return>", "_Enter")
        window2['-FOLDERPDB_INPUT-'].bind("<Tab>", "_Enter")

        if event2 == "Exit" or event2 == Sg.WIN_CLOSED:
            break

        if event2 == '-FOLDERPDB_CLOSE-':
            break

        if event2 == '-FOLDERPDB_INPUT-':
            pdbfolder_dict_options['-FOLDERPDB_INPUT-'] = window2['-FOLDERPDB_INPUT-'].get()
            list_pdb_files = sorted(glob.glob(os.path.join(pdbfolder_dict_options['-FOLDERPDB_INPUT-'], "*.pdb")))
            npdb_files = len(list_pdb_files)
            window2['-FOLDERPDB_NPDBS-'].update('Number of PDBs in folder = {}'.format(npdb_files))
    window2.close()


# =============================================================================
def open_advance_window_cluster_qmrmsd(loc, open_advance_window_cluster_qmrmsd):

    text = [Sg.Text('RMSD QM clustering method.', pad=(100, 0),
                    tooltip="Cluster the conformers in accordance with the heavy atoms RMSD to the minimum.")]

    layout = \
        [Sg.Text('Cutoff RMSD Cluster QM (A):', size=(30, 1), enable_events=True),
         Sg.Input(key='-CUTOFF_RMSD_QM-', size=(12, 1),
                  tooltip='RMSD Cutoff to clusterize the conformers after QM calculations.',
                  enable_events=True, default_text=open_advance_window_cluster_qmrmsd['-CUTOFF_RMSD_QM-'])]
    layout1 = \
        [Sg.Text('Cutoff Energy Threshold(kcal/mol):', size=(30, 1), enable_events=True),
         Sg.Input(key='-CUTOFF_ENERGY_QM-', size=(12, 1),
                  tooltip='Energy Cutoff to clusterize the conformers after QM calculations.',
                  enable_events=True, default_text=open_advance_window_cluster_qmrmsd['-CUTOFF_ENERGY_QM-'])]

    layout.append(layout1)

    row_buttons = \
        [Sg.Button('Default Data', key='-RMSD_DEFAULT_DATA-', disabled=False, size=(20, 1), enable_events=True),
         Sg.Button('Close', key='-RMSD_CLOSE-', disabled=False, size=(20, 1), enable_events=True)]

    window2 = Sg.Window("Clustering method RMSD", [text, layout, row_buttons], modal=True, location=loc,
                        background_color='lightblue', size=(400, 150))

    first = True
    while True:

        event2, values2 = window2.read()
        if event2 == "Exit" or event2 == Sg.WIN_CLOSED or event2 == '-RMSD_CLOSE-':
            break

        if event2 == '-RMSD_DEFAULT_DATA-':
            window2['-CUTOFF_RMSD_QM-'].update(value=open_advance_window_cluster_qmrmsd['-CUTOFF_RMSD_QM-'])
            window2['-CUTOFF_ENERGY_QM-'].update(value=open_advance_window_cluster_qmrmsd['-CUTOFF_ENERGY_QM-'])

        if event2 == '-CUTOFF_RMSD_QM-':
            open_advance_window_cluster_qmrmsd['-CUTOFF_RMSD_QM-'] = window2['-CUTOFF_RMSD_QM-'].get()

        if event2 == '-CUTOFF_ENERGY_QM-':
            open_advance_window_cluster_qmrmsd['-CUTOFF_ENERGY_QM-'] = window2['-CUTOFF_ENERGY_QM-'].get()

    window2.close()


# =============================================================================
def open_window_log(window, loc):
    layout = [[Sg.Text("View Log", justification='c', size=(500, 1))],
              [Sg.Multiline(size=(360, 60), key='-LOG_WINDOW-', horizontal_scroll=True)]
              ]
    newloc = (loc[0] + 100, loc[1] + 30)
    window3 = Sg.Window("View Log", layout, modal=True, location=newloc,
                        background_color='lightblue', size=(1250, 700), finalize=True)

    filelog = os.path.join(window['-LOCAL_DIR-'].get(), window['-FILENAME_LOG-'].get())
    with open(filelog, 'r') as fin:
        lines = ""
        for iline in fin.readlines():
            lines += iline
        window3['-LOG_WINDOW-'].update(value=lines)

    while True:
        event3, values3 = window3.read()
        if event3 == "Exit" or event3 == Sg.WIN_CLOSED:
            break


# =============================================================================
def open_window_results(window, vmd_path):
    def draw_figure(canvas, figure):
        figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
        figure_canvas_agg.draw()
        figure_canvas_agg.get_tk_widget().pack(side='top', fill='both', expand=True)
        return figure_canvas_agg

    fulldbpath = os.path.join(window['-LOCAL_DIR-'].get(), window['-DATABASE_NAME-'].get())
    fulllogpath = os.path.join(window['-LOCAL_DIR-'].get(), window['-FILENAME_LOG-'].get())

    # Get data from the database
    db = sqlite3.connect(fulldbpath)
    df = pd.read_sql_query("SELECT * from 'qm_jobs'", db)
    headings = df.columns.values.tolist()
    data = df.values.tolist()
    newdata = []
    for item in data:
        d = list(map(str, item))
        newdata.append(d)
    del data

    # Find names of the lowest energy conformers in each cluster
    icluster = 1
    listcluster = []
    delta_elist = []
    clusterlist = []
    with open(fulllogpath) as flog:
        stringfile = flog.read()
        try:
            tmp = stringfile.split('## Optimized conformers:\n')[1]
            result = tmp.split('## Optimized conformers')[0].split("\n")[:-1]
        except IndexError:
            result = None
        if result is not None:
            for item in result:
                listcluster.append("Cluster {0:03d}: {1:s}".format(icluster, item.replace("\t", "")))
                pattern = item.replace('\t\t\t', '').replace('_allign.mol2', '')
                clusterlist.append(icluster)
                delta_elist.append(df.loc[df["name_job"] == pattern]['DeltaE'].values[0])
                icluster += 1

    # Create plot with QM energies
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    # fig, ax = plt.subplots()
    ax.bar(clusterlist, delta_elist)
    ax.set_ylabel(r'$\Delta$E (kcal/mol)')
    ax.set_xlabel('# Cluster')
    ax.set_title('Relative energy (kcal/mol)')

    pad1 = ((5, 5), (5, 5))
    pad2 = (70, 0)
    pad3 = ((50, 0), (10, 0))
    pad4 = ((10, 0), (10, 0))
    layout = [[Sg.Text('Database results (file: {})'.format(window['-DATABASE_NAME-'].get()), size=(1000, 1),
                       background_color='white', pad=pad1, justification='center', text_color='blue', )],
              [Sg.Table(newdata, headings=headings, justification='center',
                        key='-DBTABLE-', vertical_scroll_only=False, num_rows=12, pad=pad2)],
              [Sg.Canvas(key='-CANVAS-', pad=pad3, size=(480, 450)),
               Sg.Column([[Sg.Text("Lowest energy conformer of each cluster.", text_color='blue',
                                   background_color='white', justification='center', size=(51, 1))],
                          [Sg.Listbox(listcluster, size=(50, 10), background_color='white')],
                          [Sg.Text("VMD Path", text_color='blue',
                                   background_color='white', justification='center', size=(51, 1))],
                          [Sg.Input(key='-INPUTVMDPATH-', background_color='white',
                                    size=(51, 1), default_text=vmd_path)],
                          [Sg.Button("BROWSE to VMD exe", key='-BUTTONBROWSEVMD-', enable_events=True),
                           Sg.Button("Run VMD", key='-BUTTONRUNVMD-', enable_events=True)]
                          ],
                         background_color='white', size=(480, 450), pad=pad4),
               ],
              ]

    loc = window.current_location()
    newloc = (loc[0] + 300, loc[1] - 50)
    window4 = Sg.Window("Results Windows", layout, finalize=True, location=newloc, size=(1200, 800),
                        background_color='white')
    draw_figure(window4['-CANVAS-'].TKCanvas, fig)

    while True:

        event4, values4 = window4.read()
        if event4 == Sg.WINDOW_CLOSED:
            break

        if event4 == "-BUTTONBROWSEVMD-":
            loc = window.current_location()
            x, y = window.size
            newloc2 = (loc[0] + (x / 2.) * 0.5, loc[1] + (y / 2.) * 0.5)
            filename = Sg.popup_get_file('Find path to vmd binary file', no_window=False,
                                         location=newloc2, file_types=(("ALL Files", "*"),))
            if filename is not None and len(filename) > 2:
                vmd_path = filename
                window4['-INPUTVMDPATH-'].update(filename)

        if event4 == "-BUTTONRUNVMD-":
            tclpath = os.path.join(window['-LOCAL_DIR-'].get(), window['-PATTERN-'].get() + "_g16_results",
                                   "QM_optimized_conformers.tcl")
            cmd = "{} -e {}".format(vmd_path, tclpath)
            os.system(cmd)

    window4.close()
