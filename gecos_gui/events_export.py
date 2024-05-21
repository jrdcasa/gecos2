import os
from gecos_gui.events_check import popup_error


# =============================================================================
def export_jsonfile_from_gui(window, filename, rdkit_dict_options,
                             openbabel_dict_options, pass_encrypted_file,
                             save=True):
    """
    Create a json file from the GUI containing the keywords for GeCos

    Args:
        window: window instance
        filename: json file
        rdkit_dict_options:
        openbabel_dict_options:
        pass_encrypted_file:
        save: Save file

    """

    # If filename None, the CANCEL button has been pushed
    if filename is None or len(filename) == 0:
        return False

    keywords_to_write_str = {'moleculefile': '-MOLECULE_INPUT-',
                             'nameserver': '-NAME_SERVER-',
                             'username': '-USER_NAME-',
                             'keyfile': '-KEY_SSH_FILE-',
                             'encrypted_passwd': None,
                             'partition': '-SLURM_PART-',
                             'nodemaster': '-NODE_MASTER-',
                             'partitionmaster': '-SLURM_PART_MASTER-',
                             'localdir': '-LOCAL_DIR-',
                             'remotedir': '-REMOTE_DIR-',
                             'pattern': '-PATTERN-',
                             'databasename': '-DATABASE_NAME-',
                             'filelog': '-FILENAME_LOG-',
                             'g16_key': '-G16_KEYWORDS-',
                             'exec_rmsddock': '-DOCKRMSDPACK-',
                             'exec_g16': '-GAUSSIAN16PACK-',
                             'conformer_program': '-CONFPACK-'
                             }
    keywords_to_write_nostr = {'g16_nproc': '-G16_NPROC-',
                               'g16_mem': '-G16_MEM-',
                               'nconfs': '-NCONF-',
                               'cutoff_rmsd_QM': '-CUTOFF_RMSD_QM-',
                               'minimize_iterations': '-MIN_ITER_MM-',
                               'charge': '-CHARGE-',
                               'multiplicity': '-MULTIPLICITY-'
                               }
    keywords_to_write_list = {'exclude_nodes': '-EXCLUDE_NODES-'}
    keywords_to_write_combo = {'bond_perception': '-BOND_PERCEPTION-',
                               'write_gaussian': '-WRITE_GAUSSIAN-'}

    lines = "{\n"
    for key, value in keywords_to_write_str.items():
        if key.upper() == "ENCRYPTED_PASSWD":
            if len(pass_encrypted_file) > 1:
                lines += "\t\"{}\": \"{}\",\n".format(key, pass_encrypted_file)
        elif key.upper() == "NODEMASTER":
            nmaster = window[value].get()
            if len(nmaster) > 1 and nmaster.upper() != "NONE":
                lines += "\t\"{}\": \"{}\",\n".format(key, nmaster)
            else:
                lines += "\t\"{}\": \"\",\n".format(key)
        else:
            data = window[value].get()
            lines += "\t\"{}\": \"{}\",\n".format(key, data)

    for key, value in keywords_to_write_nostr.items():
        data = window[value].get()
        # For empty string assign zero value
        if not data:
            data = "0"
        lines += "\t\"{}\": {},\n".format(key, data)

    for key, value in keywords_to_write_list.items():
        data = window[value].get()
        lines += "\t\"{}\": [\n".format(key, data)
        for item in data.split(","):
            lines += "\t\t\"{}\",\n".format(item)
        lines = lines[:-2]  # Delete last comma
        lines += "\t],\n"

    for key, value in keywords_to_write_combo.items():
        data = window[value].get()
        if data:
            lines += "\t\"{}\": true,\n".format(key, data)
        else:
            lines += "\t\"{}\": false,\n".format(key, data)

    # Write advanced options only when its value is different to the default
    if window['-CONFPACK-'].get().upper() == "RDKIT":
        if rdkit_dict_options['-RDKIT_MAXATTEMPTS-'] != 1000:
            key = "rdkit_maxattempts"
            data = rdkit_dict_options['-RDKIT_MAXATTEMPTS-']
            lines += "\t\"{}\": {},\n".format(key, data)
        if rdkit_dict_options['-RDKIT_PRUNERMSTHRESH-'] != -0.01:
            key = 'rdkit_prunermsthresh'
            data = rdkit_dict_options['-RDKIT_PRUNERMSTHRESH-']
            lines += "\t\"{}\": {},\n".format(key, data)
        if not rdkit_dict_options['-RDKIT_USEEXPTORSIONANGLEPREFS-']:
            key = 'rdkit_useexptorsionangleprefs'
            lines += "\t\"{}\": {},\n".format(key, 'false')
        if not rdkit_dict_options['-RDKIT_USEBASICKNOWLEDGE-']:
            key = 'rdkit_usebasicknowledge'
            lines += "\t\"{}\": {},\n".format(key, 'false')
        if not rdkit_dict_options['-RDKIT_ENFORCECHIRALITY-']:
            key = 'rdkit_enforcechirality'
            lines += "\t\"{}\": {},\n".format(key, 'false')
        if rdkit_dict_options['-RDKIT_FFNAME-'] != "MMFF":
            key = 'rdkit_ffname'
            data = rdkit_dict_options['-RDKIT_FFNAME-']
            lines += "\t\"{}\": \"{}\",\n".format(key, data)
        if rdkit_dict_options['-RDKIT_CLUSTER_METHOD-'] != "RMSD":
            key = 'rdkit_cluster_method'
            data = rdkit_dict_options['-RDKIT_CLUSTER_METHOD-']
            lines += "\t\"{}\": \"{}\",\n".format(key, data)
        if rdkit_dict_options['-RDKIT_RMSD_THRES-'] != 2.0:
            key = 'rdkit_cluster_thres'
            data = rdkit_dict_options['-RDKIT_RMSD_THRES-']
            lines += "\t\"{}\": {},\n".format(key, data)
        lines = lines[:-2]  # Delete last comma
    elif window['-CONFPACK-'].get().upper() == "OPENBABEL":
        if openbabel_dict_options['-CONFAB_RMSD_CUTOFF-'] != 0.5:
            key = "openbabel_rmsd_cutoff_diversity"
            data = openbabel_dict_options['-CONFAB_RMSD_CUTOFF-']
            lines += "\t\"{}\": {},\n".format(key, data)
        if openbabel_dict_options['-CONFAB_ENERGY_CUTOFF-'] != 50.0:
            key = "openbabel_energy_cutoff"
            data = openbabel_dict_options['-CONFAB_ENERGY_CUTOFF-']
            lines += "\t\"{}\": {},\n".format(key, data)
        if openbabel_dict_options['-CONFAB_VERBOSE-']:
            key = 'openbabel_verbose'
            lines += "\t\"{}\": {},\n".format(key, 'false')
        if openbabel_dict_options['-CONFAB_RMSD_CUTOFF_RMSDDOCK-']:
            key = 'openbabel_rmsd_cutoff_rmsddock'
            data = openbabel_dict_options['-CONFAB_RMSD_CUTOFF_RMSDDOCK-']
            lines += "\t\"{}\": {},\n".format(key, data)
        if openbabel_dict_options['-CONFAB_FFNAME-']:
            key = 'openbabel_ffname'
            data = openbabel_dict_options['-CONFAB_FFNAME-']
            lines += "\t\"{}\": \"{}\",\n".format(key, data)
        if openbabel_dict_options['-CONFAB_ENERGY_THRESHOLD-']:
            key = 'openbabel_energy_threshold'
            data = openbabel_dict_options['-CONFAB_ENERGY_THRESHOLD-']
            lines += "\t\"{}\": {},\n".format(key, data)
        if openbabel_dict_options['-CONFAB_MAX_ENERGY_CLUSTERS-']:
            key = 'openbabel_max_energy_clusters'
            data = openbabel_dict_options['-CONFAB_MAX_ENERGY_CLUSTERS-']
            lines += "\t\"{}\": {},\n".format(key, data)
        lines = lines[:-2]  # Delete last comma
    lines += "\n}"

    # File to write
    if save and filename is not None:
        try:
            with open(filename, 'w') as f:
                f.writelines(lines)
        except FileNotFoundError:
            pass


# =============================================================================
def write_python_script_from_gui(window, filename, rdkit_dict_options, openbabel_dict_options,
                                 extractmd_dict_options, systematicgrid_dict_options, folderpdb_dict_options,
                                 cluster_qmrmsd_dict_options,
                                 pass_encrypted_file, save=True):

    # If filename None, the CANCEL button has been pushed
    if filename is None or len(filename) == 0:
        return False

    v_list_nodes = []
    for item in window['-EXCLUDE_NODES-'].get().split(","):
        v_list_nodes.append(item)

    lines = ""
    lines += "import os\n"
    lines += "import utils\n"
    lines += "import gecos\n"
    lines += "from pathlib import Path\n"
    lines += "\n"
    # KEYWORDS =====================================
    lines += "v_filename = '{}'\n".format(window['-MOLECULE_INPUT-'].get())
    lines += "v_nameserver = '{}'\n".format(window['-NAME_SERVER-'].get())
    lines += "v_username = '{}'\n".format(window['-USER_NAME-'].get())
    lines += "v_keysshfile = '{}'\n".format(window['-KEY_SSH_FILE-'].get())
    if len(pass_encrypted_file) > 1:
        lines += "v_encrypt_pass = '{}'\n".format(pass_encrypted_file)
    else:
        lines += "v_encrypt_pass = None\n"
    lines += "v_slurm_part = '{}'\n".format(window['-SLURM_PART-'].get())
    multistep = window['-BASH_EXTRAINFO-'].get()
    lines += "v_bash_extrainfo = ["
    for i in multistep.split("\n"):
        lines += "'{}',".format(i)
    lines += "]\n"
    lines += "v_list_nodes = {}\n".format(v_list_nodes)
    lines += "v_slurm_part_master = '{}'\n".format(window['-SLURM_PART_MASTER-'].get())
    if len(window['-NODE_MASTER-'].get()) > 1 and window['-NODE_MASTER-'].get().upper() != "NONE":
        lines += "v_node_master = '{}'\n".format(window['-NODE_MASTER-'].get())
    else:
        lines += "v_node_master = None\n"
    lines += "v_maxjobsslurm = {0:d}\n".format(int(window['-MAX_JOBS_SLURM-'].get()))
    lines += "v_timelimitslurm = '{}'\n".format(window['-TIME_NODES-'].get())
    lines += "v_localdir = '{}'\n".format(window['-LOCAL_DIR-'].get())
    lines += "v_remotedir = '{}'\n".format(window['-REMOTE_DIR-'].get())
    lines += "v_pattern = '{}'\n".format(window['-PATTERN-'].get())
    lines += "v_databasefullpath = '{}'\n".format(os.path.join(window['-LOCAL_DIR-'].get(),
                                                               window['-DATABASE_NAME-'].get()))
    lines += "v_fileoutputfullpath = '{}'\n".format(os.path.join(window['-LOCAL_DIR-'].get(),
                                                                 window['-FILENAME_LOG-'].get()))
    lines += "v_g16path = '{}'\n".format(window['-GAUSSIAN16PACK-'].get())
    lines += "v_g16_keywords = '{}'\n".format(window['-G16_KEYWORDS-'].get())
    multistep = window['-GAUSSIAN16_EXTRAINFO-'].get()
    lines += "v_g16_extrainfo = ["
    for i in multistep.split("\n"):
        lines += "'{}',".format(i)
    lines += "]\n"
    lines += "v_g16_ncpus = {0:d}\n".format(int(window['-G16_NPROC-'].get()))
    lines += "v_g16_mem = {0:d}\n".format(int(window['-G16_MEM-'].get()))
    lines += "v_g16_charge = {0:d}\n".format(int(window['-CHARGE-'].get()))
    lines += "v_g16_multiplicity = {0:d}\n".format(int(window['-MULTIPLICITY-'].get()))
    if window['-WRITE_GAUSSIAN-'].get():
        lines += "v_write_gaussian = True\n"
    else:
        lines += "v_write_gaussian = False\n"
    if window['-RUN_GAUSSIAN-'].get():
        lines += "v_run_gaussian = True\n"
    else:
        lines += "v_run_gaussian = False\n"

    if window['-CONFPACK-'].get().upper() == "RDKIT" or window['-CONFPACK-'].get().upper() == "OPENBABEL":
        #cJ lines += "v_nconfs = {0:d}\n".format(int(window['-NCONF-'].get()))
        #cJ lines += "v_min_iter_mm = {0:d}\n".format(int(window['-MIN_ITER_MM-'].get()))
        #cJ lines += "v_cutoff_rmsd_qm = {0:.1f}\n".format(float(window['-CUTOFF_RMSD_QM-'].get()))
        lines += "v_cutoff_rmsd_qm = {0:.1f}\n".format(float(cluster_qmrmsd_dict_options['-CUTOFF_RMSD_QM-']))
        lines += "v_energy_threshold = {0:.1f}\n".format(float(cluster_qmrmsd_dict_options['-CUTOFF_ENERGY_QM-']))
        if window['-BOND_PERCEPTION-'].get():
            lines += "v_bond_perception = True\n"
        else:
            lines += "v_bond_perception = False\n"
    else:
        #cJ lines += "v_nconfs = 0\n"
        #cJ lines += "v_min_iter_mm = 0\n"
        #cJ lines += "v_cutoff_rmsd_qm = {0:.1f}\n".format(float(window['-CUTOFF_RMSD_QM-'].get()))
        #cJ lines += "v_energy_threshold = {0:.1f}\n".format(float(window['-CUTOFF_ENERGY_QM-'].get()))
        lines += "v_cutoff_rmsd_qm = {0:.1f}\n".format(float(cluster_qmrmsd_dict_options['-CUTOFF_RMSD_QM-']))
        lines += "v_energy_threshold = {0:.1f}\n".format(float(cluster_qmrmsd_dict_options['-CUTOFF_ENERGY_QM-']))
        if window['-BOND_PERCEPTION-'].get():
            lines += "v_bond_perception = True\n"
        else:
            lines += "v_bond_perception = False\n"

    lines += "v_dockrmsdpack = '{}'\n".format(window['-DOCKRMSDPACK-'].get())
    lines += "v_confpack = '{}'\n".format(window['-CONFPACK-'].get())

    # RDKITS PARAMETERS ===================
    if window['-CONFPACK-'].get().upper() == "RDKIT":
        lines += "v_rdkit_nconfs = {0:d}\n".format(int(rdkit_dict_options['-RDKIT_NCONF-']))
        lines += "v_rdkit_min_iter_mm = {0:d}\n".format(int(rdkit_dict_options['-RDKIT_MIN_ITER_MM-']))
        lines += "v_rdkit_maxattempts = {0:d}\n".format(int(rdkit_dict_options['-RDKIT_MAXATTEMPTS-']))
        lines += "v_rdkit_prunermsthresh = {0:.3f}\n".format(float(rdkit_dict_options['-RDKIT_PRUNERMSTHRESH-']))
        if rdkit_dict_options['-RDKIT_USEEXPTORSIONANGLEPREFS-']:
            lines += "v_rdkit_useexptorsionangleprefs = True\n"
        else:
            lines += "v_rdkit_useexptorsionangleprefs = False\n"
        if rdkit_dict_options['-RDKIT_USEBASICKNOWLEDGE-']:
            lines += "v_rdkit_usebasicknowlwdge = True\n"
        else:
            lines += "v_rdkit_usebasicknowlwdge = False\n"
        if rdkit_dict_options['-RDKIT_ENFORCECHIRALITY-']:
            lines += "v_rdkit_enforcechirality = True\n"
        else:
            lines += "v_rdkit_enforcechirality = False\n"
        if rdkit_dict_options['-RDKIT_RMSD_ONLY_HEAVY-']:
            lines += "v_rdkit_rmsd_only_heavy = True\n"
        else:
            lines += "v_rdkit_rmsd_only_heavy = False\n"
        lines += "v_rdkit_cluster_method = '{}'\n".format(rdkit_dict_options['-RDKIT_CLUSTER_METHOD-'])
        lines += "v_rdkit_ffname = '{}'\n".format(rdkit_dict_options['-RDKIT_FFNAME-'])
        lines += "v_rdkit_rmsd_thres = {}\n".format(float(rdkit_dict_options['-RDKIT_RMSD_THRES-']))
        lines += "v_rdkit_energy_thres = {}\n".format(float(rdkit_dict_options['-RDKIT_ENERGY_THRES-']))
        lines += "v_rdkit_rotconst_thres = {}\n".format(float(rdkit_dict_options['-RDKIT_ROTCONST_THRES-']))
        lines += "v_rdkit_window_energy = {}\n".format(float(rdkit_dict_options['-RDKIT_WINDOW_ENERGY-']))
    # CONFAB PARAMETERS ===================
    elif window['-CONFPACK-'].get().upper() == "OPENBABEL":
        lines += "v_openbabel_nconf = {0:d}\n".format(int(openbabel_dict_options['-CONFAB_NCONF-']))
        lines += "v_openbabel_min_iter_mm = {0:d}\n".format(int(openbabel_dict_options['-CONFAB_MIN_ITER_MM-']))
        lines += "v_openbabel_rmsd_cutoff_confab = {0:f}\n".\
            format(float(openbabel_dict_options['-CONFAB_RMSD_CUTOFF-']))
        lines += "v_openbabel_energy_cutoff_confab = {0:f}\n".\
            format(float(openbabel_dict_options['-CONFAB_ENERGY_CUTOFF-']))
        if openbabel_dict_options['-CONFAB_VERBOSE-']:
            lines += "v_openbabel_verbose = True\n"
        else:
            lines += "v_openbabel_verbose = False\n"
        lines += "v_openbabel_rmsddock_confab = {0:f}\n".\
            format(float(openbabel_dict_options['-CONFAB_RMSD_CUTOFF_RMSDDOCK-']))
        lines += "v_openbabel_ffname = '{}'\n".format(openbabel_dict_options['-CONFAB_FFNAME-'])

        lines += "v_openbabel_cluster_energy_threshold = {0:f}\n".\
            format(float(openbabel_dict_options['-CONFAB_ENERGY_THRESHOLD-']))
        lines += "v_openbabel_cluster_max_number_cluster = {0:d}\n".\
            format(int(openbabel_dict_options['-CONFAB_MAX_ENERGY_CLUSTERS-']))

    elif window['-CONFPACK-'].get().upper() == "EXTRACT FROM MD FRAME":
        lines += "v_filemolextractfullpath = \"{}\"\n".format(window['-MOLECULE_INPUT-'].get())
        lines += "v_extract_method = \"{}\"\n".format(extractmd_dict_options['-METHOD_EXTRACT-'])
        lines += "v_radius_sphere = {}\n".format(float(extractmd_dict_options['-RADIUS_SPHERE-']))
        if extractmd_dict_options['-EXTRACT_MONOMER-']:
            lines += "v_extract_monomers = True\n"
            lines += "v_extract_type = \"monomer\"\n"
        else:
            lines += "v_extract_monomers = False\n"
        if extractmd_dict_options['-EXTRACT_PAIR-']:
            lines += "v_extract_pairs = True\n"
            lines += "v_extract_type = \"pair\"\n"
        else:
            lines += "v_extract_pairs = False\n"
        if extractmd_dict_options['-EXTRACT_ONLYDIFFMOL-']:
            lines += "v_extract_onlydifferentmols = True\n"
        else:
            lines += "v_extract_onlydifferentmols = False\n"
    elif window['-CONFPACK-'].get().upper() == "SYSTEMATIC GRID":
        lines += "v_sg_ndihedrals = {}\n".format(systematicgrid_dict_options['-SG_NDIHEDRALS-'])
        lines += "v_sg_listdih = {}\n".format(systematicgrid_dict_options['-SG_DIH_STEPS-'])
        lines += "v_sg_MMoptimization = {}\n".format(systematicgrid_dict_options['-SG_MM_OPTIMIZATION-'])
        lines += "v_sg_maxMMoptiter = {}\n".format(systematicgrid_dict_options['-SG_MM_MAX_ITER-'])
        lines += "v_sg_add_freeze_qm = {}\n".format(systematicgrid_dict_options['-SG_ADD_QM-'])
    elif window['-CONFPACK-'].get().upper() == 'FOLDER WITH PDB FILES':
        lines += "v_pdbfolder = '{}'\n".format(folderpdb_dict_options['-FOLDERPDB_INPUT-'])

    lines += "\n"
    lines += "cwd = Path(__file__).parent.absolute()\n"
    lines += "donepath = os.path.join(cwd, 'done')\n"

    # LOOP
    lines += "\n"
    lines += "if not os.path.isfile(v_databasefullpath) and  not os.path.isfile(donepath):\n\n"
    lines += "    log = utils.init_logger(\n" \
             "        \"Output2\",\n" \
             "        fileoutput=v_fileoutputfullpath,\n" \
             "        append=False, inscreen=True)\n\n"

    if window['-CONFPACK-'].get().upper() == "RDKIT":
        lines += "    g1 = gecos.GecosRdkit(\n" \
                 "        filename=v_filename,\n" \
                 "        total_charge=v_g16_charge,\n" \
                 "        bond_perception=v_bond_perception,\n" \
                 "        logger=log)\n"
        lines += "\n"
        lines += "    g1.generate_conformers(\n" \
                 "        v_localdir,\n" \
                 "        nconfs=v_rdkit_nconfs,\n" \
                 "        minimize_iterations=v_rdkit_min_iter_mm,\n" \
                 "        maxattempts=v_rdkit_maxattempts,\n" \
                 "        prunermsthresh=v_rdkit_prunermsthresh,\n" \
                 "        useexptorsionangleprefs=v_rdkit_useexptorsionangleprefs,\n" \
                 "        usebasicknowledge=v_rdkit_usebasicknowlwdge,\n" \
                 "        enforcechirality=v_rdkit_enforcechirality,\n" \
                 "        ff_name=v_rdkit_ffname,\n" \
                 "        cluster_method=v_rdkit_cluster_method,\n" \
                 "        rmsd_only_heavy = v_rdkit_rmsd_only_heavy,\n" \
                 "        rmsd_threshold=v_rdkit_rmsd_thres,\n" \
                 "        energy_threshold=v_rdkit_energy_thres,\n" \
                 "        rotconst_threshold=v_rdkit_rotconst_thres,\n" \
                 "        window_energy=v_rdkit_window_energy,\n" \
                 "        write_gaussian=v_write_gaussian,\n" \
                 "        pattern=v_pattern,\n" \
                 "        g16_key=v_g16_keywords,\n" \
                 "        g16_nproc=v_g16_ncpus,\n" \
                 "        g16_mem=v_g16_mem,\n" \
                 "        g16_extra_info=v_g16_extrainfo,\n" \
                 "        charge=v_g16_charge,\n" \
                 "        multiplicity=v_g16_multiplicity)\n"

    elif window['-CONFPACK-'].get().upper() == "OPENBABEL":
        lines += "    g1 = gecos.GecosPyBabel(\n" \
                 "        filename=v_filename,\n" \
                 "        exec_rmsddock=v_dockrmsdpack,\n" \
                 "        total_charge=v_g16_charge,\n" \
                 "        bond_perception=v_bond_perception,\n" \
                 "        logger=log)\n"
        lines += "\n"
        lines += "    g1.generate_conformers(\n" \
                 "            v_localdir,\n" \
                 "            nconfs=v_openbabel_nconf,\n" \
                 "            minimize_iterations=v_openbabel_min_iter_mm,\n" \
                 "            rmsd_cutoff_confab=v_openbabel_rmsd_cutoff_confab,\n" \
                 "            energy_cutoff_confab=v_openbabel_energy_cutoff_confab,\n" \
                 "            confab_verbose_confab=v_openbabel_verbose,\n" \
                 "            cutoff_rmsddock_confab=v_openbabel_rmsddock_confab,\n" \
                 "            energy_threshold_cluster=v_openbabel_cluster_energy_threshold,\n" \
                 "            max_number_cluster=v_openbabel_cluster_max_number_cluster,\n" \
                 "            ff_name=v_openbabel_ffname,\n" \
                 "            pattern=v_pattern,\n" \
                 "            write_gaussian=v_write_gaussian,\n" \
                 "            g16_key=v_g16_keywords,\n" \
                 "            g16_nproc=v_g16_ncpus,\n" \
                 "            g16_mem=v_g16_mem,\n" \
                 "            g16_extra_info=v_g16_extrainfo,\n" \
                 "            charge=v_g16_charge,\n" \
                 "            multiplicity=v_g16_multiplicity\n" \
                 "        )\n"
        lines += "\n"

    elif window['-CONFPACK-'].get().upper() == "EXTRACT FROM MD FRAME":

        lines += "    g1 = gecos.GecosExtractNeighbors(\n" \
                 "                           v_filemolextractfullpath,\n" \
                 "                           write_tcl=True,\n" \
                 "                           pattern=v_pattern,\n" \
                 "                           radius=v_radius_sphere,\n" \
                 "                           logger=log)\n"
        lines += "\n"
        lines += "    g1.extract_conformers(\n" \
                 "            v_localdir,\n" \
                 "            calc_com=True,\n" \
                 "            method=v_extract_method,\n" \
                 "            extract_type=v_extract_type,\n" \
                 "            write_gaussian=v_write_gaussian,\n" \
                 "            g16_key=v_g16_keywords,\n" \
                 "            g16_nproc=v_g16_ncpus,\n" \
                 "            g16_mem=v_g16_mem,\n" \
                 "            g16_extra_info=v_g16_extrainfo,\n" \
                 "            pattern=v_pattern,\n" \
                 "            charge=v_g16_charge,\n" \
                 "            multiplicity=v_g16_multiplicity,\n" \
                 "            onlydifferentmols=v_extract_onlydifferentmols)\n"

    elif window['-CONFPACK-'].get().upper() == "SYSTEMATIC GRID":

        lines += "    g1 = gecos.GecosSystematicGrid(\n" \
                 "                           v_filename,\n" \
                 "                           pattern=v_pattern,\n" \
                 "                           total_charge=v_g16_charge,\n" \
                 "                           bond_perception=v_bond_perception,\n" \
                 "                           logger=log)\n"
        lines += "\n"
        lines += "    g1.generate_systematic_conformers(\n" \
                 "            v_localdir,\n" \
                 "            optimize=v_sg_MMoptimization,\n" \
                 "            maxmmoptiters=v_sg_maxMMoptiter,\n" \
                 "            ndihedrals=v_sg_ndihedrals,\n" \
                 "            dih_list=v_sg_listdih,\n"  \
                 "            write_gaussian=v_write_gaussian,\n" \
                 "            g16_key=v_g16_keywords,\n" \
                 "            g16_nproc=v_g16_ncpus,\n" \
                 "            g16_mem=v_g16_mem,\n" \
                 "            g16_extra_info=v_g16_extrainfo,\n" \
                 "            pattern=v_pattern,\n" \
                 "            charge=v_g16_charge,\n" \
                 "            multiplicity=v_g16_multiplicity)\n" \

    elif window['-CONFPACK-'].get().upper() == 'FOLDER WITH PDB FILES':

        lines += "    g1 = gecos.GecosPdbFolder(\n" \
                 "                           v_pdbfolder,\n" \
                 "                           logger=log)\n"
        lines += "\n"

        lines += "    g1.generate_pdbfolder(\n" \
                  "            write_gaussian=v_write_gaussian,\n" \
                 "            g16_key=v_g16_keywords,\n" \
                 "            g16_nproc=v_g16_ncpus,\n" \
                 "            g16_mem=v_g16_mem,\n" \
                 "            g16_extra_info=v_g16_extrainfo,\n" \
                 "            pattern=v_pattern,\n" \
                 "            charge=v_g16_charge,\n" \
                 "            multiplicity=v_g16_multiplicity)\n" \


    else:
        msg = "Package {} to calculate conformers is not available.".format(window['-CONFPACK-'].get())
        popup_error(window, msg)

    lines += "\n"
    lines += "    if v_run_gaussian:\n"
    lines += "        gecos.send_qm_conformers(" \
             "\n            v_nameserver," \
             "\n            v_databasefullpath," \
             "\n            v_username," \
             "\n            v_keysshfile," \
             "\n            v_localdir," \
             "\n            v_remotedir," \
             "\n            v_g16path," \
             "\n            maxjobsslurm=v_maxjobsslurm," \
             "\n            regex='*g16*/*.com'," \
             "\n            partition=v_slurm_part," \
             "\n            exclude_nodes=v_list_nodes," \
             "\n            ncpus=v_g16_ncpus, " \
             "\n            partitionmaster=v_slurm_part_master," \
             "\n            nodemaster=v_node_master," \
             "\n            mem=v_g16_mem," \
             "\n            timelimit=v_timelimitslurm," \
             "\n            encrypted_pass=v_encrypt_pass," \
             "\n            extraslurminfo=v_bash_extrainfo," \
             "\n            logger=log)\n"
    lines += "    else:\n"
    lines += "        m = 'QM calculations will not be performed'\n"
    lines += "        print(m) if log is None else log.info(m)\n"
    lines += "        with open('done', 'w') as f:\n"
    lines += "            f.writelines('')\n"

    lines += "\n"
    lines += "else:\n"
    lines += "\n"
    lines += "    log = utils.init_logger(" \
             "\n            \"Output2\"," \
             "\n            fileoutput=v_fileoutputfullpath," \
             "\n            append=True," \
             "\n            inscreen=True)\n"
    lines += "\n"
    lines += "    v_outdir = os.path.join(v_localdir, v_pattern + '_g16_results')\n"
    lines += "\n"
    lines += "    if not os.path.isfile(donepath):\n"
    lines += "        gecos.check_qm_jobs(" \
             "\n            v_nameserver," \
             "\n            v_databasefullpath," \
             "\n            v_username," \
             "\n            v_keysshfile," \
             "\n            v_localdir," \
             "\n            v_remotedir," \
             "\n            v_outdir," \
             "\n            v_pattern," \
             "\n            v_dockrmsdpack," \
             "\n            encrypted_pass=v_encrypt_pass," \
             "\n            cutoff_rmsd=v_cutoff_rmsd_qm," \
             "\n            energy_threshold=v_energy_threshold," \
             "\n            logger=log)\n"

    lines += "\nprint(\"Job Done!!!\")\n"

    # File to write
    if save and filename is not None:
        try:
            with open(filename, 'w') as f:
                f.writelines(lines)
        except FileNotFoundError:
            pass


# =============================================================================
def write_python_script_prop_from_gui(window, filename,
                                      dict_properties, save=True):

    # If filename None, the CANCEL button has been pushed
    if filename is None or len(filename) == 0:
        return False

    lines = ""
    lines += "import os\n"
    lines += "import utils\n"
    lines += "import utils.gaussian16\n"
    lines += "import gecos\n"
    lines += "\n"

    for ikey, value in dict_properties.items():

        if isinstance(value, str):
            if ikey == 'p_node_master':
                if len(window['-NODE_MASTER-'].get()) > 1 and window['-NODE_MASTER-'].get().upper() != "NONE":
                    lines += "p_node_master = '{}'\n".format(window['-NODE_MASTER-'].get())
                else:
                    lines += "p_node_master = None\n"
            elif ikey == 'p_encrypt_pass':
                if len(value) > 1:
                    lines += "p_encrypt_pass = '{}'\n".format(value)
                else:
                    lines += "p_encrypt_pass = None\n"
            else:
                lines += "{} = '{}'\n".format(ikey, dict_properties[ikey])
        elif isinstance(value, list) and ikey == 'p_list_nodes':
            v_list_nodes = []
            for item in window['-EXCLUDE_NODES-'].get().split(","):
                v_list_nodes.append(item)
            lines += "{} = {}\n".format(ikey, v_list_nodes)
        elif isinstance(value, list) and ikey == 'p_fileproplist':
            p_fileproplist = []
            for item in window['-LISTBOX_MOL2-'].Values:
                p_fileproplist.append(item)
            lines += "{} = {}\n".format(ikey, p_fileproplist)
        elif isinstance(value, list) and ikey == 'p_bash_extrainfo':
            lines += "{} = {}\n".format(ikey, value)
        elif isinstance(value, int):
            lines += "{} = {}\n".format(ikey, dict_properties[ikey])
        elif ikey == 'p_run_gaussian':
            if value:
                lines += "{} = True\n".format(ikey)
            else:
                lines += "{} = False\n".format(ikey)
    if dict_properties['p_g16_keywords'].count("freq") == 1 and dict_properties['p_g16_keywords'].count("wfn") == 0:
        p_pattern = "FREQ"
    elif dict_properties['p_g16_keywords'].count("wfn") == 1:
        p_pattern = "WFN"
    else:
        p_pattern = "SP"
    lines += "{} = '{}'\n".format('p_pattern', p_pattern)

    regex = '*_{}_*.com'.format(p_pattern)
    lines += "p_dockrmsdpack = '{}'\n".format(window['-DOCKRMSDPACK-'].get())
    lines += "p_maxjobsslurm = {0:d}\n".format(int(window['-MAX_JOBS_SLURM-'].get()))
    lines += "p_cutoffenergy = {0:f}\n".format(float(window['-QM_PROP_CUTOFF_ENERGY-'].get()))

    # LOOP
    lines += "\n"
    lines += "if not os.path.isfile(p_databasefullpath):\n\n"

    lines += "    log = utils.init_logger(\n" \
             "        \"Output3\",\n" \
             "        fileoutput=p_fileoutputfullpath,\n" \
             "        append=False, inscreen=False)\n\n"

    lines += "    utils.print_header(log)\n\n"

    lines += "    utils.gaussian16.write_gaussian_from_mol2("\
             "\n            p_fileproplist," \
             "\n            p_localdir_mol2conf, " \
             "\n            p_localdir, " \
             "\n            pattern=p_pattern, " \
             "\n            g16_key=p_g16_keywords," \
             "\n            g16_nproc=p_ncpus," \
             "\n            g16_mem=p_mem," \
             "\n            charge=p_charge," \
             "\n            multiplicity=p_multiplicity)\n\n"

    lines += "    if p_run_gaussian:\n"
    lines += "        gecos.send_qm_conformers(" \
             "\n            p_nameserver," \
             "\n            p_databasefullpath," \
             "\n            p_username," \
             "\n            p_keysshfile," \
             "\n            p_localdir,"\
             "\n            p_remotedir," \
             "\n            p_g16path," \
             "\n            maxjobsslurm = p_maxjobsslurm,"\
             "\n            regex='{}'," \
             "\n            partition=p_slurm_part," \
             "\n            exclude_nodes=p_list_nodes," \
             "\n            ncpus=p_ncpus, " \
             "\n            partitionmaster=p_slurm_part_master," \
             "\n            nodemaster=p_node_master," \
             "\n            mem=p_mem," \
             "\n            encrypted_pass=p_encrypt_pass,"\
             "\n            extraslurminfo=p_bash_extrainfo,"\
             "\n            logger=log)\n".format(regex)

    lines += "\n"
    lines += "else:\n"
    lines += "\n"
    lines += "    log = utils.init_logger(" \
             "\n            \"Output3\"," \
             "\n            fileoutput=p_fileoutputfullpath," \
             "\n            append=True," \
             "\n            inscreen=False)\n"

    lines += "    p_outdir = os.path.join(p_localdir, p_pattern + '_g16_results')\n"
    lines += "    p_cutoff_rmsd_qm = 10.0\n"
    lines += "\n"
    lines += "    gecos.check_qm_jobs(" \
             "\n            p_nameserver," \
             "\n            p_databasefullpath," \
             "\n            p_username," \
             "\n            p_keysshfile," \
             "\n            p_localdir," \
             "\n            p_remotedir," \
             "\n            p_outdir," \
             "\n            p_pattern," \
             "\n            p_dockrmsdpack," \
             "\n            encrypted_pass=p_encrypt_pass," \
             "\n            cutoff_rmsd=p_cutoff_rmsd_qm," \
             "\n            logger=log)\n"

    lines += 'print ("Job Done!!!!")'

    # File to write
    if save and filename is not None:
        try:
            with open(filename, 'w') as f:
                f.writelines(lines)
        except FileNotFoundError:
            pass
