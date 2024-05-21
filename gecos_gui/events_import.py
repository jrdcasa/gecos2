import re
import ast


# =============================================================================
def import_pythonfile_to_gui(window, filename, rdkit_dict_options,
                             openbabel_dict_options, extractmd_dict_options,
                             systematicgrid_dict_options, pdbfolder_dict_options):

    """
    Use a json file containing the keywords for GeCos and load the parameters in the GUI

    Args:
        window: window instance
        filename: python file
        rdkit_dict_options:
        openbabel_dict_options
        extractmd_dict_options
        systematicgrid_dict_options
        pdbfolder_dict_options

    """

    # If filename None, the CANCEL button has been pushed
    if filename is None or len(filename) == 0:
        return False

    # Open python file
    var_to_dict_str = {
        'v_filename': '-MOLECULE_INPUT-',
        'v_nameserver': '-NAME_SERVER-',
        'v_username': '-USER_NAME-',
        'v_keysshfile': '-KEY_SSH_FILE-',
        'v_slurm_part': '-SLURM_PART-',
        'v_slurm_part_master': '-SLURM_PART_MASTER-',
        'v_localdir': '-LOCAL_DIR-',
        'v_remotedir': '-REMOTE_DIR-',
        'v_pattern ': '-PATTERN-',
        'v_g16_keywords': '-G16_KEYWORDS-',
        'v_dockrmsdpack': '-DOCKRMSDPACK-',
        'v_g16path': '-GAUSSIAN16PACK-',
        'v_timelimitslurm': '-TIME_NODES-'
    }

    var_to_encrypt = {'v_encrypt_pass': '-ENCRYPT_PASS_FILE-',
                      'v_node_master': '-NODE_MASTER-',
                      }

    var_to_dict_combobox = {
        'v_confpack': '-CONFPACK-'
    }

    var_to_dict_list = {
        'v_list_nodes': '-EXCLUDE_NODES-',
    }

    var_to_dict_multipleline = {
        'v_g16_extrainfo': '-GAUSSIAN16_EXTRAINFO-',
        'v_bash_extrainfo': "-BASH_EXTRAINFO-"
    }

    var_to_dict_pathlast = {
        'v_databasefullpath': '-DATABASE_NAME-',
        'v_fileoutputfullpath': '-FILENAME_LOG-'
    }

    var_to_dict_intfloat = {
        'v_g16_ncpus': '-G16_NPROC-',
        'v_g16_mem': '-G16_MEM-',
        'v_g16_charge': '-CHARGE-',
        'v_g16_multiplicity': '-MULTIPLICITY-',
        #cJ 'v_cutoff_rmsd_qm': '-CUTOFF_RMSD_QM-',
        #cJ 'v_energy_threshold': '-CUTOFF_ENERGY_QM-',
        'v_maxjobsslurm': '-MAX_JOBS_SLURM-'
    }

    var_to_dict_boolean = {
        'v_write_gaussian': '-WRITE_GAUSSIAN-',
        'v_run_gaussian': '-RUN_GAUSSIAN-',
        'v_bond_perception': '-BOND_PERCEPTION-',
    }

    var_to_dict_rdkit_str = {
        'v_rdkit_cluster_method': '-RDKIT_CLUSTER_METHOD-',
        'v_rdkit_ffname': '-RDKIT_FFNAME-'
    }
    var_to_dict_rdkit_boolean = {
        'v_rdkit_useexptorsionangleprefs': '-RDKIT_USEEXPTORSIONANGLEPREFS-',
        'v_rdkit_usebasicknowlwdge': '-RDKIT_USEBASICKNOWLEDGE-',
        'v_rdkit_enforcechirality': '-RDKIT_ENFORCECHIRALITY-',
        'v_rdkit_rmsd_only_heavy': '-RDKIT_RMSD_ONLY_HEAVY-'
    }

    var_to_dict_rdkit_intfloat = {
        'v_rdkit_maxattempts': '-RDKIT_MAXATTEMPTS-',
        'v_rdkit_prunermsthresh': '-RDKIT_PRUNERMSTHRESH-',
        'v_rdkit_rmsd_thres': '-RDKIT_RMSD_THRES-',
        'v_rdkit_nconfs': '-RDKIT_NCONF-',
        'v_rdkit_min_iter_mm': '-RDKIT_MIN_ITER_MM-',
        'v_rdkit_energy_thres': '-RDKIT_ENERGY_THRES-',
        'v_rdkit_rotconst_thres': '-RDKIT_ROTCONST_THRES-',
        'v_rdkit_window_energy': '-RDKIT_WINDOW_ENERGY-',
    }

    var_to_dict_openbabel_intfloat = {
        'v_openbabel_nconf': '-CONFAB_NCONF-',
        'v_openbabel_min_iter_mm': '-CONFAB_MIN_ITER_MM-',
        'v_openbabel_rmsd_cutoff_confab': '-CONFAB_RMSD_CUTOFF-',
        'v_openbabel_energy_cutoff_confab': '-CONFAB_ENERGY_CUTOFF-',
        'v_openbabel_rmsddock_confab': '-CONFAB_RMSD_CUTOFF_RMSDDOCK-',
        'v_openbabel_cluster_energy_threshold': '-CONFAB_ENERGY_THRESHOLD-',
        'v_openbabel_cluster_max_number_cluster': '-CONFAB_MAX_ENERGY_CLUSTERS-'
    }

    var_to_dict_openbabel_boolean = {
        'v_openbabel_verbose': '-CONFAB_VERBOSE-'
    }

    var_to_dict_openbabel_str = {
        'v_openbabel_ffname': '-CONFAB_FFNAME-'
    }

    var_to_dict_extractmd_str = {
        'v_filemolextractfullpath': '-MOLECULE_INPUT-',
        'v_extract_method': '-METHOD_EXTRACT-',
    }

    var_to_dict_extractmd_intfloat = {
        'v_radius_sphere': '-RADIUS_SPHERE-',
    }

    var_to_dict_extractmd_boolean = {
        'v_extract_monomers': '-EXTRACT_MONOMER-',
        'v_extract_pairs': '-EXTRACT_PAIR-',
        'v_extract_onlydifferentmols': '-EXTRACT_ONLYDIFFMOL-'
    }

    var_to_dict_systematicgrid_intfloat = {
        'v_sg_ndihedrals': '-SG_NDIHEDRALS-',
        'v_sg_maxMMoptiter': '-SG_MM_MAX_ITER-'
    }

    var_to_dict_systematicgrid_list = {
        'v_sg_listdih': '-SG_DIH_STEPS-',
    }

    var_to_dict_systematicgrid_boolean = {
        'v_sg_MMoptimization': '-SG_MM_OPTIMIZATION-',
        'v_sg_add_freeze_qm': '-SG_ADD_QM-',
    }

    var_to_dict_pdbfolder_str = {
        'v_pdbfolder': '-FOLDERPDB_INPUT-',
    }

    with open(filename, 'r') as fpython:
        data = fpython.read()

        for ikey in var_to_dict_str.keys():
            regex = ikey+".*="
            if re.search(regex, data) is not None:
                str2 = re.search(ikey+".*=.*", data).group(0)
                # value = re.search("(?:'|\").*(?:'|\")", str2).group(0)
                value = re.search("(?:['\"]).*(?:['\"])", str2).group(0)
                value = value.replace("'", "")
                value = value.replace("\"", "")
                label = var_to_dict_str[ikey]
                window[label].update(value)

        for ikey in var_to_encrypt.keys():
            regex = ikey + ".*="
            if re.search(regex, data) is not None:
                str2 = re.search(ikey + ".*=.*", data).group(0)
                tmp = re.search("(?:['\"]).*(?:['\"])", str2)
                if tmp is not None:
                    value = tmp.group(0)
                    value = value.replace("'", "")
                    value = value.replace("\"", "")
                    label = var_to_encrypt[ikey]
                    window[label].update(value)

        for ikey in var_to_dict_combobox.keys():
            regex = ikey+".*="
            if re.search(regex, data) is not None:
                str2 = re.search(ikey+".*=.*", data).group(0)
                value = re.search("(?:['\"]).*(?:['\"])", str2).group(0)
                value = value.replace("'", "")
                value = value.replace("\"", "")
                label = var_to_dict_combobox[ikey]
                window[label].update(value)
                # TODO: This must be revised, the names have been changed
                # if value.upper() == 'EXTRACT FROM MD FRAME':
                #     window['-NCONF-'].update(disabled=True, text_color='grey')
                #     window['-MIN_ITER_MM-'].update(disabled=True, text_color='grey')
                #     window['-CUTOFF_RMSD_QM-'].update(disabled=False)
                #     window['-CUTOFF_ENERGY_QM-'].update(disabled=False)

        for ikey in var_to_dict_list.keys():
            regex = ikey+".*="
            if re.search(regex, data) is not None:
                str2 = re.search(ikey+".*=.*", data).group(0)
                # value = re.search('(?:\'|").*(?:\'|")', str2).group(0)
                value = re.search("(?:['\"]).*(?:['\"])", str2).group(0)
                value = value.replace("'", "")
                value = value.replace("\"", "")
                label = var_to_dict_list[ikey]
                window[label].update(value)

        for ikey in var_to_dict_pathlast.keys():
            regex = ikey+".*="
            if re.search(regex, data) is not None:
                str2 = re.search(ikey+".*=.*", data).group(0)
                # value = re.search('(?:\'|").*(?:\'|")', str2).group(0)
                value = re.search("(?:['\"]).*(?:['\"])", str2).group(0)
                value = value.replace("'", "")
                value = value.replace("\"", "")
                value = value.split("/")[-1]
                label = var_to_dict_pathlast[ikey]
                window[label].update(value)

        for ikey in var_to_dict_intfloat.keys():
            regex = ikey+".*="
            if re.search(regex, data) is not None:
                value = re.search(ikey+".*=.*", data).group(0)
                value = value.split("=")[-1]
                label = var_to_dict_intfloat[ikey]
                window[label].update(value)

        for ikey in var_to_dict_boolean.keys():
            regex = ikey+".*="
            if re.search(regex, data) is not None:
                value = re.search(ikey+".*=.*", data).group(0)
                value = value.replace("'", "")
                value = value.replace("\"", "")
                value = value.split("=")[-1]
                if value.find('True') != -1:
                    valbol = True
                else:
                    valbol = False
                label = var_to_dict_boolean[ikey]
                window[label].update(value=valbol)

        for ikey in var_to_dict_rdkit_str.keys():
            regex = ikey+".*="
            if re.search(regex, data) is not None:
                str2 = re.search(ikey+".*=.*", data).group(0)
                # value = re.search('(?:\'|").*(?:\'|")', str2).group(0)
                value = re.search("(?:['\"]).*(?:['\"])", str2).group(0)
                value = value.replace("'", "")
                value = value.replace("\"", "")
                label = var_to_dict_rdkit_str[ikey]
                rdkit_dict_options[label] = value

        for ikey in var_to_dict_rdkit_intfloat.keys():
            regex = ikey+".*="
            if re.search(regex, data) is not None:
                value = re.search(ikey+".*=.*", data).group(0)
                value = value.split("=")[-1]
                label = var_to_dict_rdkit_intfloat[ikey]
                rdkit_dict_options[label] = value

        for ikey in var_to_dict_rdkit_boolean.keys():
            regex = ikey+".*="
            if re.search(regex, data) is not None:
                value = re.search(ikey+".*=.*", data).group(0)
                value = value.replace("'", "")
                value = value.replace("\"", "")
                value = value.split("=")[-1]
                if value.find('True') != -1:
                    valbol = True
                else:
                    valbol = False
                label = var_to_dict_rdkit_boolean[ikey]
                rdkit_dict_options[label] = valbol

        for ikey in var_to_dict_openbabel_intfloat.keys():
            regex = ikey+".*="
            if re.search(regex, data) is not None:
                value = re.search(ikey+".*=.*", data).group(0)
                value = value.split("=")[-1]
                label = var_to_dict_openbabel_intfloat[ikey]
                openbabel_dict_options[label] = value

        for ikey in var_to_dict_openbabel_boolean.keys():
            regex = ikey+".*="
            if re.search(regex, data) is not None:
                value = re.search(ikey+".*=.*", data).group(0)
                value = value.replace("'", "")
                value = value.replace("\"", "")
                value = value.split("=")[-1]
                if value.find('True') != -1:
                    valbol = True
                else:
                    valbol = False
                label = var_to_dict_openbabel_boolean[ikey]
                openbabel_dict_options[label] = valbol

        for ikey in var_to_dict_openbabel_str.keys():
            regex = ikey+".*="
            if re.search(regex, data) is not None:
                str2 = re.search(ikey+".*=.*", data).group(0)
                # value = re.search('(?:\'|").*(?:\'|")', str2).group(0)
                value = re.search("(?:['\"]).*(?:['\"])", str2).group(0)
                value = value.replace("'", "")
                value = value.replace("\"", "")
                label = var_to_dict_openbabel_str[ikey]
                openbabel_dict_options[label] = value

        for ikey in var_to_dict_extractmd_str.keys():
            regex = ikey+".*="
            if re.search(regex, data) is not None:
                str2 = re.search(ikey+".*=.*", data).group(0)
                # value = re.search('(?:\'|").*(?:\'|")', str2).group(0)
                value = re.search("(?:['\"]).*(?:['\"])", str2).group(0)
                value = value.replace("'", "")
                value = value.replace("\"", "")
                label = var_to_dict_extractmd_str[ikey]
                # print(label, ikey, var_to_dict_extractmd_str[ikey], value)
                extractmd_dict_options[label] = value

        for ikey in var_to_dict_pdbfolder_str.keys():
            regex = ikey+".*="
            if re.search(regex, data) is not None:
                str2 = re.search(ikey+".*=.*", data).group(0)
                # value = re.search('(?:\'|").*(?:\'|")', str2).group(0)
                value = re.search("(?:['\"]).*(?:['\"])", str2).group(0)
                value = value.replace("'", "")
                value = value.replace("\"", "")
                label = var_to_dict_pdbfolder_str[ikey]
                # print(label, ikey, var_to_dict_extractmd_str[ikey], value)
                pdbfolder_dict_options[label] = value

        for ikey in var_to_dict_extractmd_boolean.keys():
            regex = ikey+".*="
            if re.search(regex, data) is not None:
                value = re.search(ikey+".*=.*", data).group(0)
                value = value.replace("'", "")
                value = value.replace("\"", "")
                value = value.split("=")[-1]
                if value.find('True') != -1:
                    valbol = True
                else:
                    valbol = False
                label = var_to_dict_extractmd_boolean[ikey]
                extractmd_dict_options[label] = valbol

        for ikey in var_to_dict_extractmd_intfloat.keys():
            regex = ikey+".*="
            if re.search(regex, data) is not None:
                value = re.search(ikey+".*=.*", data).group(0)
                value = value.split("=")[-1]
                label = var_to_dict_extractmd_intfloat[ikey]
                extractmd_dict_options[label] = value

        for ikey in var_to_dict_systematicgrid_intfloat.keys():
            regex = ikey+".*="
            if re.search(regex, data) is not None:
                value = re.search(ikey+".*=.*", data).group(0)
                value = value.split("=")[-1]
                label = var_to_dict_systematicgrid_intfloat[ikey]
                systematicgrid_dict_options[label] = int(value)

        for ikey in var_to_dict_systematicgrid_list.keys():
            regex = ikey+".*="
            if re.search(regex, data) is not None:
                str2 = re.search(ikey+".*=.*", data).group(0)
                token_list = str2.split("=")[-1].lstrip()
                label = var_to_dict_systematicgrid_list[ikey]
                systematicgrid_dict_options[label] = ast.literal_eval(token_list)

        for ikey in var_to_dict_systematicgrid_boolean.keys():
            regex = ikey+".*="
            if re.search(regex, data) is not None:
                value = re.search(ikey+".*=.*", data).group(0)
                value = value.replace("'", "")
                value = value.replace("\"", "")
                value = value.split("=")[-1]
                if value.find('True') != -1:
                    valbol = True
                else:
                    valbol = False
                label = var_to_dict_systematicgrid_boolean[ikey]
                systematicgrid_dict_options[label] = valbol

        for ikey in var_to_dict_multipleline.keys():
            regex = ikey+".*="
            if re.search(regex, data) is not None:
                str2 = re.search(ikey+".*=.*", data).group(0)
                value = re.search("(?:['\"]).*(?:['\"])", str2).group(0)
                value = value.replace("'", "")
                value = value.replace("\"", "")
                token_list = value.split(",")
                line = ""
                for itoken in token_list:
                    line += itoken.lstrip()+"\n"
                label = var_to_dict_multipleline[ikey]
                window[label].update(value=line)


# =============================================================================
def import_pythonfileprops_to_gui(window, filename):

    # If filename None, the CANCEL button has been pushed
    if filename is None or len(filename) == 0:
        return False

    var_to_dict_str = {
        'p_nameserver': '-NAME_SERVER-',
        'p_username': '-USER_NAME-',
        'p_keysshfile': '-KEY_SSH_FILE-',
        'p_slurm_part': '-SLURM_PART-',
        'p_slurm_part_master': '-SLURM_PART_MASTER-',
        'p_g16path': '-GAUSSIAN16PACK-',
        'p_g16_keywords': '-KEYWORD_LINE-',
        'p_localdir_mol2conf': '-QM_PROP_MOL2LOCAL_DIR-',
        'p_localdir': '-QM_PROP_LOCAL_DIR-',
        'p_remotedir': '-QM_PROP_REMOTE_DIR-',
        'p_env_combo': '-ENV_COMBO-'
    }

    var_to_dict_boolean = {
        'p_run_gaussian': '-CHECKBOX_RUN_GAUSSIAN_OPT_PROP-'
    }

    var_to_encrypt = {'v_encrypt_pass': '-ENCRYPT_PASS_FILE-',
                      'v_node_master': '-NODE_MASTER-'
                      }

    var_to_dict_list = {
        'p_list_nodes': '-EXCLUDE_NODES-',
    }

    var_to_dict_list_multiline = {
        'p_bash_extrainfo': '-BASH_EXTRAINFO-'
    }

    var_to_dict_intfloat = {
        'p_ncpus': '-G16_NPROC-',
        'p_mem': '-G16_MEM-',
        'p_charge': '-CHARGE-',
        'p_multiplicity': '-MULTIPLICITY-',
        'p_maxjobsslurm': '-MAX_JOBS_SLURM-',
        'p_cutoffenergy': '-QM_PROP_CUTOFF_ENERGY-'
    }

    var_to_dict_pathlast = {
        'p_databasefullpath': '-INPUT_DATABASE_PROP-',
        'p_fileoutputfullpath': '-INPUT_LOG_PROP-'
    }

    with open(filename, 'r') as fpython:
        data = fpython.read()

        for ikey in var_to_dict_str.keys():
            regex = ikey+".*="
            if re.search(regex, data) is not None:
                str2 = re.search(ikey+".*=.*", data).group(0)
                # value = re.search("(?:'|\").*(?:'|\")", str2).group(0)
                value = re.search("(?:['\"]).*(?:['\"])", str2).group(0)
                value = value.replace("'", "")
                value = value.replace("\"", "")
                label = var_to_dict_str[ikey]
                window[label].update(value)

        for ikey in var_to_dict_intfloat.keys():
            regex = ikey+".*="
            if re.search(regex, data) is not None:
                value = re.search(ikey+".*=.*", data).group(0)
                value = value.split("=")[-1]
                label = var_to_dict_intfloat[ikey]
                window[label].update(value)

        for ikey in var_to_dict_pathlast.keys():
            regex = ikey+".*="
            if re.search(regex, data) is not None:
                str2 = re.search(ikey+".*=.*", data).group(0)
                # value = re.search('(?:\'|").*(?:\'|")', str2).group(0)
                value = re.search("(?:['\"]).*(?:['\"])", str2).group(0)
                value = value.replace("'", "")
                value = value.replace("\"", "")
                value = value.split("/")[-1]
                label = var_to_dict_pathlast[ikey]
                window[label].update(value)

        for ikey in var_to_dict_list.keys():
            regex = ikey+".*="
            if re.search(regex, data) is not None:
                str2 = re.search(ikey+".*=.*", data).group(0)
                # value = re.search('(?:\'|").*(?:\'|")', str2).group(0)
                value = re.search("(?:['\"]).*(?:['\"])", str2).group(0)
                value = value.replace("'", "")
                value = value.replace("\"", "")
                label = var_to_dict_list[ikey]
                window[label].update(value)

        for ikey in var_to_dict_list_multiline.keys():
            regex = ikey+".*="
            if re.search(regex, data) is not None:
                str2 = re.search(ikey+".*=.*", data).group(0)
                # value = re.search('(?:\'|").*(?:\'|")', str2).group(0)
                value = re.search("(?:['\"]).*(?:['\"])", str2).group(0)
                value = value.replace("'", "")
                lll = value.split(",")
                valuenew = ""
                for item in lll:
                    valuenew += item.strip()+"\n"
                label = var_to_dict_list_multiline[ikey]
                window[label].update(valuenew)

        for ikey in var_to_encrypt.keys():
            regex = ikey + ".*="
            if re.search(regex, data) is not None:
                str2 = re.search(ikey + ".*=.*", data).group(0)
                # value = re.search("(?:'|\").*(?:'|\")", str2).group(0)
                tmp = re.search("(?:['\"]).*(?:['\"])", str2)
                if tmp is not None:
                    value = tmp.group(0)
                    value = value.replace("'", "")
                    value = value.replace("\"", "")
                    label = var_to_encrypt[ikey]
                    window[label].update(value)

        for ikey in var_to_dict_boolean.keys():
            regex = ikey+".*="
            if re.search(regex, data) is not None:
                value = re.search(ikey+".*=.*", data).group(0)
                value = value.replace("'", "")
                value = value.replace("\"", "")
                value = value.split("=")[-1]
                if value.find('True') != -1:
                    valbol = True
                else:
                    valbol = False
                label = var_to_dict_boolean[ikey]
                window[label].update(valbol)

        # Extract parameters from the Keyword_line
        tokens = window['-KEYWORD_LINE-'].get().split()
        tokens_str = ' '.join(tokens)
        try:
            method, basisset = tokens[1].split("/")
        except ValueError:
            method = tokens[1]

        try:
            solvent_model = [string for string in tokens if "scrf" in string][0].split("(")[-1][0:-1]
        except (Exception,):
            solvent_model = 'None'

        try:
            solvent = [string for string in tokens if "solvent" in string][0].split("=")[-1][0:-1]
        except (Exception,):
            solvent = 'None'
        if solvent == '':
            solvent = 'None'

        window['-INPUT_METHOD_PROP-'].update(value=method)
        try:
            window['-INPUT_BASISSET_PROP-'].update(value=basisset)
        except UnboundLocalError:
            window['-INPUT_BASISSET_PROP-'].update(value="")
        window['-COMBO_MODELSOLVENT_PROP-'].update(value=solvent_model)
        window['-COMBO_SOLVENT_PROP-'].update(value=solvent)

        if tokens_str.count("freq") == 1 and tokens_str.count("wfn") == 0:
            window['-RADIO_FREQ-'].update(True)
            if tokens_str.count("scale") == 1:
                v = float([string for string in tokens if "scale" in string][0].split("=")[-1])
                window['-INPUT_SCALEVAL_FREQ-'].update(value=v, disabled=False, text_color='black')
            if tokens_str.count("anharm") == 1:
                window['-CHECKBOX_ANHAR_FREQ-'].update(value=True, disabled=False, text_color='black')
        elif tokens_str.count("wfn") == 1:
            window['-RADIO_ESPWFN-'].update(True)
        else:
            window['-RADIO_SP-'].update(True)

        window['-G16_KEYWORDS-'].update(window['-KEYWORD_LINE-'].get())

        window['-PROP_HIDEINPUTSCRIPT-'].update(value=filename)
