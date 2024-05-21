import PySimpleGUI as Sg
import paramiko
import os
import sys
from collections import defaultdict
from socket import gaierror, error
from passwd_encrypt.passwd_encrypt import pw_encrypt_msg, pw_decrypt_msg


# =============================================================================
def popup_error(window, msg):
    loc = window.current_location()
    x, y = window.size
    newloc = (loc[0] + x / 2., loc[1] + y / 2.)
    # window.disappear()
    Sg.popup(msg, title='ERROR', grab_anywhere=True, location=newloc,
             background_color='red', text_color='white')
    # window.reappear()


# =============================================================================
def popup_msg(window, msg):
    loc = window.current_location()
    x, y = window.size
    newloc = (loc[0] + x / 2., loc[1] + y / 2.)
    # window.disappear()
    Sg.popup(msg, title='INFO', grab_anywhere=True, location=newloc,
             background_color='green', text_color='white', )
    # window.reappear()


# =============================================================================
def popup_msg_run(window, msg):

    loc = window.current_location()
    newloc = (loc[0] + 100., loc[1] + 100.)
    # window.disappear()
    Sg.popup(msg, title='Running', grab_anywhere=True, location=newloc,
             background_color='blue', text_color='white', non_blocking=True, keep_on_top=True)
    # window.reappear()


# =============================================================================
def check_files_localdir(window, files_localpath_string):

    """
    Check files in the localdir
    These arguments are defined in the General Inputs window of the GUI

    :arg:
        window --> PySimpleGUI Window object
        str_databases_outputs --> A list ['-MOLECULE_INPUT-', '-LOCAL_DIR-', '-DOCKRMSDPACK-']

    :return:
        A boolean
    """

    for ikey in files_localpath_string:
        pathfile = window[ikey].get()
        if not os.path.exists(pathfile) and ikey == "-DOCKRMSDPACK-":
            is_default_exe_found = False
            for isyspath in sys.path:
                path = os.path.join(isyspath, "thirdparty/dockrmsd.x")
                if os.path.isfile(path):
                    window[ikey].update(path)
                    is_default_exe_found = True
                    break
            if not is_default_exe_found:
                msg = "{}: {} does not exist\n".format(ikey, pathfile)
                msg += "and default dockrmsd cannot find in sys.path"
                popup_error(window, msg)
                return False
        elif not os.path.exists(pathfile):
            msg = "{}: {} does not exist".format(ikey, pathfile)
            popup_error(window, msg)
            return False

    return True


# =============================================================================
def check_server_stuffs(window, dict_options, pass_encrypted_file):

    """
    Check argument related to the server.
    These arguments are defined in the General Inputs window of the GUI

    :arg:
        window --> PySimpleGUI Window object
        keyfile --> from Key SSH file file
        pass_encrypted_file

    :return:
    A tuple --> (server_object, boolean)
    """

    nameserver = dict_options["nameserver"]
    username = dict_options["username"]
    keyfile = dict_options["keyfile"]

    # Check server stuffs =====================================================
    server = paramiko.SSHClient()
    server.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    try:
        msg_key = "Key found, but password authentication is needed."
        key = paramiko.RSAKey.from_private_key_file(keyfile)
    except paramiko.ssh_exception.SSHException:
        msg = "Not a valid RSA private key file:\nkeyfile {}\n  ".format(keyfile)
        popup_error(window, msg)
        return None, False
    except FileNotFoundError:
        msg = "Key file is not found"
        popup_error(window, msg)
        return None, False

    try:
        server.connect(nameserver, username=username, pkey=key, timeout=30)
    except (gaierror, paramiko.ssh_exception.NoValidConnectionsError):
        msg = "Unable to connect to {} ".format(nameserver)
        popup_error(window, msg)
        return None, False
    except (paramiko.ssh_exception.AuthenticationException, UnboundLocalError):
        msg = msg_key + "\n" + "Enter password to connect remote server:"
        loc = window.current_location()
        x, y = window.size
        newloc = (loc[0] + (x / 2.)*0.5, loc[1] + (y / 2.)*0.5)

        if not os.path.isfile(pass_encrypted_file):
            msg_pubkey = "Insert the path for a RSA public key:"
            public_key = Sg.popup_get_text(msg_pubkey, location=newloc)
            try:
                path_to_public_key = os.path.split(public_key)[0]
                password = Sg.popup_get_text(msg, location=newloc, password_char='*')
                pass_encryped_file = os.path.join(path_to_public_key, "passwd_encrypted.bin")
                pw_encrypt_msg(public_key, password, fout_name=pass_encryped_file)
                popup_msg(window, "Encripted password is stored in {}".format(pass_encryped_file))
                window['-ENCRYPT_PASS_FILE-'].update(pass_encryped_file)
            except (TypeError, AttributeError, FileNotFoundError):
                msg = "Something wrong with authentication. Check key path and/or username"
                popup_error(window, msg)
                return None, False
        else:
            pass_encryped_file = window['-ENCRYPT_PASS_FILE-'].get()
            password = pw_decrypt_msg(keyfile, pass_encryped_file)

        try:
            server.connect(nameserver, username=username, password=password, timeout=30)
        except paramiko.ssh_exception.AuthenticationException:
            msg = "Authentication problem with password:\nusername {}\n Try again".format(username)
            popup_error(window, msg)
            return None, False
        except error:
            msg = "Timeout to connect (30 seconds)\n"
            msg += "Server name {}\n".format(nameserver)
            popup_error(window, msg)
            return None, False

    except AttributeError:
        msg = "There is some problem with the key file ({})\n".format(key)
        msg += "Server name {}\n".format(nameserver)
        popup_error(window, msg)
        return None, False

    except error:
        msg = "Timeout to connect (30 seconds)\n"
        msg += "Server name {}\n".format(nameserver)
        popup_error(window, msg)
        return None, False

    return server, True


# =============================================================================
def check_partition_stuffs(window, dict_options, server):

    """
    Check partition name into the server.
    These arguments are defined in the General Inputs window of the GUI

    :arg:
        window --> PySimpleGUI Window object
        keyfile --> from Key SSH file file
        server --> An server object

    :return:
        A boolean
    """

    partition = dict_options["partition"]
    partitionmaster = dict_options["partitionmaster"]

    for ipart in [partition, partitionmaster]:
        try:
            command = "sinfo -s | egrep {}".format(ipart)
            stdin, stdout, stderr = server.exec_command(command, timeout=30)  # Non-blocking call
        except paramiko.ssh_exception.SSHException:
            msg = "SLURM Partition {} cannot be checked\n  ".format(ipart)
            msg += "Try again!!!!"
            popup_error(window, msg)
            return False

        out_txt = stdout.read().decode('utf8')
        if len(out_txt) < 1:
            msg = "Partition {} does not exist\n  ".format(ipart)
            popup_error(window, msg)
            return False

    return True


# =============================================================================
def check_nodemaster(window, dict_options, server):

    """
    Check partition name into the server.
    These arguments are defined in the General Inputs window of the GUI

    :arg:
        window --> PySimpleGUI Window object
        keyfile --> from Key SSH file file
        server --> An server object

    :return:
        A boolean
    """

    nodemaster = dict_options["nodemaster"]
    partitionmaster = dict_options["partitionmaster"]

    if len(nodemaster) < 1 or nodemaster is None:
        return True

    command = "scontrol show node {} | egrep 'not found'".format(nodemaster)
    stdin, stdout, stderr = server.exec_command(command, timeout=10)  # Non-blocking call
    out_txt = stdout.read().decode('utf8')
    err_txt = stderr.read().decode('utf8')
    if len(out_txt) > 1 or len(err_txt) > 1:
        msg = "Nodemaster {} does not exist\n in partition {}\n  ".format(nodemaster, partitionmaster)
        msg += "\nOutput txt: {}".format(out_txt)
        msg += "\nError  txt: {}".format(err_txt)
        popup_error(window, msg)
        return False

    return True


# =============================================================================
def check_remotedir(window, dict_options, server):

    """
    Check remote dir in the remote server.
    These arguments are defined in the General Inputs window of the GUI

    :arg:
        window --> PySimpleGUI Window object
        keyfile --> from Key SSH file file
        server --> An server object

    :return:
        A boolean
    """

    remotedir = dict_options["remotedir"]
    nameserver = dict_options["nameserver"]

    command = "ls -d {}".format(remotedir)
    stdin, stdout, stderr = server.exec_command(command, timeout=10)  # Non-blocking call
    err_txt = stderr.read().decode('utf8')
    if len(err_txt) > 1 or len(remotedir) == 0:
        msg = "Remotedir {} does not exist\n in server {}\n  ".format(remotedir, nameserver)
        msg += "\nError txt: {}".format(err_txt)
        popup_error(window, msg)
        return False

    return True


# =============================================================================
def check_g16_stuffs(window, dict_options, server):

    """
    Check path of g16 in the remote server and keywords length
    These arguments are defined in the General Inputs window of the GUI

    :arg:
        window --> PySimpleGUI Window object
        keyfile --> from Key SSH file file
        server --> An server object

    :return:
        A boolean
    """

    g16path = dict_options["g16path"]
    nameserver = dict_options["nameserver"]
    # Check gaussian keywords
    keywords = window['-G16_KEYWORDS-'].get()
    if len(keywords) < 5:
        msg = "G16 keywords seems to be incorrect\n  "
        popup_error(window, msg)
        return False

    # Check gaussian path
    command = "ls -d {}".format(g16path)
    stdin, stdout, stderr = server.exec_command(command, timeout=10)  # Non-blocking call
    err_txt = stderr.read().decode('utf8')
    if len(err_txt) > 1 or len(g16path) == 0:
        msg = "G16 path {} does not exist\n in server {}\n  ".format(g16path, nameserver)
        msg += "\nError txt: {}".format(err_txt)
        popup_error(window, msg)
        return False

    return True


# =============================================================================
def check_dblog_pattern(window, str_databases_outputs):

    """
    Check path of g16 in the remote server and keywords length
    These arguments are defined in the General Inputs window of the GUI

    :arg:
        window --> PySimpleGUI Window object
        str_databases_outputs --> A list ['-PATTERN-', '-DATABASE_NAME-', '-FILENAME_LOG-']
        server --> An server object

    :return:
        A boolean
    """

    for ikey in str_databases_outputs:
        string = window[ikey].get()
        if not string:
            msg = "Field {} cannot be empty\n  ".format(ikey)
            popup_error(window, msg)
            return False
    return True


# =============================================================================
def check_types(window, keys_input_int_labels, keys_input_float_labels):

    for ilabel in keys_input_int_labels:
        try:
            int(window[ilabel].get())
        except ValueError:
            msg = "{}: {} must be an integer.\n  ".format(ilabel, window[ilabel].get())
            popup_error(window, msg)
            return False

    for ilabel in keys_input_float_labels:
        try:
            float(window[ilabel].get())
        except ValueError:
            msg = "{}: {} must be a float.\n  ".format(ilabel, window[ilabel].get())
            popup_error(window, msg)
            return False

    return True


# =============================================================================
def check_prop_conformers(window, keys_properties, dict_properties):

    """
    Check property for conformers
    These arguments are defined in the Conformer Properties window of the GUI

    :arg:
        window --> PySimpleGUI Window object

    :return:
        A boolean
    """

    dict_properties['p_nameserver'] = str(window['-NAME_SERVER-'].get())
    dict_properties['p_username'] = str(window['-USER_NAME-'].get())
    dict_properties['p_keysshfile'] = str(window['-KEY_SSH_FILE-'].get())
    dict_properties['p_encrypt_pass'] = str(window['-ENCRYPT_PASS_FILE-'].get())
    dict_properties['p_slurm_part'] = str(window['-SLURM_PART-'].get())
    dict_properties['p_list_nodes'] = [window['-EXCLUDE_NODES-'].get()]
    dict_properties['p_slurm_part_master'] = str(window['-SLURM_PART_MASTER-'].get())
    dict_properties['p_node_master'] = str(window['-NODE_MASTER-'].get())
    dict_properties['p_g16path'] = str(window['-GAUSSIAN16PACK-'].get())
    dict_properties['p_ncpus'] = int(window['-G16_NPROC-'].get())
    dict_properties['p_mem'] = int(window['-G16_MEM-'].get())
    dict_properties['p_charge'] = int(window['-CHARGE-'].get())
    dict_properties['p_multiplicity'] = int(window['-MULTIPLICITY-'].get())
    dict_properties['p_remotedir'] = str(window['-QM_PROP_REMOTE_DIR-'].get())
    dict_properties['p_fileproplist'] = [window['-LISTBOX_MOL2-'].Values]
    dict_properties['p_dockrmsdpack'] = window['-DOCKRMSDPACK-'].get()
    dict_properties['p_bash_extrainfo'] = window["-BASH_EXTRAINFO-"].get().split("\n")
    dict_properties["p_env_combo"] = window["-ENV_COMBO-"].get()
    dict_properties["p_run_gaussian"] = window["-CHECKBOX_RUN_GAUSSIAN_OPT_PROP-"].get()

    for ikey in keys_properties:
        if ikey == '-QM_PROP_LOCAL_DIR-':
            d = window[ikey].get()
            if not os.path.isdir(d):
                msg = "{}: {} must exist in the local server.\n  ".format(ikey, d)
                popup_error(window, msg)
                return False
            else:
                dict_properties["p_localdir"] = str(window[ikey].get())
        elif ikey == '-QM_PROP_MOL2LOCAL_DIR-':
            d = window[ikey].get()
            if not os.path.isdir(d):
                msg = "{}: {} must exist in the local server.\n  ".format(ikey, d)
                popup_error(window, msg)
                return False
            else:
                dict_properties["p_localdir_mol2conf"] = str(window[ikey].get())
        elif ikey == '-KEYWORD_LINE-':
            dict_properties["p_g16_keywords"] = window[ikey].get()
            if len(dict_properties["p_g16_keywords"]) < 5:
                msg = "Gaussian 16 keywords are not correct.\n {} \n".format(str(dict_properties["p_g16_keywords"]))
                popup_error(window, msg)
                return False
        elif ikey == '-INPUT_DATABASE_PROP-':
            p = window['-INPUT_DATABASE_PROP-'].get()
            if len(p) < 1:
                msg = "{} cannot be empty.\n  ".format(ikey)
                popup_error(window, msg)
                return False
            d = os.path.join(window['-QM_PROP_LOCAL_DIR-'].get(), p)
            dict_properties["p_databasefullpath"] = os.path.join(str(d))
        elif ikey == '-INPUT_LOG_PROP-':
            p = window['-INPUT_LOG_PROP-'].get()
            if len(p) < 1:
                msg = "{} cannot be empty.\n  ".format(ikey)
                popup_error(window, msg)
                return False
            d = os.path.join(window['-QM_PROP_LOCAL_DIR-'].get(), p)
            dict_properties["p_fileoutputfullpath"] = os.path.join(str(d))

        elif ikey == '-DOCKRMSDPACK-':
            is_default_exe_found = False
            for isyspath in sys.path:
                path = os.path.join(isyspath, "thirdparty/dockrmsd.x")
                if os.path.isfile(path):
                    window[ikey].update(path)
                    is_default_exe_found = True
                    break
            if not is_default_exe_found:
                pathfile = window["-DOCKRMSDPACK-"].get()
                msg = "{}: {} does not exist\n".format(ikey, pathfile)
                msg += "and default dockrmsd cannot find in sys.path"
                popup_error(window, msg)
                return False

    return True


# =============================================================================
def check_prop_types(window, keys_input_int_labels):

    for ilabel in keys_input_int_labels:
        try:
            int(window[ilabel].get())
        except ValueError:
            msg = "{}: {} must be an integer.\n  ".format(ilabel, window[ilabel].get())
            popup_error(window, msg)
            return False

    return True


# =============================================================================
def check_extract_options(window, keys_extract):

    dict_extract = defaultdict()
    keys_extract = ['-MOLECULE_INPUT_EXTRACT-', '-METHOD_EXTRACT-', '-RADIUS_SPHERE-']

    for ikey in keys_extract:
        if ikey == '-MOLECULE_INPUT_EXTRACT-':
            # Check if file exists
            if not os.path.isfile(window[ikey].get()):
                msg = "File {} must be exist.\n  ".format(window[ikey].get())
                popup_error(window, msg)
                return None, False
            else:
                dict_extract['v_filemolextractfullpath'] = window[ikey].get()
                dirtmp = os.path.split(window[ikey].get())
                dict_extract['v_fileoutputfullpath'] = os.path.join(dirtmp[0], "extract_mol.log")
        elif ikey == '-METHOD_EXTRACT-':
            dict_extract['v_extractmethod'] = window[ikey].get()
        elif ikey == '-RADIUS_SPHERE-':
            ivalue = window[ikey].get()
            try:
                dict_extract['v_radiussphere'] = float(ivalue)
            except ValueError:
                msg = "Radius must be a real number.\n "
                popup_error(window, msg)
                return None, False

    return dict_extract, True
