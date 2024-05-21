import argparse
import json
import os
import glob


# =============================================================================
def parse_arguments():

    desc = """ Generation of conformers (GeCos) """

    parser = argparse.ArgumentParser(description=desc)
    group1 = parser.add_mutually_exclusive_group(required=True)
    group1.add_argument("-j", "--json", dest="jsonfile",
                        help="A json file containing the parameters for a simulation using GeCos",
                        action="store", metavar="JSON_FILE")
    group1.add_argument("-p", "--python", dest="pythonfile",
                        help="A python script for running a simulation using GeCos",
                        action="store", metavar="PYTHON_FILE")

    args = parser.parse_args()

    return args


# =============================================================================
def read_json(filename):

    with open(filename, 'r') as fson:
        data = json.load(fson)

    filename_python = filename.split(".")[0]+".py"
    create_python_script(filename_python, data)

    return filename_python


# =============================================================================
def create_python_script(filename, keywords_dict, save=True):

    lines = ""
    lines += "import os\n"
    lines += "import utils\n"
    lines += "import gecos\n"
    lines += "\n"

    # KEYWORDS =====================================
    lines += "v_filename = '{}'\n".format(keywords_dict["moleculefile"])
    lines += "v_nameserver = '{}'\n".format(keywords_dict["nameserver"])
    lines += "v_username = '{}'\n".format(keywords_dict["username"])
    lines += "v_keysshfile = '{}'\n".format(keywords_dict["keyfile"])
    try:
        lines += "v_encrypt_pass = '{}'\n".format(keywords_dict["encrypted_passwd"])
    except KeyError:
        lines += "v_encrypt_pass = None\n"
    lines += "v_slurm_part = '{}'\n".format(keywords_dict["partition"])
    lines += "v_list_nodes = {}\n".format(keywords_dict["exclude_nodes"])
    lines += "v_slurm_part_master = '{}'\n".format(keywords_dict["partitionmaster"])
    lines += "v_node_master = '{}'\n".format(keywords_dict["nodemaster"])
    lines += "v_localdir = '{}'\n".format(keywords_dict["localdir"])
    lines += "v_remotedir = '{}'\n".format(keywords_dict["remotedir"])
    lines += "v_pattern = '{}'\n".format(keywords_dict["pattern"])
    lines += "v_databasefullpath = '{}'\n".format(os.path.join(keywords_dict["localdir"],
                                                               keywords_dict["databasename"]))
    lines += "v_fileoutputfullpath = '{}'\n".format(os.path.join(keywords_dict["localdir"],
                                                                 keywords_dict["filelog"]))
    lines += "v_g16path = '{}'\n".format(keywords_dict["exec_g16"])
    lines += "v_g16_keywords = '{}'\n".format(keywords_dict["g16_key"])
    lines += "v_ncpus = {0:d}\n".format(int(keywords_dict["g16_nproc"]))
    lines += "v_mem = {0:d}\n".format(int(keywords_dict["g16_mem"]))
    lines += "v_charge = {0:d}\n".format(int(keywords_dict["charge"]))
    lines += "v_multiplicity = {0:d}\n".format(int(keywords_dict["multiplicity"]))
    if keywords_dict["write_gaussian"]:
        lines += "v_write_gaussian = True\n"
    else:
        lines += "v_write_gaussian = False\n"
    lines += "v_nconfs = {0:d}\n".format(int(keywords_dict["nconfs"]))
    lines += "v_min_iter_mm = {0:d}\n".format(int(keywords_dict["minimize_iterations"]))
    lines += "v_cutoff_rmsd_qm = {0:.1f}\n".format(float(keywords_dict["cutoff_rmsd_QM"]))

    if keywords_dict["bond_perception"]:
        lines += "v_bond_perception = True\n"
    else:
        lines += "v_bond_perception = False\n"

    lines += "v_dockrmsdpack = '{}'\n".format(keywords_dict["exec_rmsddock"])
    lines += "v_confpack = '{}'\n".format(keywords_dict["conformer_program"])

    # RDKITS PARAMETERS ===================
    if keywords_dict["conformer_program"].upper() == "RDKIT":
        try:
            lines += "v_rdkit_maxattempts = {0:d}\n".format(int(keywords_dict["rdkit_maxattempts"]))
        except KeyError:
            lines += "v_rdkit_maxattempts = 1000\n"
        try:
            lines += "v_rdkit_prunermsthresh = {0:.3f}\n".format(float(keywords_dict["rdkit_prunermsthresh"]))
        except KeyError:
            lines += "v_rdkit_prunermsthresh = -0.01\n"
        try:
            if keywords_dict["rdkit_useexptorsionangleprefs"]:
                lines += "v_rdkit_useexptorsionangleprefs = True\n"
            else:
                lines += "v_rdkit_useexptorsionangleprefs = False\n"
        except KeyError:
            lines += "v_rdkit_useexptorsionangleprefs = True\n"
        try:
            if keywords_dict['rdkit_usebasicknowlwdge']:
                lines += "v_rdkit_usebasicknowlwdge = True\n"
            else:
                lines += "v_rdkit_usebasicknowlwdge = False\n"
        except KeyError:
            lines += "v_rdkit_usebasicknowlwdge = True\n"
        try:
            if keywords_dict["rdkit_enforcechirality"]:
                lines += "v_rdkit_enforcechirality = True\n"
            else:
                lines += "v_rdkit_enforcechirality = False\n"
        except KeyError:
            lines += "v_rdkit_enforcechirality = True\n"
        try:
            lines += "v_rdkit_cluster_method = '{}'\n".format(keywords_dict['rdkit_cluster_method'])
        except KeyError:
            lines += "v_rdkit_cluster_method = 'rmsd'\n"
        try:
            lines += "v_rdkit_ffname = '{}'\n".format(keywords_dict["rdkit_ffname"])
        except KeyError:
            lines += "v_rdkit_ffname = 'MMFF'\n"
        try:
            lines += "v_rdkit_cluster_thres = {}\n".format(float(keywords_dict["rdkit_cluster_thres"]))
        except KeyError:
            lines += "v_rdkit_cluster_thres = 1.0\n"
    # CONFAB PARAMETERS ===================
    elif keywords_dict["conformer_program"].upper() == "OPENBABEL":
        try:
            lines += "v_openbabel_rmsd_cutoff_confab = {0:f}\n".\
             format(float(keywords_dict["openbabel_rmsd_cutoff_confab"]))
        except KeyError:
            lines += "v_openbabel_rmsd_cutoff_confab = 0.5\n"
        try:
            lines += "v_openbabel_energy_cutoff_confab = {0:f}\n".\
                 format(float(keywords_dict["openbabel_energy_cutoff_confab"]))
        except KeyError:
            lines += "v_openbabel_energy_cutoff_confab = 50.0\n"
        try:
            if keywords_dict['openbabel_verbose']:
                lines += "v_openbabel_verbose = True\n"
            else:
                lines += "v_openbabel_verbose = False\n"
        except KeyError:
            lines += "v_openbabel_verbose = False\n"
        try:
            lines += "v_openbabel_rmsddock_confab = {0:f}\n".\
                format(float(keywords_dict["openbabel_rmsddock_confab"]))
        except KeyError:
            lines += "v_openbabel_rmsddock_confab = 2.0\n"
        try:
            lines += "v_openbabel_ffname = '{}'\n".format(keywords_dict["openbabel_ffname"])
        except KeyError:
            lines += "v_openbabel_ffname = 'MMFF'\n"
        try:
            lines += "v_openbabel_cluster_energy_threshold = {0:f}\n".\
                format(float(keywords_dict["openbabel_cluster_energy_threshold"]))
        except KeyError:
            lines += "v_openbabel_cluster_energy_threshold = 99999.0\n"
        try:
            lines += "v_openbabel_cluster_max_number_cluster = {0:d}\n".\
                format(int(keywords_dict['openbabel_cluster_max_number_cluster']))
        except KeyError:
            lines += "v_openbabel_cluster_max_number_cluster = 100\n"

    # LOOP
    lines += "\n"
    lines += "if not os.path.isfile(v_databasefullpath):\n\n"
    lines += "    log = utils.init_logger(\n" \
             "        \"Output2\",\n" \
             "        fileoutput=v_fileoutputfullpath,\n" \
             "        append=False, inscreen=False)\n\n"

    if keywords_dict["conformer_program"].upper() == "RDKIT":
        lines += "    g1 = gecos.GecosRdkit(\n" \
                 "        filename=v_filename,\n" \
                 "        total_charge=v_charge,\n" \
                 "        bond_perception=v_bond_perception,\n" \
                 "        logger=log)\n"
        lines += "\n"
        lines += "    g1.generate_conformers(\n" \
                 "        v_localdir,\n" \
                 "        nconfs=v_nconfs,\n" \
                 "        minimize_iterations=v_min_iter_mm,\n" \
                 "        maxattempts=v_rdkit_maxattempts,\n" \
                 "        prunermsthresh=v_rdkit_prunermsthresh,\n" \
                 "        useexptorsionangleprefs=v_rdkit_useexptorsionangleprefs,\n" \
                 "        usebasicknowledge=v_rdkit_usebasicknowlwdge,\n" \
                 "        enforcechirality=v_rdkit_enforcechirality,\n" \
                 "        ff_name=v_rdkit_ffname,\n" \
                 "        cluster_method=v_rdkit_cluster_method,\n" \
                 "        cluster_threshold=v_rdkit_cluster_thres,\n" \
                 "        write_gaussian=v_write_gaussian,\n" \
                 "        pattern=v_pattern,\n" \
                 "        g16_key=v_g16_keywords,\n" \
                 "        g16_nproc=v_ncpus,\n" \
                 "        g16_mem=v_mem,\n" \
                 "        charge=v_charge,\n" \
                 "        multiplicity=v_multiplicity)\n"
    elif keywords_dict["conformer_program"].upper() == "OPENBABEL":
        lines += "    g1 = gecos.GecosPyBabel(\n" \
                 "        filename=v_filename,\n" \
                 "        exec_rmsddock=v_dockrmsdpack,\n" \
                 "        total_charge=v_charge,\n" \
                 "        bond_perception=v_bond_perception,\n" \
                 "        logger=log)\n"
        lines += "\n"
        lines += "    g1.generate_conformers(\n" \
                 "        v_localdir,\n" \
                 "        nconfs=v_nconfs,\n" \
                 "        minimize_iterations=v_min_iter_mm,\n" \
                 "        rmsd_cutoff_confab=v_openbabel_rmsd_cutoff_confab,\n" \
                 "        energy_cutoff_confab=v_openbabel_energy_cutoff_confab,\n" \
                 "        confab_verbose_confab=v_openbabel_verbose,\n" \
                 "        cutoff_rmsddock_confab=v_openbabel_rmsddock_confab,\n" \
                 "        energy_threshold_cluster=v_openbabel_cluster_energy_threshold,\n" \
                 "        max_number_cluster=v_openbabel_cluster_max_number_cluster,\n" \
                 "        ff_name=v_openbabel_ffname,\n" \
                 "        pattern=v_pattern,\n" \
                 "        write_gaussian=v_write_gaussian,\n" \
                 "        g16_key=v_g16_keywords,\n" \
                 "        g16_nproc=v_ncpus,\n" \
                 "        g16_mem=v_mem,\n" \
                 "        charge=v_charge,\n" \
                 "        multiplicity=v_multiplicity\n" \
                 "        )\n"

        lines += "\n"
    else:
        msg = "Package {} to calculate conformers is not available.".format(keywords_dict["conformer_program"])
        print(msg)

    lines += "\n"
    lines += "    gecos.send_qm_conformers(" \
             "\n            v_nameserver," \
             "\n            v_databasefullpath," \
             "\n            v_username," \
             "\n            v_keysshfile," \
             "\n            v_localdir," \
             "\n            v_remotedir," \
             "\n            v_g16path,"\
             "\n            regex='*g16*/*.com'," \
             "\n            partition=v_slurm_part," \
             "\n            exclude_nodes=v_list_nodes," \
             "\n            ncpus=v_ncpus, " \
             "\n            partitionmaster=v_slurm_part_master," \
             "\n            nodemaster=v_node_master," \
             "\n            mem=v_mem," \
             "\n            encrypted_pass=v_encrypt_pass,"\
             "\n            logger=log)\n"
    lines += "\n"
    lines += "else:\n"
    lines += "\n"
    lines += "    log = utils.init_logger(" \
             "\n            \"Output2\"," \
             "\n            fileoutput=v_fileoutputfullpath," \
             "\n            append=True," \
             "\n            inscreen=False)\n"
    lines += "\n"
    lines += "    v_outdir = os.path.join(v_localdir, v_pattern + '_g16_results')\n"
    lines += "\n"
    lines += "    gecos.check_qm_jobs(" \
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
def main_app():

    opts = parse_arguments()

    # Read json and create a python script to run gecos
    filename_python = ""
    if opts.jsonfile:
        if opts.jsonfile.find(".json") != -1:
            filename_python = read_json(opts.jsonfile)
        else:
            print("File {} seems to be in a format different to json.".format(opts.jsonfile))
            exit()
    if opts.pythonfile:
        if opts.pythonfile.find(".py") != -1:
            filename_python = opts.pythonfile
        else:
            print("File {} seems to be not a python script.".format(opts.pythonfile))
            exit()

    local_dir = os.getcwd()
    try:
        database_name = glob.glob("*.db")[0]
    except IndexError:
        database_name = "None"
    print("Local_dir    : {}".format(local_dir))
    print("Database Name: {}".format(database_name))

    # Run Gecos
    fulldbpath = os.path.join(local_dir, database_name)
    donepath = os.path.join(local_dir, "done")
    if os.path.isfile(donepath):
        msg = "GeCos calculation seems to be finished."
        print(msg)
    else:
        if not os.path.isfile(fulldbpath):
            msg = "Database is not available. GeCos will run conformer search."
            print(msg)
        else:
            msg = "Database exists. GeCos will check the calculations."
            print(msg)

    os.system("python " + filename_python)


# =============================================================================
if __name__ == "__main__":
    main_app()
