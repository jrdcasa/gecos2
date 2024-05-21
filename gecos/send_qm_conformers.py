import os
import utils
import datetime
from collections import defaultdict
import pandas as pd
import glob


# #################################################################################################################
def send_qm_conformers(nameserver, databasename, username, keyfile,
                       localdir, remotedir, g16exec, maxjobsslurm=50, regex="*g16*/*.com",
                       partition=None, partitionmaster=None, nodemaster=None,
                       exclude_nodes=None, ncpus=None, mem=None, timelimit=None, encrypted_pass=None,
                       extraslurminfo=None, logger=None):
    """

    Args:
        nameserver:
        databasename:
        username:
        keyfile:
        regex: Regular expression for the gaussian inputs
            (example: IsoP_???_gaussian.com)
        localdir:
        remotedir:
        g16exec:
        maxjobsslurm:
        partition:
        partitionmaster:
        nodemaster:
        exclude_nodes:
        ncpus:
        mem:
        timelimit:
        encrypted_pass (str): File containing the encrypted password usign key_file.pub
        extraslurminfo:
        logger:

    Returns:

    """

    m = "\t\t**************** OPTIMIZE CONFORMERS USING GAUSSIAN16 ***************"
    print(m) if logger is None else logger.info(m)

    # Server connection
    server = utils.ServerSlurmBasic(nameserver, databasename,
                                    username, keyfile, encrypted_pass=encrypted_pass, logger=None)
    server.connection(None)

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\t\t1. Connected to {} ({})\n".format(nameserver, now)
    m += "\t\t\t   Database {} created".format(databasename)
    print(m) if logger is None else logger.info(m)

    # List of the qm input files
    j = os.path.join(localdir, regex)
    com_list = sorted(glob.glob(j))

    # Send the QM calculations in the server
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\t\t2. Prepare SLURM script for each conformer ({}).".format(now)
    print(m) if logger is None else logger.info(m)
    utils.prepare_slurm_script_g16(com_list, g16exec, partition=partition,
                                   exclude_nodes=exclude_nodes,
                                   ncpus=ncpus, mem=mem, timelimit=timelimit, extraslurminfo=extraslurminfo)

    # Send the files to the server
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\t\t3. Send input files (com and sh) to server. ({})\n".format(now)
    m += "\t\t\t   Prepare and send full.sh script.\n".format(now)
    m += "\t\t\t   Update database.".format(now)
    print(m) if logger is None else logger.info(m)
    server.send_input_files_to_server(com_list, localdir, remotedir, maxjobsslurm)

    # Send calculations to slurm
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\t\t4. Send QM calculations to server. ({})".format(now)
    print(m) if logger is None else logger.info(m)
    server.send_qm_remote_calc(remotedir, partitionmaster, nodemaster, timelimit=timelimit, memlimit=mem)

    # Close server connection
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\t\t5. Close connection".format(now)
    print(m) if logger is None else logger.info(m)
    server.close_connection()

    m = "\t\t**************** OPTIMIZE CONFORMERS USING GAUSSIAN16 ***************\n"
    print(m) if logger is None else logger.info(m)


# ##################################################################################################################
def check_qm_jobs(nameserver, databasename, username, keyfile,
                  localdir, remotedir, outdir, pattern, exec_rmsddock, encrypted_pass=None,
                  energy_threshold=99999, cutoff_rmsd=0.0, logger=None):

    """

    Args:
        nameserver:
        databasename:
        username:
        keyfile:
        localdir:
        remotedir:
        outdir:
        pattern:
        exec_rmsddock:
        encrypted_pass (str): File containing the encrypted password usign key_file.pub
        energy_threshold:
        cutoff_rmsd:
        logger:

    Returns:

    """

    start = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")

    # Server connection
    server = utils.ServerSlurmBasic(nameserver, databasename,
                                    username, keyfile, encrypted_pass=encrypted_pass, logger=logger)
    server.connection(None)

    # Check qm jobs in the server
    completed = server.server_check_qm_jobs(localdir, remotedir, outdir, exec_rmsddock,
                                            energy_threshold=energy_threshold, cutoff_rmsd=cutoff_rmsd)
    nconfs = server._db_base.number_of_rows('qm_jobs')

    # Close server connection
    server.close_connection()
    if completed == -2:
        write_final_results(localdir, remotedir, server, pattern, energy_threshold, cutoff_rmsd, logger)
    elif completed == -1:
        return True
    else:
        m = "\t\t\t :::: Checking start: {}\n".format(start)
        m += "\t\t\t :::: Completed {} of {} jobs\n".format(completed, nconfs)
        end = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m += "\t\t\t :::: Checking end: {}\n\n".format(end)
        print(m) if logger is None else logger.info(m)
    return False


# ===========================================================================================
def write_final_results(localdir, remotedir, server, pattern, energy_threshold, cutoff_rmsd, logger=None):

    m = "\t\t**************** SUMMARY ***************\n"
    print(m) if logger is None else logger.info(m)

    # Convert database to panda's dataframe
    df = pd.read_sql_query("SELECT * from 'qm_jobs'", server._db_base._con)
    df = df.sort_values(["Cluster", "DeltaE"], ascending=[True, True])
    server._db_base.close_connection()

    # Files
    m = "\t\tLocaldir                      : {}\n".format(localdir)
    m += "\t\tDatabase                      : {}\n".format(server._db_base._dbname)
    m += "\t\tTrajectory file MM conformers : {}\n".format(pattern+"_conf_min_trj.pdb")
    m += "\t\tDirectory with QM files       : {}\n".format(pattern+"_g16_results/")
    print(m) if logger is None else logger.info(m)

    # Write table
    m = "\t\t{0:1s} {1:^35s} {2:^19s} {3:^8s} {4:^17s} {5:^10s}\n".format("#", "Name_job", "DeltaE_QM(kcal/mol)",
                                                                          "RMSD(A)", "RMSD_incluster(A)", "Cluster_QM")
    m += "\t\t# =============================================================================================\n"

    cd = defaultdict(list)
    for ind in df.index:
        try:
            line = "\t\t{0:^37s} {1:^19.2f} {2:^8.3f} {3:^17.3f} {4:^10d}\n".format(df['name_job'][ind],
                                                                                    df['DeltaE'][ind],
                                                                                    df['RMSD'][ind],
                                                                                    df['RMSD_INCLUSTER'][ind],
                                                                                    int(df['Cluster'][ind]))
            cd[df['Cluster'][ind]].append(df['name_job'][ind])
        except (Exception, ):
            line = "\t\t{0:^37s} {1:^19.2f} {2:^8.3f} {3:^17.3f} {4:^10d}\n".format(df['name_job'][ind],
                                                                                    999999.9,
                                                                                    -0.1,
                                                                                    -0.1,
                                                                                    -1)

        m += line
    print(m) if logger is None else logger.info(m)

    with open(os.path.join(localdir, "data_cluster.dat"), 'w') as f:
        line = "# CLuster  RMSD(A)   Delta-Energy(kcal/mol)\n"
        f.writelines(line)
        for idx in range(len(df)):
            try:
                line = "{0:6d}   {1:6d}   {2:10.3f}   {3:10.3f}\n".\
                    format(idx, int(df.Cluster[idx]), df.RMSD[idx], df.DeltaE[idx])
                f.writelines(line)
            except ValueError:
                continue
            except TypeError:
                continue
    m = "\t\t## Optimized conformers:\n"
    list_name_conformers = []
    for key, value in cd.items():
        m += "\t\t\t{}\n".format(value[0]+"_allign.mol2")
        list_name_conformers.append(value[0]+"_allign.mol2")
    m += "\t\t## Optimized conformers\n"

    m += "\n\t\t## Statistics QM calculations.\n"
    ll = [i for i in df.status_job if i == "COMPLETED"]
    nn = [i for i in df.status_job if i == "FAILED"]
    m += "\t\t   Total launched jobs    : {}\n".format(len(df))
    m += "\t\t   Total COMPLETED jobs   : {}\n".format(len(ll))
    m += "\t\t   Total FAILED    jobs   : {}\n".format(len(nn))
    m += "\t\t   Conformer energy range : 0.0 to {0:.2f} kcal/mol\n".format(float(max(df.DeltaE)))
    m += "\t\t   Total number of cluster: {}\n".format(len(cd))
    m += "\t\t       Energy threshold for cluster ({0:.2f} kcal/mol) and\n" \
         "\t\t       RMSD threshold ({1:.2f} kcal/mol) for clustering\n".format(energy_threshold, cutoff_rmsd)
    m += "\t\t## Statistics QM calculations.\n"

    # Write tcl to visualize the optimize conformers
    tclfullpath = os.path.join(localdir, pattern+"_g16_results", "QM_optimized_conformers.tcl")
    write_tcl_optimize_conformers(list_name_conformers, tclfullpath)

    m += "\n\t\t\t-------- ATTENTION --------\n"
    m += "\t\t\tYou should remove the files in {}\n".format(remotedir)
    m += "\t\t\tfrom the server {}\n".format(server._nameserver)
    m += "\t\t\t-------- ATTENTION --------\n\n"
    m += "\t\t**************** SUMMARY ***************\n"
    print(m) if logger is None else logger.info(m)

    end = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\n\t\t Job Finnished: {}\n".format(end)
    print(m) if logger is None else logger.info(m)


# ===========================================================================================
def write_tcl_optimize_conformers(list_mol2, pathfiletcl):

    lines = "proc newRep { sel type color rep imol} {\n"
    lines += "    mol selection $sel\n"
    lines += "    mol representation $type\n"
    lines += "    mol addrep $imol\n"
    lines += "    mol showrep $imol $rep on\n"
    lines += "    mol modcolor $rep $imol $color\n"
    lines += "}\n"
    lines += "\n"
    lines += "set dir \"{}\"\n".format(os.path.split(pathfiletcl)[0])
    lines += "\n"
    lines += "display projection orthographic\n"
    lines += "axes location off\n"
    lines += "color Display Background white\n"
    lines += "display depthcue off\n"
    lines += "\n"

    lines += "set listFiles {}\n"
    for item in list_mol2:
        lines += "lappend listFiles [list $dir/{}]\n".format(item)
    lines += "\n"

    lines += "\n"
    lines += "foreach ifile $listFiles {\n    mol addfile $ifile type mol2\n    set imol [molinfo top]\n}"
    lines += "\n"
    lines += 'set imol_ref [molinfo top]\n'
    lines += 'mol rename $imol_ref "OptimizedConformers"\n'
    lines += 'mol delrep 0 $imol_ref\n'
    lines += 'set rep1 0\n'
    lines += 'newRep "all" "CPK" "Name" $rep1 $imol_ref\n'
    lines += 'animate goto start\n'
    lines += '\n'

    with open(pathfiletcl, 'w') as ftcl:
        ftcl.writelines(lines)
