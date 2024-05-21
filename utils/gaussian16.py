import glob
import os
import subprocess
import datetime
import numpy as np
from openbabel import openbabel as ob
from collections import defaultdict
from utils.atomic_data import atomic_number, element_vdw_truhlar_radius


# ===========================================================================================
def prepare_slurm_script_g16(listfiles, g16exec, partition=None, exclude_nodes=None,
                             ncpus=None, mem=None, timelimit=None, extraslurminfo=None):
    """

    Args:
        listfiles:
        g16exec: example: /opt/gaussian/g16_legacy
        partition:
        exclude_nodes:
        ncpus:
        mem:
        timelimit:
        extraslurminfo:

    Returns:

    """

    g16path = os.path.split(g16exec)[0]

    ll = "#!/bin/bash\n"
    if partition is not None:
        ll += "#SBATCH --partition={}\n".format(partition)
    if exclude_nodes is not None:
        lnodes = ""
        for item in exclude_nodes:
            lnodes += item + ", "
        lnodes = lnodes[:-2]
        ll += "#SBATCH --exclude=\"{}\"\n".format(lnodes)
    if ncpus is not None:
        ll += "#SBATCH --cpus-per-task={}\n".format(ncpus)
    if mem is not None:
        ll += "#SBATCH --mem={}M\n".format(mem)
    if timelimit is not None:
        ll += "#SBATCH --time={}\n".format(timelimit)

    # Send all files in listfiles from localdir to remotedir
    jobindex = 0
    for ifile in listfiles:
        base = os.path.splitext(ifile)[0]
        inputfile_com = os.path.split(ifile)[-1]
        ll2 = "#SBATCH --job-name={}\n".format(inputfile_com.split(".")[0])
        ll2 += "\n"

        try:
            for item in extraslurminfo:
                ll2 += "{}\n".format(item)
            ll2 += "g16 {}\n".format(inputfile_com)
            jobindex += 1
        except TypeError:
            ll2 += "g16legacy_root={}\n".format(g16path)
            ll2 += 'GAUSS_SCRDIR="$TMPDIR"\n'
            ll2 += "source $g16legacy_root/bsd/g16.profile\n"
            ll2 += "export g16legacy_root GAUSS_SCRDIR\n"

        with open(base + ".sh", 'w') as fin:
            fin.writelines(ll)
            fin.writelines(ll2)


# ===========================================================================================
def generate_bashscript_send_slurm(localdir, maxjobsslurm=50):

    """Generate the script **full_send.sh** in order to send jobs to a SLURM server

    Args:
        localdir (str): Path to store files in the local server.
        maxjobsslurm (str): Number of maximum jobs send to slurm

    """

    localfile = localdir+"/"+"full_send.sh"
    with open(localfile, 'w') as f:
        ll = '#!/bin/bash\n'
        ll += "\n"
        ll += "# NJOBS          --> Number of jobs sent to the slurm system\n"
        ll += "# MAXJOBSINSLURM --> Maximum number of jobs in the slurm system\n"
        ll += "# JOBSEND        --> Number of jobs finished in the slurm system\n"
        ll += "# TOTALJOBS      --> Jobs to be sent to the slurm system\n"
        ll += "# jobs.txt       --> Info of the jobs sent or finished in the slurm system\n"
        ll += '\n'
        ll += 'MAXJOBSINSLURM={}\n'.format(maxjobsslurm)
        ll += '\n'
        ll += 'NJOBS=`squeue -h |wc -ll`\n'
        ll += '\n'
        ll += 'COM=(`ls *.com`)\n'
        ll += 'LENGTH=${#COM[@]}\n'
        ll += '\n'

        ll += 'if [[ ! -e ./jobs.txt ]]; then\n'
        ll += '    echo -n >./jobs.txt\n'
        ll += 'fi\n'
        ll += '\n'
        ll += 'index=0\n'
        ll += 'while [ ${index} -lt ${LENGTH} ]; do\n'
        ll += '\n'
        ll += '    current=${COM[$index]}\n'
        ll += '\n'
        ll += '    if [[ $NJOBS -lt $MAXJOBSINSLURM ]]; then\n'
        ll += '        base="${COM[$index]%.*}"\n'
        ll += '        sbatch ${base}.sh  1 > tmp.txt\n'
        ll += '        jobid=`awk \'{print $NF}\' tmp.txt`\n'
        ll += '        echo "${jobid} ${base} ${base}.log" >>./jobs.txt\n'
        ll += '        rm tmp.txt\n'
        ll += '        index=`echo "$index+1" | bc -l`\n'
        ll += '        echo "NEW `date` --> JOBSEND: ${index}, TOTALJOBS: ${TOTALJOBS}, ${base}"\n'
        ll += '    else\n'
        ll += '        # Each 60 seconds checks the jobs\n'
        ll += '        sleep 60\n'
        ll += '        echo  "WAIT `date` --> JOBSEND: ${JOBSEND}, TOTALJOBS: ${TOTALJOBS}"\n'
        ll += '    fi\n'
        ll += '\n'
        ll += '    NJOBS=`squeue -h |wc -ll`\n'
        ll += '    TOTALJOBS=`ls -ld ./*_[0-9]*com |wc -ll`\n'
        ll += 'done\n'
        ll += '\necho \'Jobs Done!!!!\'\n'

        txt = ll

        f.writelines(txt)


# ===========================================================================================
def generate_bashscript_check_jobs(qm_engine, localdir, inputname="check_remote_dir.sh"):

    """Generate a bash script to check the jobs in the server. The name of the bash script is given by
    the parameter inputname

    Args:
        qm_engine (str):Name of the QM package (gaussian or nwchem)
        localdir (str): Path to store files in the local server.
        inputname (str):Name of the script

    """
    extract_energy = ""
    extract_time = ""
    if qm_engine.lower() == "gaussian":
        extract_energy = r"E\("
        extract_time = "`egrep 'Elapsed' $output  | awk '{print $3*86400+$5*3600+$7*60+$9}' | tail -1`\n"
        # is_complete_calc = "Normal term"
    elif qm_engine.lower() == "nwchem":
        extract_energy = "Total DFT energy"
        extract_time = "`egrep wall $output |tail -1|awk '{print $6}'|sed 's/.$//'`\n"
        # is_complete_calc = "Total times"
    elif qm_engine.lower() == "gamess":
        # is_complete_calc = "EXECUTION OF GAMESS TERMINATED"
        extract_energy = "TOTAL ENERGY ="
        extract_time = "`egrep 'TOTAL WALL CLOCK TIME' $idir/$output | tail -1 | awk '{print $5'}`"

    if localdir[-1] == "/":
        localfile = localdir + inputname
    else:
        localfile = localdir + "/" + inputname

    with open(localfile, 'w') as f:
        ll = '#!/bin/bash\n'
        ll += 'rm -f "summary.txt"\n'
        ll += 'rm -f "summary_energy.txt"\n'
        ll += 'echo -n >summary.txt\n'
        ll += 'echo -n >summary_energy.txt\n'
        ll += 'input="jobs.txt"\n'
        ll += 'index=0\n'
        ll += 'if [[ -f "$input" ]]; then\n'
        ll += '    while IFS= read -r line\n'
        ll += '    do\n'
        ll += '        pid=`echo "$line" | awk \'{print $1}\'`\n'
        ll += '        idir=`echo "$line" | awk \'{print $2}\'`\n'
        ll += '        output=`echo "$line" | awk \'{print $3}\'`\n'
        ll += '        p=`sacct -j $pid --format=JobID,JobName%60,State| egrep -v "batch|----|JobID" | head -1`\n'
        ll += '        [[ ! -z `echo $p` ]] && echo $p >>summary.txt\n'
        ll += '        if [[ ! -z `echo $p | egrep "COMPLETED"` ]]; then\n'
        if qm_engine.lower() == "gamess":
            ll += '        en=`egrep "{}"  $output | awk \'{{print $4}}\'`\n'.format(extract_energy)
        else:
            ll += '            en=`egrep "{}"  $output | tail -1 |awk \'{{print $5}}\'`\n'.format(extract_energy)
            ll += '            enmp2=`egrep "UMP2" $output | tail -1 | awk \'{{print $6}}\'|tr \'D\' \'E\'`\n'
            ll += '            [[ -z "$enmp2" ]] && en=$en || tmp=$enmp2\n'
            ll += '            [[ ! -z "$tmp" ]] && en=$(printf "%.14f" $tmp)\n'
        ll += '            time_s={}\n'.format(extract_time)
        ll += '            if [ $index -eq 0 ]; then\n'
        ll += '                 e0=$en\n'
        ll += '            fi\n'
        ll += '            index=`echo $index+1| bc -l`\n'
        ll += '            erel=`echo $e0 $en | awk \'{printf "%.3f", ($2-$1)*627.5095}\'`\n'
        ll += '            echo "$idir $pid $en $erel $time_s" >>summary_energy.txt\n'
        ll += '         #echo "$idir $pid $en $erel"\n'
        ll += '        fi\n'
        ll += '    done < "$input"\n'
        ll += '\n'
        ll += '    [[ ! -z `egrep "PENDING|RUNNING" summary.txt` ]]  || echo -n >done\n'
        ll += 'fi\n'
        ll += '#echo -n >done\n'
        txt = ll
        f.writelines(txt)


# ===========================================================================================
def get_optimized_coordinates(localdir):

    # All log files
    logfiles = sorted(glob.glob(os.path.join(localdir, "*.log")))

    # Open movie file
    fmovie = open(os.path.join(localdir, "Conformers_opt_noallign.xyz"), 'w')

    for ilog in logfiles:
        # Find last_index for each log file
        with open(ilog, 'r') as fin:
            lines = fin.readlines()
            index = 0
            for iline in lines:
                if iline.find("Input orientation:") != -1 or iline.find("Standard orientation") != -1:
                    index_last = index
                if iline.find("NAtoms=") != -1:
                    natoms = int(iline.split()[1])
                index += 1
        # Write xyz file for the optimized structure
        basename = ilog.split(".")[0]
        xyzname = basename+".xyz"
        with open(xyzname, 'w') as fout:
            fout.writelines("{}\n".format(natoms))
            fout.writelines("{}\n".format(xyzname))
            fmovie.writelines("{}\n".format(natoms))
            fmovie.writelines("{}\n".format(xyzname))

            jindex = index_last+5
            for iatom in range(0, natoms):
                ll = lines[jindex].split()
                fout.writelines("{} {} {} {}\n".format(ll[1], ll[3], ll[4], ll[5]))
                fmovie.writelines("{} {} {} {}\n".format(ll[1], ll[3], ll[4], ll[5]))
                jindex += 1

        # Create mol2 file
        obconversion = ob.OBConversion()
        obconversion.SetInAndOutFormats('xyz', 'mol2')
        obmol = ob.OBMol()
        obconversion.ReadFile(obmol, xyzname)
        obconversion.WriteFile(obmol, basename+".mol2")

    fmovie.close()


# ===========================================================================================
def cluster_optimized_coordinates(edict, localdir, exec_rmsddock, cutoff=1.0,
                                  energy_threshold=99999, maximum_number_clusters=100):

    """

    Args:
        edict:
        localdir:
        exec_rmsddock:
        cutoff:
        energy_threshold:
        maximum_number_clusters:

    Returns:

    """

    # Conversion
    hart2kcal = 627.5095

    # Calculate deltaE and RMSD
    isfirst = True
    mol1 = None
    rmsd_dict = defaultdict()
    deltaE_dict = defaultdict()
    for item in sorted(edict, key=edict.get, reverse=False):
        if isfirst:
            isfirst = False
            ob1 = ob.OBConversion()
            mol1 = ob.OBMol()
            ref_name_mol2 = os.path.join(localdir, item+".mol2")
            ob1.ReadFile(mol1, ref_name_mol2)
            namefile = os.path.join(localdir, item+"_allign.mol2")
            ob1.WriteFile(mol1, namefile)
            deltaE_dict[item] = 0.0
            emin = edict[item]
        else:
            ob2 = ob.OBConversion()
            mol2 = ob.OBMol()
            namefile = os.path.join(localdir, item+".mol2")
            ob2.ReadFile(mol2, namefile)
            # IncludeH = False; symmetry = False
            a = ob.OBAlign(False, False)
            a.SetRefMol(mol1)
            a.SetTargetMol(mol2)
            a.Align()
            a.UpdateCoords(mol2)
            tmp_name_mol2 = os.path.join(localdir, item+"_allign.mol2")
            ob2.WriteFile(mol2, tmp_name_mol2)
            deltaE_dict[item] = (edict[item] - emin)*hart2kcal

    # Clusterize
    # Building the clusters
    cluster_dict = defaultdict(dict)
    icluster = 0
    index = 0
    rmsd_dict_incluster = defaultdict()
    # deltaE_dict is sorted from lowest to hightest energies.
    nitems = len(deltaE_dict)
    for item in deltaE_dict:
        if index % (int(np.floor(nitems*0.1))+1) == 0:
            now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
            print("Clustering element {} out of {} ({})".format(index, nitems, now))

        energy = deltaE_dict[item]

        if icluster == 0:
            threshold = energy + energy_threshold
            icluster += 1
            cluster_dict[icluster] = {"seed": index, "lowest_energy": energy, "highest_energy": energy,
                                      "nelements": 0, "pairs": [], "files": [], "pattern": []}
            cluster_dict[icluster]["pairs"].append([0.000, energy])
            cluster_dict[icluster]["files"].append(item+"_allign.mol2")
            cluster_dict[icluster]["pattern"].append(item)
            cluster_dict[icluster]["nelements"] += 1
            rmsd_dict[item] = 0.000
            rmsd_dict_incluster [item] = 0.000
            min_energy = energy

        # elif energy < threshold:
        else:
            tmp_name_mol2 = os.path.join(localdir, item+"_allign.mol2")
            energy_current = deltaE_dict[item]
            found = False
            for i in range(1, icluster + 1):
                idx = cluster_dict[i]["seed"]
                ref_name_mol2 = os.path.join(localdir, cluster_dict[i]["files"][0])
                energy_ref = deltaE_dict[cluster_dict[i]["pattern"][0]]
                # Take into account Hs
                # cmd = '{} {} {} -h -c -s'.format(exec_rmsddock, ref_name_mol2, tmp_name_mol2)
                # Do not take into account Hs
                cmd = '{} {} {} -c -s'.format(exec_rmsddock, ref_name_mol2, tmp_name_mol2)
                # RMSD with the seed of the cluster[i]
                try:
                    rmsd_dock = float(subprocess.check_output(cmd, shell=True).decode())
                except subprocess.CalledProcessError:
                    try:
                        cmd = '{} {} {} -s'.format(exec_rmsddock, ref_name_mol2, tmp_name_mol2)
                        rmsd_dock = float(subprocess.check_output(cmd, shell=True).decode())
                    except subprocess.CalledProcessError:
                        # This exception can occur if the Template and query don't have the same bonding network
                        # so rmsd_dock cannot be calculated. Arbitrarily a value of 1000 is assigned.
                        rmsd_dock = 1000
                rmsd_dict_incluster[item] = rmsd_dock
                # Calculate the rmsd with the seed of the cluster 0 (lowest energy)
                min_name_mol2 = os.path.join(localdir, cluster_dict[1]["files"][0])
                # Do not take into account Hs. -c --> Assume atomic correspondence between files
                try:
                    cmd = '{} {} {} -c -s'.format(exec_rmsddock, min_name_mol2, tmp_name_mol2)
                    rmsd_dict[item] = float(subprocess.check_output(cmd, shell=True).decode())
                except subprocess.CalledProcessError:
                    try:
                        cmd = '{} {} {} -s'.format(exec_rmsddock, min_name_mol2, tmp_name_mol2)
                        rmsd_dict[item] = float(subprocess.check_output(cmd, shell=True).decode())
                    except subprocess.CalledProcessError:
                        # This exception can occur if the Template and query don't have the same bonding network
                        # so rmsd_dock cannot be calculated. Arbitrarily a value of 1000 is assigned.
                        rmsd_dict[item] = 1000
                if float(rmsd_dock) < cutoff and abs(energy_current-energy_ref) < energy_threshold:
                    cluster_dict[i]["pairs"].append([rmsd_dock, energy])
                    cluster_dict[i]["files"].append(item+"_allign.mol2")
                    cluster_dict[i]["pattern"].append(item)
                    cluster_dict[i]["highest_energy"] = energy
                    cluster_dict[i]["nelements"] += 1
                    found = True
                    break

            if not found:
                if icluster < maximum_number_clusters:
                    icluster += 1
                    cluster_dict[icluster] = {"seed": index, "lowest_energy": energy, "highest_energy": energy,
                                              "nelements": 0, "pairs": [], "files": [], "pattern": []}
                    cluster_dict[icluster]["files"].append(item+"_allign.mol2")
                    cluster_dict[icluster]["pattern"].append(item)
                    cluster_dict[icluster]["nelements"] += 1
                    cluster_dict[icluster]["pairs"].append([rmsd_dock, energy])
                    rmsd_dict_incluster[item] = 0.000
                else:
                    if cluster_dict[maximum_number_clusters + 1]:
                        cluster_dict[maximum_number_clusters + 1]["pairs"].append([rmsd_dock, energy])
                        cluster_dict[maximum_number_clusters + 1]["files"].append(item+"_allign.mol2")
                        cluster_dict[maximum_number_clusters + 1]["pattern"].append(item)
                        cluster_dict[maximum_number_clusters + 1]["highest_energy"] = energy
                        cluster_dict[maximum_number_clusters + 1]["nelements"] += 1
                        rmsd_dict_incluster[item] = 0.000
                    else:
                        cluster_dict[maximum_number_clusters + 1] = {"seed": index, "lowest_energy": energy,
                                                                     "highest_energy": energy,
                                                                     "nelements": 0, "pairs": [],
                                                                     "files": [], "pattern": []}
                        cluster_dict[maximum_number_clusters + 1]["files"].append(item+"_allign.mol2")
                        cluster_dict[maximum_number_clusters + 1]["pattern"].append(item+"_allign.mol2")
                        cluster_dict[maximum_number_clusters + 1]["nelements"] += 1
                        cluster_dict[maximum_number_clusters + 1]["pairs"].append([rmsd_dock, energy])
                        rmsd_dict_incluster[item] = 0.000
        # else:
        #     pass

        index += 1

    return cluster_dict, deltaE_dict, rmsd_dict, rmsd_dict_incluster


# ===========================================================================================
def write_gaussian_from_mol2(fileproplist, dirmol2, localdir, pattern="QM",
                             g16_key="#p 6-31g* mp2", g16_mem=4000,
                             g16_nproc=4, charge=0, multiplicity=1, logger=None):

    idx = 0
    for imol2 in fileproplist:
        fmol2path = os.path.join(dirmol2, imol2)

        obconversion = ob.OBConversion()
        obconversion.SetInAndOutFormats('mol2', 'xyz')

        obmol = ob.OBMol()
        xyz_obmol = obconversion.ReadFile(obmol, fmol2path)
        print(fmol2path)
        if not xyz_obmol:
            m = "Something is wrong with Gaussian16 input files 2\n"
            m += "Gaussian input files are not written"
            print(m) if logger is None else logger.error(m)
            return None

        xyz_list = obconversion.WriteString(obmol).split("\n")

        try:
            iconf_seed = imol2.split(".")[0]
        except:
            iconf_seed = idx

        fname = os.path.join(localdir, "{0:s}_{1:s}_g16.com".format(iconf_seed, pattern))
        fchk_name = "{0:s}_{1:s}_gaussian.chk".format(iconf_seed, pattern)

        with open(fname, 'w') as f:
            f.writelines("%chk={}\n".format(fchk_name))
            f.writelines("%nproc={}\n".format(g16_nproc))
            f.writelines("%mem={}Mb\n".format(g16_mem))
            f.writelines("{}\n".format(g16_key))
            f.writelines("\nConformer name {0:s}.\n".format(iconf_seed))
            f.writelines("\n")
            f.writelines("{0:1d} {1:1d}\n".format(charge, multiplicity))
            for line in xyz_list[2:]:
                f.writelines(line + "\n")
            if pattern.upper() == "WFN":
                # Write Connectivity matrix ===============================
                nbonds = obmol.NumBonds()
                connect = defaultdict(list)
                bondorder = defaultdict(list)
                atomic_numbers_set = []
                for ibond in range(nbonds):
                    i_obbond = obmol.GetBondById(ibond)
                    iatom = i_obbond.GetBeginAtomIdx()
                    jatom = i_obbond.GetEndAtomIdx()
                    iatom_atnum = i_obbond.GetBeginAtom().GetAtomicNum()
                    jatom_atnum = i_obbond.GetEndAtom().GetAtomicNum()
                    atomic_numbers_set.append(iatom_atnum)
                    atomic_numbers_set.append(jatom_atnum)
                    bo = i_obbond.GetBondOrder()
                    connect[iatom].append(jatom)
                    bondorder[iatom].append(bo)
                    connect[jatom].append(iatom)
                    bondorder[jatom].append(bo)
                for iat, value in sorted(connect.items()):
                    for jat in value:
                        connect[jat].remove(iat)

                atomic_numbers_set = set(atomic_numbers_set)

                line = ""
                for ikey, value in sorted(connect.items()):
                    line += "{0:4d} ".format(ikey)
                    jdx = 0
                    for item in value:
                        line += "{0:4d} {1:3.1f} ".format(item, bondorder[ikey][jdx])
                        jdx += 1
                    line += "\n"
                f.writelines(line + "\n")

                # Write radii
                line = ""
                for item in atomic_numbers_set:
                    for key, value in atomic_number.items():
                        if value == item:
                            line += "{0:s} {1:4.2f}\n".format(key, element_vdw_truhlar_radius[key])
                f.writelines(line + "\n")

                # Write wfn file name
                fwfn_name = "{0:s}_{1:s}_g16.wfn".format(iconf_seed, pattern)
                f.writelines(fwfn_name + "\n")
                f.writelines("\n")

        idx += 1
