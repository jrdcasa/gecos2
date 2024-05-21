import pandas as pd
import subprocess
import datetime
import os
import utils
from _collections import defaultdict
from openbabel import openbabel as ob


class GecosPyBabel:

    """
    This class computes conformers of a molecule using the algorithm ConFab [#]_
    implemeted in openbabel


    ** References **
    .. [#] Boyle et al, Confab - Systematic generation of diverse low-energy conformers", Journal of Cheminformatics 3,
     Article number: 8 (2011), https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-3-8
    """

    # =========================================================================
    def __init__(self, filename, exec_rmsddock, total_charge=0, bond_perception=False, logger=None):

        """

        Args:
            filename (str): Path to the input file
            total_charge (int) : Total charge of the molecule.
            bond_perception (bool): If True the bonds are calculated using the ``indigo-bondorders`` software
                (located in thirdparty/indigo-bondorder).
            logger (logger): An instance of the Logger class

        TODO:
            Implementation of the bond_perception with indigo.

        """

        self._logger = logger
        self._charge = total_charge
        self._mol_pybabel = None
        self._energy_conf = []
        self._isheader_print = False
        self._confs_to_write = 0
        self._exec_rmsddock = exec_rmsddock

        #self._df_conformers = pd.DataFrame(columns=['Conformations', 'IniEnergy', 'OptEnergy', 'RMSD', 'Cluster'])
        self._df_conformers = None

        file_format = filename.strip().split('.')[-1]

        # Report
        if not self._isheader_print:
            utils.print_header(self._logger)
            self._isheader_print = True

        m = "\t\t******** CREATION OF THE GECOS PYBABEL OBJECT ********\n\n"
        m += "\t\tConformers will be generated with openbabel library" \
             " (https://openbabel.org/wiki/Python, https://openbabel.github.io/api/3.0/index.shtml)\n"
        m += "\n"

        obconversion = ob.OBConversion()
        obconversion.SetInFormat(file_format)
        obmol = ob.OBMol()

        mols = []
        mols_in_file = obconversion.ReadFile(obmol, filename)
        while mols_in_file:
            mols.append(obmol)
            obmol = ob.OBMol()
            mols_in_file = obconversion.Read(obmol)

        # Only the first molecule is read in the `filename`
        # if there are more than one molecule in the input
        try:
            self._mol_pybabel = mols[0]
        except IndexError:
            self._mol_pybabel = None

        if self._mol_pybabel is not None:
            m += "\t\tMolecular file seed  : {}\n".format(filename)
            m += "\t\tBond Perception: {}\n".format(bond_perception)
            m += "\n"
        else:
            m += "\t\t\t WARNING!!!! Molecule might not be correctly setup.\n"
            m += "\t\t\t WARNING!!!! Openbabel molecule seems to be empty (self._mol_pybabel is None).\n"
            m += "\t\t\t WARNING!!!! conformers cannot be generated.\n"
            m += "\t\t\t WARNING!!!! Revise both input parameters and/or {} inputfile .\n".format(filename)
            m += "\n"

        m += "\t\tFormula          = {}\n".format(self._mol_pybabel.GetFormula())
        m += "\t\tNumber of atoms  = {}\n".format(self._mol_pybabel.NumAtoms())
        m += "\t\tMolecular weight = {0:.2f} g/mol\n".format(self._mol_pybabel.GetMolWt())
        m += "\t\tRotable bonds    = {} (does not include CH3 groups)\n".format(self._mol_pybabel.NumRotors())
        m += "\t\tNumber of rings  = {}\n".format(len(self._mol_pybabel.GetSSSR()))
        m += "\t\t******** CREATION OF THE GECOS PYBABEL OBJECT ********\n\n"
        print(m) if self._logger is None else self._logger.info(m)

        pass

    # =========================================================================
    def _mm_calc_energy(self, obmol=None, ff_name="MMFF"):

        if ff_name.upper() == "MMFF":
            ff = ob.OBForceField.FindForceField("mmff94")
        elif ff_name.upper() == "UFF":
            ff = ob.OBForceField.FindForceField("uff")
        elif ff_name.upper() == "GAFF":
            ff = ob.OBForceField.FindForceField("gaff")

        if obmol is None:
            ff.Setup(self._mol_pybabel)
        else:
            ff.Setup(obmol)
        energy = ff.Energy()

        return energy

    # =========================================================================
    def generate_conformers(self, localdir, nconfs=100000, minimize_iterations=0,
                            rmsd_cutoff_confab=0.5, energy_cutoff_confab=50.0,
                            confab_verbose_confab=False, cutoff_rmsddock_confab=2.0,
                            energy_threshold_cluster=99999, max_number_cluster=100, ff_name="MMFF",
                            write_gaussian=True, g16_key="#p 6-31g* mp2", g16_mem=4000,
                            g16_nproc=4, pattern="conformers", charge=0, multiplicity=1,
                            g16_extra_info="", isclustering=True, debug_flag=False):

        m = "\t\t**************** GENERATE CONFORMERS ***************\n"
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m += "\t\t 1. Generating conformers with openbabel ({})".format(nconfs, now)
        print(m) if self._logger is None else self._logger.info(m)
        if ff_name.upper() == "MMFF":
            ff = ob.OBForceField.FindForceField("mmff94")
        elif ff_name.upper() == "UFF":
            ff = ob.OBForceField.FindForceField("uff")
        elif ff_name.upper() == "GAFF":
            ff = ob.OBForceField.FindForceField("gaff")

        # Check if the force field is working with the molecule.
        try:
            assert (ff.Setup(self._mol_pybabel))
        except AssertionError:
            m = "\t\t Force field {} have problems with this molecule\n".format(ff_name)
            m += "\t\t Try other force field: GAFF, MMFF or UFF\n".format(ff_name)
            m += "\t\t PROGRAM STOPS!!!!!"
            print(m) if self._logger is None else self._logger.info(m)
            exit()

        nrotors = self._mol_pybabel.NumRotors()
        m = "\t\t\tMolecule = {}\n".format(self._mol_pybabel.GetTitle())
        m += "\t\t\tNumber of rotatable bonds = {}\n".format(str(nrotors))
        m += "\t\t\tForce Field = {}\n".format(ff.Description())
        m += "\t\t\tConformer Searching = {}".format("ConFab")
        print(m) if self._logger is None else self._logger.info(m)

        # Algorithm ConFab in:
        # Boyle et al, Journal of Cheminformatics 3, Article number: 8 (2011)
        # "Confab - Systematic generation of diverse low-energy conformers"
        # https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-3-8
        ff.DiverseConfGen(rmsd_cutoff_confab, nconfs, energy_cutoff_confab, confab_verbose_confab)
        ff.GetConformers(self._mol_pybabel)

        # No include the input conformation in the output
        self._confs_to_write = self._mol_pybabel.NumConformers() - 1

        m = "\t\t\tOptions: rmsd_cutoff = {} A, nconfs = {}\n".format(rmsd_cutoff_confab, nconfs)
        m += "\t\t\tOptions: energy_cutoff = {} kcal/mol, confab_verbose = {}\n".\
            format(energy_cutoff_confab, confab_verbose_confab)
        m += "\t\t\t ** CONFAB has generated {} conformers **".format(self._confs_to_write)
        print(m) if self._logger is None else self._logger.info(m)

        # Calculate mm energy, optimize and write noopt and opt conformers to two files.
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\t\t 2. Minimizing {} conformers. " \
            "Max_iters = {} ({})\n".format(self._confs_to_write, minimize_iterations, now)
        if minimize_iterations > 0:
            m += "\t\t\t Using {} forcefield".format(ff_name)
        else:
            m += "\t\t\t Minimization is not performed."
            minimize_iterations = 0
        print(m) if self._logger is None else self._logger.info(m)

        optimized_strings_mols = self._mm_calc_energy(ff, localdir, pattern=pattern,
                                                      out_format='mol2',
                                                      minimize_iterations=minimize_iterations, debug_flag=debug_flag)

        # Clustering MM conformers
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\t\t 3. Cluster Conformers ({})".format(now)
        print(m) if self._logger is None else self._logger.info(m)
        if isclustering:
            MMclusters = self._cluster_conformers(optimized_strings_mols, pattern=pattern,
                                                  cutoff_rmsddock=cutoff_rmsddock_confab,
                                                  energy_threshold=energy_threshold_cluster,
                                                  maximum_number_clusters=max_number_cluster,
                                                  ff_name=ff_name)
        else:
            MMclusters = defaultdict(dict)
            fname = "{}_info_conformers.dat".format(pattern)
            with open(fname, 'w') as f:
                f.write(self._df_conformers.to_string())

        # DEBUG write all conformers ==========================
        if debug_flag:
            filepdb_name = os.path.join(localdir, "{}_allconf_ob_nomin_trj.pdb".format(pattern))
            self.write_all_conformers_to_pdb(filepdb_name)
        # DEBUG write all conformers ==========================

        m = "\t\t\t{} conformers generated grouped in {} clusters (rmsd_threshold= {} angstroms)".\
            format(self._confs_to_write, len(MMclusters), cutoff_rmsddock_confab)
        print(m) if self._logger is None else self._logger.info(m)

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\t\t 4. Get structure of minimum energy for each cluster ({})".format(now)
        print(m) if self._logger is None else self._logger.info(m)

        filemol2_name = os.path.join(localdir, "{}_cluster_ob_optMM_trj.mol2".format(pattern))
        lowest_cluster_strings_mols = self._write_min_cluster_conformers_to_mol2(filemol2_name, MMclusters,
                                                                                 optimized_strings_mols)
        # DEBUG write all conformers ==========================
        if debug_flag:
            self._write_allstructures_clusters_to_mol2(MMclusters, optimized_strings_mols)
        # DEBUG ===============================================

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\t\t 5. Write Conformers to MOL2 ({})\n".format(now)
        m += "\t\t\tFilename: {}".format(filemol2_name)
        print(m) if self._logger is None else self._logger.info(m)

        if write_gaussian:
            self._write_gaussian(localdir, optimized_strings_mols, MMclusters, pattern=pattern, g16_key=g16_key,
                                 g16_mem=g16_mem, g16_nproc=g16_nproc, g16_extra_info=g16_extra_info, charge=charge, multiplicity=multiplicity)

        m = "\t\t**************** GENERATE CONFORMERS ***************\n"
        print(m) if self._logger is None else self._logger.info(m)

        m = "\n\t\t**************** CONFORMER MM ENERGIES ***************\n"
        m += "\n\t\tEnergies in kcal/mol, RMSD in angstroms."
        m += "\n\t\tRMSD threshold = {0:5.3f} angstroms.".format(rmsd_cutoff_confab)
        m += "\n\t\tNumber of clusters = {0:4d}".format(len(MMclusters))
        print(m) if self._logger is None else self._logger.info(m)
        m = "\t\tCluster Conformer_ID  {}_energy  Relative_energy " \
            "Highest_energy    RMSD   nelements\n".format(ff_name)
        separator = "\t\t" + len(m) * "-" + "\n"
        m += separator
        icluster = 1
        minEnergy = MMclusters[1]['lowest_energy']
        for icluster, ival in MMclusters.items():

            iconf = ival['seed']
            energy_abs = ival['lowest_energy']
            m += "\t\t {0:^5d}  {1:^12d}  {2:^12.2f}  {3:^12.2f}  {4:^14.2f}  {5:^10.3f}  {6:^6d}\n". \
                format(icluster, iconf, energy_abs,
                       energy_abs - minEnergy, MMclusters[icluster]['highest_energy'],
                       MMclusters[icluster]['pairs'][0][0],
                       MMclusters[icluster]['nelements'])
        m += separator
        m += "\t\t**************** CONFORMER MM ENERGIES ***************\n"
        print(m) if self._logger is None else self._logger.info(m)

        return True

    # =========================================================================
    @staticmethod
    def _write_allstructures_clusters_to_mol2(MMclusters, optimized_strings):

        folder = "full_mm_ob_optMM_clusters"

        try:
            os.mkdir(folder)
        except FileExistsError:
            import shutil
            shutil.rmtree(folder)
            os.mkdir(folder)

        for key, value in MMclusters.items():
            fname = os.path.join(folder, "cluster_{0:05d}.mol2".format(key))
            with open(fname, "a") as f:
                for item in value["elements"]:
                    imol = optimized_strings[item]
                    f.write(imol)

    # =========================================================================
    @staticmethod
    def _write_min_cluster_conformers_to_mol2(filename,
                                              MMcluster, optimized_strings_mols):

        """
        Write the conformer of minimum energy for each cluster in PDB trajectory

        Args:
            filename (str): Name of the PDB file
            MMcluster (defaultdict): A dict containing info og the MM clusterization
            optimized_strings_mols (list): List containing the mol2 string of all MM optimized conformers.

        Returns:
            A list containing the index of the lowest-energy conformer for each cluster

        """

        conf_id_list = []
        lowest_cluster_strings_mols = []

        with open(filename, 'w') as f:
            for key, value in MMcluster.items():
                iseed = value['seed']
                imol = optimized_strings_mols[iseed]
                lowest_cluster_strings_mols.append(imol)
                f.write(imol)

        return lowest_cluster_strings_mols

    # =========================================================================
    def _mm_calc_energy(self, ff, localdir, pattern="conformer", out_format="xyz",
                        minimize_iterations=0, debug_flag=False):

        # WRITE conformers without optimization ===================================================
        obconversion = ob.OBConversion()
        obconversion.SetOutFormat(out_format)
        nooptimized_strings_mol = []

        for iconf in range(self._confs_to_write):
            self._mol_pybabel.SetConformer(iconf)
            nooptimized_strings_mol.append(obconversion.WriteString(self._mol_pybabel))

        fnameout_noopt = "{}_allconf_ob_noopt_trj.{}".format(pattern, out_format)

        with open(fnameout_noopt, 'w') as f:
            for input_string in nooptimized_strings_mol:
                f.write(input_string)

        # READ conformers without optimizations and store them in mols_opt ========================
        obconversion = ob.OBConversion()
        obconversion.SetInAndOutFormats(out_format, out_format)
        obmol = ob.OBMol()
        mols_opt = []
        mols_in_file = obconversion.ReadFile(obmol, fnameout_noopt)
        while mols_in_file:
            mols_opt.append(obmol)
            obmol = ob.OBMol()
            mols_in_file = obconversion.Read(obmol)
        obconversion.CloseOutFile()

        # WRITE optimized conformers
        optimized_strings_mols = []
        iconf = 0
        m = "\t\t# Conformation E_initial E_final (in kcal/mol)\n"
        m += "\t\t# ============================================\n"
        for iobmol in mols_opt:
            ff.Setup(iobmol)
            e1 = ff.Energy()
            if minimize_iterations > 0:
                ff.ConjugateGradients(minimize_iterations)
            ff.GetCoordinates(iobmol)
            e2 = ff.Energy()
            optimized_strings_mols.append(obconversion.WriteString(iobmol))
            m = "\t\t{0:>6d}  {1:>10.3f}  {2:>10.3f}".format(iconf, e1, e2)
            # Deprecated pandas
            # self._df_conformers = self._df_conformers.append({'Conformations': iconf,
            #                                                   'IniEnergy': e1,
            #                                                   'OptEnergy': e2,
            #                                                   'Cluster': iconf}, ignore_index=True)
            if self._df_conformers is None:
                self._df_conformers = pd.DataFrame({'Conformations': [iconf],
                                                    'IniEnergy': [e1],
                                                    'OptEnergy': [e2],
                                                    'Cluster': [iconf]})
            else:
                new_row = pd.DataFrame({'Conformations': [iconf],
                                        'IniEnergy': [e1],
                                        'OptEnergy': [e2],
                                        'Cluster': [iconf]})
                self._df_conformers = pd.concat([self._df_conformers, new_row], ignore_index=True)

            iconf += 1

        self._df_conformers = self._df_conformers.sort_values('OptEnergy')

        if debug_flag:
            # WRITE conformers with optimization
            fnameout_opt = "{}_allconf_ob_optMM_trj.{}".format(pattern, out_format)
            with open(fnameout_opt, 'w') as f:
                for output_string in optimized_strings_mols:
                    f.write(output_string)

        return optimized_strings_mols

    # =========================================================================
    def _cluster_conformers(self, optimized_strings_mols, pattern="Mol", cutoff_rmsddock=2.0,
                            energy_threshold=99999,
                            maximum_number_clusters=100, ff_name=None):

        # Get the index of the lowest-energy structure
        idx_ref_str = int(self._df_conformers[['OptEnergy']].idxmin()[-1])
        ref_name_mol2 = "ref_tmp.mol2"
        with open(ref_name_mol2, 'w') as f:
            f.write(optimized_strings_mols[idx_ref_str])

        # Loop over the rest of conformations and calculate rmsd
        for iconf in range(len(optimized_strings_mols)):
            tmp_name_mol2 = "target_tmp.mol2"
            with open(tmp_name_mol2, 'w') as f:
                f.write(optimized_strings_mols[iconf])

            if self._exec_rmsddock is not None:
                cmd = '{} {} {} -s'.format(self._exec_rmsddock, ref_name_mol2, tmp_name_mol2)
                rmsd_dock = subprocess.check_output(cmd, shell=True)
                self._df_conformers.at[iconf, 'RMSD'] = rmsd_dock.decode()

        # Building the clusters
        cluster = defaultdict(dict)
        icluster = 0

        for index in self._df_conformers.index:

            energy = self._df_conformers['OptEnergy'][index]
            iconf = self._df_conformers['Conformations'][index]

            if icluster == 0:
                threshold = energy + energy_threshold
                icluster += 1
                cluster[icluster] = {"seed": index, "lowest_energy": energy, "highest_energy": energy,
                                     "nelements": 0, "elements": [], "pairs": [], "files": []}
                cluster[icluster]["pairs"].append([0.000, energy])
                cluster[icluster]["files"].append("")
                cluster[icluster]["nelements"] += 1
                cluster[icluster]["elements"].append(index)
                min_energy = energy
                self._df_conformers.at[index, 'Cluster'] = icluster

            elif energy < threshold:
                tmp_name_mol2 = "target_tmp.mol2"
                with open(tmp_name_mol2, 'w') as f:
                    f.write(optimized_strings_mols[index])
                found = False
                for i in range(1, icluster + 1):
                    idx = cluster[i]["seed"]
                    ref_name_mol2 = "ref_tmp.mol2"
                    with open(ref_name_mol2, 'w') as f:
                        f.write(optimized_strings_mols[int(idx)])

                    cmd = '{} {} {} -h -s'.format(self._exec_rmsddock, ref_name_mol2, tmp_name_mol2)
                    rmsd_dock = float(subprocess.check_output(cmd, shell=True).decode())

                    if float(rmsd_dock) < cutoff_rmsddock:
                        cluster[i]["pairs"].append([rmsd_dock, energy])
                        cluster[i]["files"].append("")
                        cluster[i]["highest_energy"] = energy
                        cluster[i]["nelements"] += 1
                        cluster[i]["elements"].append(index)
                        self._df_conformers.at[index, 'Cluster'] = i
                        found = True
                        break

                if not found:
                    if icluster < maximum_number_clusters:
                        icluster += 1
                        cluster[icluster] = {"seed": index, "lowest_energy": energy, "highest_energy": energy,
                                             "nelements": 0, "elements": [], "pairs": [], "files": []}
                        cluster[icluster]["files"].append(index)
                        cluster[icluster]["nelements"] += 1
                        cluster[icluster]["elements"].append(index)
                        cluster[icluster]["pairs"].append([rmsd_dock, energy])
                        self._df_conformers.at[index, 'Cluster'] = icluster
                    else:
                        if cluster[maximum_number_clusters+1]:
                            cluster[maximum_number_clusters+1]["pairs"].append([rmsd_dock, energy])
                            cluster[maximum_number_clusters+1]["files"].append(index)
                            cluster[maximum_number_clusters+1]["highest_energy"] = energy
                            cluster[maximum_number_clusters+1]["nelements"] += 1
                            cluster[maximum_number_clusters+1]["elements"].append(index)
                            self._df_conformers.at[index, 'Cluster'] = icluster
                        else:
                            cluster[maximum_number_clusters+1] = {"seed": iconf, "lowest_energy": energy, "highest_energy": energy,
                                                                  "nelements": 0, "pairs": [], "files": []}
                            cluster[maximum_number_clusters+1]["files"].append(index)
                            cluster[maximum_number_clusters+1]["nelements"] += 1
                            cluster[maximum_number_clusters+1]["elements"].append(index)
                            cluster[maximum_number_clusters+1]["pairs"].append([rmsd_dock, energy])
                            self._df_conformers.at[index, 'Cluster'] = icluster
            else:
                pass

        os.remove("ref_tmp.mol2")
        os.remove("target_tmp.mol2")

        return cluster

    # =========================================================================
    def _write_gaussian(self, localdir, optimized_strings_mols, cluster, pattern="QM", g16_key="#p 6-31g* mp2",
                        g16_mem=4000, g16_nproc=4, g16_extra_info='', charge=0, multiplicity=1):

        # Create directory for gaussian inputs
        parent_dir = os.path.join(localdir, "{}".format("{}_g16_results/").format(pattern))
        if not os.path.isdir(parent_dir):
            os.mkdir(parent_dir)

        for key, item in cluster.items():
            iconf_seed = int(item['seed'])
            energy = item['lowest_energy']

            obconversion = ob.OBConversion()
            obconversion.SetInAndOutFormats('mol2', 'xyz')
            obmol = ob.OBMol()
            xyz_obmol = obconversion.ReadString(obmol, optimized_strings_mols[iconf_seed])
            if not xyz_obmol:
                m = "Something is wrong with Gaussian16 input files 1\n"
                m += "Gaussian input files are not written"
                print(m) if self._logger is None else self._logger.error(m)
                return None

            xyz_list = obconversion.WriteString(obmol).split("\n")

            fname = os.path.join(parent_dir, "{0:s}_{1:03d}_gaussian.com".format(pattern, iconf_seed))
            fchk_name = "{0:s}_{1:03d}_gaussian.chk".format(pattern, iconf_seed)
            with open(fname, 'w') as f:
                f.writelines("%chk={}\n".format(fchk_name))
                f.writelines("%nproc={}\n".format(g16_nproc))
                f.writelines("%mem={}Mb\n".format(g16_mem))
                f.writelines("{}\n".format(g16_key))
                f.writelines("\nConformer number {0:03d}. Energy MM = {1:.3f} kcal/mol\n".format(iconf_seed, energy))
                f.writelines("\n")
                f.writelines("{0:1d} {1:1d}\n".format(charge, multiplicity))
                for line in xyz_list[2:]:
                    f.writelines(line+"\n")

                for item2 in g16_extra_info:
                    if not item2.strip() == '':
                        f.writelines(item2)
                        f.writelines("\n")
                f.writelines("\n")

    # =========================================================================
    @staticmethod
    def write_all_conformers_to_mol2(filename, optimized_string_mols):

        with open(filename, 'w') as f:
            for imol in optimized_string_mols:
                f.write(imol)

