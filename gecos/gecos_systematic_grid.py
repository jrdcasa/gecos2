from gecos.gecos_rdkit import GecosRdkit
import rdkit.Chem.Draw
import rdkit.Chem.rdMolTransforms
import rdkit.Chem.AllChem
from rdkit import RDLogger
import itertools
import datetime
import os
from collections import defaultdict
import warnings
warnings.simplefilter("ignore", Warning)


class GecosSystematicGrid(GecosRdkit):

    # #######################################################################
    def __init__(self, filename, pattern="QM",
                 total_charge=0, bond_perception=False, logger=None):

        GecosRdkit.__init__(self, filename, total_charge, bond_perception, logger)
        self._pattern = pattern
        self._ndihedrals = 0
        self._dih_list = []     # Indices starts at zero (rdkit)
        self._dih_delta_deg = []
        self._dih_only_ttt = []  # If 1 only trans, gauche+ and gauche- states are generated

    # =========================================================================
    def _cluster_fake_conformers(self, rms_threshold=1.0, energy_threshold=99999.0,
                                 maximum_number_clusters=100):

        """
        Clustering conformers.

        The value of `dmat` is the lower half matrix of RMSD. The conformers will be
        alligned to the first conformer (reference). The rmsd is calculated using all atoms.
        The algorithm used in the clustering is the Ramos algorithm (not published)

        Args:
            conf_prop_dict (dict):
            rms_threshold (float): Threshold for the clustering algorithm in angstroms.
            energy_threshold (float):

        Returns:
            A tuple of tuples containing information about the cluster.
            Example -> rms_cluster = ( (0,1,2), (3) ) Two clusters, the first one contains the conformers 0, 1, 2.

        """

        # Get heavy atoms indexes
        list_indices_heavy_atoms = []
        for iatom in self.mol_rdkit.GetAtoms():
            if iatom.GetAtomicNum() != 1:
                list_indices_heavy_atoms.append(iatom.GetIdx())

        # Allign the conformers and returns the RMS matrix of the conformers of a molecule.
        # As a side-effect, the
        # conformers will be aligned to the first conformer (i.e. the reference) and will
        # left in the aligned state
        # dmat for 5 conformers [0, 1, 2, 3, 4]: (nrms= 4+3+2+1 = 10)
        #   [a,               [0-1,
        #   b, c               0-2, 1-2
        #   d, e, f            0-3, 1-3, 2-3
        #   g, h, i, j]        0-4, 1-4, 2-4, 3-4]
        # This matrix can be used directly in the Butina clustering algorithm
        # Recover the index for a given pair:
        #   Pair(label1,label2) --> Take the max (label1,label2) = max
        #   index = (sum(max-1 to 1)) + label2
        #   Example:
        #   Pair(2,3) --> index = (2 + 1) + 2 = 5  --> OK
        #   Pair(1,0) --> index = (0) + 0     = 0  --> OK
        #   Pair(3,4) --> index = (3+2+1) + 3 = 9  --> OK
        # dmat = rdkit.Chem.AllChem.GetConformerRMSMatrix(self.mol_rdkit,
        #                                                 atomIds=list_indices_heavy_atoms,
        #                                                 prealigned=False)

        n_conformers = self.mol_rdkit.GetNumConformers()

        # My clustering method (Ramos algorithm). The conformers are presorted
        # by optimized MM energy in the self._df_conformers dataframe
        cluster = defaultdict(dict)
        icluster = 0

        for index in self._df_conformers.index:

            energy = self._df_conformers['OptEnergy'][index]
            iconf = self._df_conformers['Conformations'][index]

            # print("{0:d} {1:d} of {2:d} (# clusters {3:d})".format(int(index), int(iconf), n_conformers, icluster))

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
                found = False
                iconf_ref = cluster[1]["seed"]
                iconf_target = index
                rmsd_noh = self.getconformerrms(iconf_ref, iconf_target,
                                                atomids=list_indices_heavy_atoms, align=True)
                if not found:
                    if icluster < maximum_number_clusters:
                        icluster += 1
                        cluster[icluster] = {"seed": index, "lowest_energy": energy, "highest_energy": energy,
                                             "nelements": 0, "elements": [], "pairs": [], "files": []}
                        cluster[icluster]["files"].append(index)
                        cluster[icluster]["nelements"] += 1
                        cluster[icluster]["elements"].append(index)
                        cluster[icluster]["pairs"].append([rmsd_noh, energy])
                        self._df_conformers.at[index, 'Cluster'] = icluster
                    else:
                        if cluster[maximum_number_clusters+1]:
                            cluster[maximum_number_clusters+1]["pairs"].append([rmsd_noh, energy])
                            cluster[maximum_number_clusters+1]["files"].append(index)
                            cluster[maximum_number_clusters+1]["highest_energy"] = energy
                            cluster[maximum_number_clusters+1]["nelements"] += 1
                            cluster[maximum_number_clusters+1]["elements"].append(index)
                            self._df_conformers.at[index, 'Cluster'] = icluster
                        else:
                            cluster[maximum_number_clusters+1] = {"seed": index, "lowest_energy": energy,
                                                                  "highest_energy": energy,
                                                                  "elements": [], "nelements": 0,
                                                                  "pairs": [], "files": []}
                            cluster[maximum_number_clusters+1]["files"].append(index)
                            cluster[maximum_number_clusters+1]["nelements"] += 1
                            cluster[maximum_number_clusters+1]["elements"].append(index)
                            cluster[maximum_number_clusters+1]["pairs"].append([rmsd_noh, energy])
                            self._df_conformers.at[index, 'Cluster'] = icluster
            else:
                pass

        return cluster

    # #######################################################################
    def generate_systematic_conformers(self, localdir, optimize,
                                       maxmmoptiters=1000, ndihedrals=None,
                                       dih_list=None, cluster_method="RMSD", cluster_threshold=0,
                                       g16_nproc=4, g16_mem=4000, g16_key="#p 6-31g* mp2", g16_extra_info="",
                                       charge=0, multiplicity=1, write_gaussian=True, pattern="conformers"):

        # Disable warning messages from RdKit
        RDLogger.DisableLog('rdApp.*')

        # Correct indices to start at zero
        self._ndihedrals = ndihedrals
        for item in dih_list:
            self._dih_list.append([item[0]-1, item[1]-1, item[2]-1, item[3]-1])
            self._dih_delta_deg.append(item[4])
            self._dih_only_ttt.append(item[5])

        # Calculate the number of steps for each angle
        angle_steps = defaultdict(list)
        for idx_idih in range(0, self._ndihedrals):
            if self._dih_only_ttt[idx_idih] == 0:
                nsteps = int(360./self._dih_delta_deg[idx_idih])
                iangle = 0
                for istep in range(0, nsteps):
                    angle_steps[idx_idih].append(iangle)     # degrees
                    iangle += self._dih_delta_deg[idx_idih]
            else:
                nsteps = 3
                angle_steps[idx_idih] = [60, 180, 300]

        # Generate the conformers
        all_lists = []
        for idih, list_idih in angle_steps.items():
            all_lists.append(list_idih)
        allcombinations = list(itertools.product(*all_lists))
        nconfs = len(allcombinations)

        m = "\t\t**************** GENERATE CONFORMERS ***************\n"
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m += "\t\t 1. Generating {} conformers by systematic approach ({})".format(nconfs, now)
        print(m) if self._logger is None else self._logger.info(m)

        # Set dihedral angles and add conformer to rdkit
        confid = 0
        energy = []
        for icomb_angle in allcombinations:
            conf = self.mol_rdkit.GetConformer(confid)
            confid += 1
            for idih in range(len(icomb_angle)):
                at1 = self._dih_list[idih][0]
                at2 = self._dih_list[idih][1]
                at3 = self._dih_list[idih][2]
                at4 = self._dih_list[idih][3]
                rdkit.Chem.rdMolTransforms.SetDihedralDeg(conf, at1, at2, at3, at4, icomb_angle[idih])
            self.mol_rdkit.AddConformer(conf, assignId=True)

        # Add dihedral constraint
        molcopyprop = rdkit.Chem.AllChem.MMFFGetMoleculeProperties(self.mol_rdkit)
        firstopt = True
        conformerpropsdict = defaultdict()
        # The first conformer is not taken into account
        for iconf in range(1, self.mol_rdkit.GetNumConformers()):
            ff = rdkit.Chem.AllChem.MMFFGetMoleculeForceField(self.mol_rdkit, molcopyprop, confId=iconf)
            conf = self.mol_rdkit.GetConformer(id=iconf)
            for idih in range(self._ndihedrals):
                at1 = self._dih_list[idih][0]
                at2 = self._dih_list[idih][1]
                at3 = self._dih_list[idih][2]
                at4 = self._dih_list[idih][3]
                angle = rdkit.Chem.rdMolTransforms.GetDihedralDeg(conf, at1, at2, at3, at4)
                ff.MMFFAddTorsionConstraint(at1, at2, at3, at4, False, angle - .1, angle + .1, 10000.0)
            # Constrained optimization using MMFF4
            if iconf % int(50.) == 0:
                m = "\t\t\t MM Optimization {} of {}".format(iconf, nconfs)
                print(m) if self._logger is None else self._logger.info(m)
            if optimize:
                now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
                if firstopt:
                    m = "\t\t 2. Optimizing {} conformers using MMFF. ({})".format(nconfs, now)
                    print(m) if self._logger is None else self._logger.info(m)
                    firstopt = False

                ff.Minimize(maxIts=maxmmoptiters)
            else:
                now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
                if firstopt:
                    m = "\t\t 2. Conformers have not been optimized. " \
                         "Strong atom overlapping is expected.({})".format(nconfs, now)
                    print(m) if self._logger is None else self._logger.info(m)
                    firstopt = False

        # The first conformer is not taken into account
        for conf_id in range(1, self.mol_rdkit.GetNumConformers()):
            # Calculate energy and minimize
            props = self.mm_calc_energy(conf_id, ff_name="MMFF", minimize_iterations=0)
            conformerpropsdict[conf_id-1] = props
        self._df_conformers = self._df_conformers.sort_values('OptEnergy')

        # cluster the MM conformers
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\t\t 3. Cluster Conformers ({}) using {} based method".format(now, cluster_method)
        print(m) if self._logger is None else self._logger.info(m)
        mm_clusters = self._cluster_fake_conformers(cluster_threshold, maximum_number_clusters=nconfs)

        pattern_single = os.path.splitext(os.path.split(self._filename)[-1])[0]
        fpdboutput = os.path.join(localdir, "{}_allconf_systgrid.pdb".format(pattern_single))
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\t\t 4. Write Conformers to PDB ({})\n".format(now)
        m += "\t\t\tFilename: {}".format(fpdboutput)
        print(m) if self._logger is None else self._logger.info(m)
        self.write_all_conformers_to_pdb(fpdboutput)

        if write_gaussian:
            now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
            m = "\t\t 5. Write Conformers to Gaussian16 inputs ({})\n".format(now)
            m += "\t\t\tInputfiles in {}_g16_results/ directory".format(pattern)
            print(m) if self._logger is None else self._logger.info(m)
            list_conf_id = [i for i in range(0, nconfs)]
            self.write_min_cluster_conformers_to_g16(localdir, pattern, conformerpropsdict, list_conf_id,
                                                     g16_nproc, g16_mem, g16_key, g16_extra_info, charge, multiplicity)

        m = "\t\t**************** GENERATE CONFORMERS ***************\n"
        print(m) if self._logger is None else self._logger.info(m)

        m = "\t\t**************** CONFORMER MM ENERGIES ***************"
        m += "\n\t\tEnergies in kcal/mol, RMSD in angstroms."
        m += "\n\t\tRMSD threshold = {0:5.3f} angstroms.".format(cluster_threshold)
        m += "\n\t\tNumber of clusters = {0:4d}".format(len(mm_clusters))
        m += "\n\t\tNumber of dihedrals to systematically scan = {}".format(self._ndihedrals)
        for idih in range(self._ndihedrals):
            at1 = self._dih_list[idih][0]+1
            at2 = self._dih_list[idih][1]+1
            at3 = self._dih_list[idih][2]+1
            at4 = self._dih_list[idih][3]+1
            m += "\n\t\t\tDihedral {}: {}, {}, {}, {}".format(idih, at1, at2, at3, at4)
        m += "\n\t\t\t\t(Indices start at 1)"
        print(m) if self._logger is None else self._logger.info(m)
        m = "\t\tCluster Conformer_ID  {}_energy  Relative_energy " \
            "RMSD   dih_values\n".format("MMFF")
        separator = "\t\t"+len(m)*"-"+"\n"
        m += separator
        min_energy = mm_clusters[1]['lowest_energy']
        for icluster, ival in mm_clusters.items():
            angles_str = ""
            iconf = ival['seed']
            conf = self.mol_rdkit.GetConformer(id=iconf)
            for idih in range(self._ndihedrals):
                at1 = self._dih_list[idih][0]
                at2 = self._dih_list[idih][1]
                at3 = self._dih_list[idih][2]
                at4 = self._dih_list[idih][3]
                angle = rdkit.Chem.rdMolTransforms.GetDihedralDeg(conf, at1, at2, at3, at4)
                angles_str += "{0:.1f}, ".format(angle)
            angles_str = angles_str[:-2]
            energy_abs = ival['lowest_energy']
            m += "\t\t {0:^5d}  {1:^12d}  {2:^12.2f}  {3:^12.2f}  {4:^10.3f} {5:s} \n".\
                 format(icluster, iconf, energy_abs,
                        energy_abs-min_energy,
                        mm_clusters[icluster]['pairs'][0][0], angles_str)
        m += separator
        m += "\t\t**************** CONFORMER MM ENERGIES ***************\n"
        print(m) if self._logger is None else self._logger.info(m)

        return True
