import numpy as np
import os
import shutil
from collections import defaultdict
from MDAnalysis import transformations as trans
import MDAnalysis as mda
import utils
from openbabel import openbabel as ob
import datetime
import warnings
warnings.simplefilter("ignore", Warning)


class GecosExtractNeighbors(object):

    __slots__ = ["_com_list", "_filename", "_universe", "_nmols_array", "_nmols", "_logger",
                 "_com_filename_pdb", "_neighbor_sphere_com", "_radius", "_isheader_print",
                 "_mol_pybabel", "_pattern", "_nmonomers", "_npairs", "_ntrimers",
                 "_extract_monomer", "_extract_pair", "_extract_trimer", "_imol1_ref", "_imol2_ref", "_ipair_ref"]

    # #######################################################################
    def __init__(self, filename, write_tcl=True,
                 radius=7.0, pattern="QM", logger=None):

        self._com_list = []
        self._filename = filename
        self._universe = mda.Universe(filename)
        self._nmols = None
        self._radius = radius
        self._isheader_print = False
        self._mol_pybabel = None
        self._pattern = pattern
        self._nmonomers = 0
        self._npairs = 0
        self._ntrimers = 0
        self._extract_monomer = True
        self._extract_pair = False
        self._extract_trimer = False
        self._imol1_ref = None
        self._imol2_ref = None
        self._ipair_ref = None

        # Openbabel object to allign and clusterize the conformers
        file_format = filename.strip().split('.')[-1]
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

        self._logger = logger
        self._com_filename_pdb = None
        self._neighbor_sphere_com = defaultdict(list)

        if not self._isheader_print:
            utils.print_header(self._logger)
            self._isheader_print = True

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\t\t******** CREATION OF THE GECOS EXTRACT NEIGHBORS OBJECT ********\n\n"
        m += "\t\tConformers will be extracted from a frame of a MD simulation.\n"
        m += "\t\t Find number of molecules in the frame.({})\n".format(now)
        self._findmolecules()
        m += "\n"
        m += "\t\tFormula             = {}\n".format(self._mol_pybabel.GetFormula())
        m += "\t\tNumber of atoms     = {}\n".format(self._mol_pybabel.NumAtoms())
        m += "\t\tMolecular weight    = {0:.2f} g/mol\n".format(self._mol_pybabel.GetMolWt())
        m += "\t\tNumber of rings     = {}\n".format(len(self._mol_pybabel.GetSSSR()))
        m += "\t\tNumber of molecules = {}\n".format(self._nmols)
        m += "\t\t******** CREATION OF THE GECOS EXTRACT NEIGHBORS OBJECT ********\n\n"
        print(m) if self._logger is None else self._logger.info(m)

    # =========================================================================
    def extract_conformers(self, localdir, calc_com=True, method="sphere", extract_type="monomer",
                           write_gaussian=True, g16_key="#p 6-31g* mp2",
                           g16_mem=4000, g16_extra_info="",
                           g16_nproc=4, pattern="conformers", charge=0, multiplicity=1,
                           isclustering=True, debug_flag=False, onlydifferentmols=False):
        """

        :param localdir:
        :param calc_com:
        :param method:
        :param extract_type:
        :param write_gaussian:
        :param g16_key:
        :param g16_mem:
        :param g16_nproc:
        :param pattern:
        :param charge:
        :param multiplicity:
        :param isclustering:
        :param debug_flag:
        :param onlydifferentmols
        :param g16_extra_info (str):

        :return:
        """

        if extract_type.upper() == "MONOMER":
            self._extract_monomer = True
            self._extract_pair = False
            self._extract_trimer = False
        elif extract_type.upper() == "PAIR":
            self._extract_monomer = False
            self._extract_pair = True
            self._extract_trimer = False
        elif extract_type.upper() == "TRIMER":
            self._extract_monomer = False
            self._extract_pair = False
            self._extract_trimer = True

        m = "\t\t**************** EXTRACT CONFORMERS ***************\n"
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m += "\t\t 1. Unwrapping molecules to calculate center of masses ({})".format(now)
        print(m) if self._logger is None else self._logger.info(m)
        # Unwrap trajectory:
        ag = self._universe.atoms
        transform = trans.unwrap(ag)
        self._universe.trajectory.add_transformations(transform)

        # COM must be done with the coordinates unwrapped.
        # Both the self._com_list and the coordinates of the universe
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\t\t 2. Calculating center of masses ({})".format(now)
        print(m) if self._logger is None else self._logger.info(m)
        if calc_com:
            self._com()

        # FOR DEBUG ================
        # ag = self._universe.atoms
        # ag.write("com_moved.pdb")
        # FOR DEBUG ================

        # Write a tmp pdb and redefine the universe in the the tajectory.
        # This is neccesary because MDAnalysis does not allow to make
        # two transformation to a given universe.
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\t\t 3. Wrapping atoms to visualize (wrapper_atoms.pdb) ({})".format(now)
        print(m) if self._logger is None else self._logger.info(m)

        ag.write("wrapped_atoms.pdb")
        self._universe = mda.Universe("wrapped_atoms.pdb")
        ag = self._universe.atoms
        transform = trans.wrap(ag)
        ag.write("wrapped_atoms.pdb")

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\t\t 4. Extracting conformations from frame ({})".format(now)
        print(m) if self._logger is None else self._logger.info(m)
        if method.upper() == "SPHERE_COM":
            self._find_sphere_com(radius=self._radius, write_gaussian=write_gaussian, g16_key=g16_key,
                                  g16_mem=g16_mem, g16_nproc=g16_nproc, pattern=pattern,
                                  g16_extra_info=g16_extra_info,
                                  charge=charge, multiplicity=multiplicity, onlydifferentmols=onlydifferentmols)
        elif method.upper() == "SPHERE_MINATOM":
            self._find_sphere_minatom(radius=self._radius, write_gaussian=write_gaussian, g16_key=g16_key,
                                      g16_mem=g16_mem, g16_nproc=g16_nproc, pattern=pattern,
                                      g16_extra_info=g16_extra_info,
                                      charge=charge, multiplicity=multiplicity, onlydifferentmols=onlydifferentmols)
        else:
            m = "\t\t Method to extract {} is not allowed.".format(method.upper())
            print(m) if self._logger is None else self._logger.error(m)
            exit()

        m = "\t\t**************** EXTRACT CONFORMERS ***************\n"
        print(m) if self._logger is None else self._logger.info(m)

    # #######################################################################
    def _findmolecules(self):

        self._nmols_array = defaultdict(list)
        for idx, frg in enumerate(self._universe.atoms.fragments):
            self._nmols_array[idx] = list(frg.indices)
        self._nmols = len(self._universe.atoms.fragments)

        return None

    # #######################################################################
    def _com(self):

        for ikey, imol in self._nmols_array.items():
            seltxt = "index "
            for iat in imol:
                seltxt += "{0:d} ".format(iat)
            g = self._universe.select_atoms(seltxt)
            self._com_list.append(g.center_of_mass())

        # Wrap COM
        tmp_com = []
        box_x = self._universe.dimensions[0]
        box_y = self._universe.dimensions[1]
        box_z = self._universe.dimensions[2]
        for imol, icom in enumerate(self._com_list):
            imgx = np.floor(icom[0]/box_x)*box_x
            imgy = np.floor(icom[1]/box_y)*box_y
            imgz = np.floor(icom[2]/box_z)*box_z
            x = icom[0] - imgx
            y = icom[1] - imgy
            z = icom[2] - imgz
            # It the COM is moved, we need also to move the molecule coordinates.
            image = [imgx, imgy, imgz]
            if any(image) != 0.0:
                for iatom in self._nmols_array[imol]:
                    self._universe.coord[iatom][0] -= imgx
                    self._universe.coord[iatom][1] -= imgy
                    self._universe.coord[iatom][2] -= imgz
            tmp_com.append([x, y, z])
        self._com_list = tmp_com

    # #######################################################################
    def _write_com_pdb(self, filename_com_pdb="coords_com.pdb", wrap=False):

        """
        Write a pdb file to check the center of mass.
        Adapted from MDAnalysis software (https://www.mdanalysis.org/)
        """

        self._com_filename_pdb = filename_com_pdb

        fmt = {
            'ATOM': (
                "ATOM  {serial:5d} {name:<4s}{altLoc:<1s}{resName:<4s}"
                "{chainID:1s}{resSeq:4d}{iCode:1s}"
                "   {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}{occupancy:6.2f}"
                "{tempFactor:6.2f}      {segID:<4s}{element:>2s}\n"),
            'HETATM': (
                "HETATM{serial:5d} {name:<4s}{altLoc:<1s}{resName:<4s}"
                "{chainID:1s}{resSeq:4d}{iCode:1s}"
                "   {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}{occupancy:6.2f}"
                "{tempFactor:6.2f}      {segID:<4s}{element:>2s}\n"),
            'REMARK': "REMARK     {0}\n",
            'COMPND': "COMPND    {0}\n",
            'HEADER': "HEADER    {0}\n",
            'TITLE': "TITLE     {0}\n",
            'MODEL': "MODEL     {0:>4d}\n",
            'NUMMDL': "NUMMDL    {0:5d}\n",
            'ENDMDL': "ENDMDL\n",
            'END': "END\n",
            'CRYST1': ("CRYST1{0:9.3f}{1:9.3f}{2:9.3f}"
                       "{3:7.2f}{4:7.2f}{5:7.2f} "
                       "{6:<11s}{7:4d}\n"),
            'CONECT': "CONECT{0}\n"
        }

        a = self._universe.dimensions[0]
        b = self._universe.dimensions[1]
        c = self._universe.dimensions[2]
        alpha = self._universe.dimensions[3]
        beta = self._universe.dimensions[4]
        gamma = self._universe.dimensions[5]

        atom_kind_molecule_label = None
        spacegroup = "P -1"
        zvalue = 1
        radtodeg = 180/np.pi
        with open(filename_com_pdb, 'w') as fpdb:

            fpdb.write(fmt['REMARK'].format('Created with polyanogro(J.Ramos)'))
            fpdb.write(fmt['CRYST1'].format(a, b, c, alpha, beta, gamma,
                                            spacegroup, zvalue))

            for idx in range(self._nmols):

                resname = "{0:03d}".format(idx)

                if atom_kind_molecule_label is None:
                    fpdb.write(fmt['HETATM'].format(
                        serial=idx+1,
                        name='C',
                        altLoc=" ",
                        resName=resname,
                        chainID=" ",
                        resSeq=idx,
                        iCode=" ",
                        pos=[i for i in self._com_list[idx]],
                        occupancy=1.0,
                        tempFactor=1.0,
                        segID="    ",
                        element='C'
                    ))

            fpdb.write('END\n')

        self._write_tcl_atoms_com()

    # # #######################################################################
    def _find_sphere_com(self, radius=None,
                         write_gaussian=True, g16_key="#p 6-31g* mp2",
                         g16_mem=4000, g16_nproc=4, pattern="conformers",
                         g16_extra_info="",
                         charge=0, multiplicity=1, onlydifferentmols=False):

        """
        Find other com's inside a sphere of this radius
        Args:
            radius: in angstroms
            extract_monomers
            extract_pairs
        Returns:

        """

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\t\t    4.a Using method SPHERE_COM with radius {} angstroms ({})".format(self._radius, now)
        print(m) if self._logger is None else self._logger.info(m)

        if self._extract_monomer:
            now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
            m = "\t\t    4.b Extracting monomers ({})".format(now)
            print(m) if self._logger is None else self._logger.info(m)
            self._extract_monomers(write_gaussian=write_gaussian, g16_key=g16_key,
                                   g16_mem=g16_mem, g16_nproc=g16_nproc, pattern=pattern,
                                   g16_extra_info=g16_extra_info,
                                   charge=charge, multiplicity=multiplicity)

        elif self._extract_pair:
            now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
            m = "\t\t    4.b Extracting pairs ({})".format(now)
            if onlydifferentmols:
                m += "\n\t\t      Only pairs containing different molecules are extracted from the trajectory."
            print(m) if self._logger is None else self._logger.info(m)
            self._extract_pairs_sphere_com(self._radius, write_gaussian=write_gaussian, g16_key=g16_key,
                                           g16_mem=g16_mem, g16_nproc=g16_nproc, pattern=pattern,
                                           charge=charge, multiplicity=multiplicity,
                                           g16_extra_info=g16_extra_info,
                                           onlydifferentmols=onlydifferentmols)
            m = "\t\t        Number of pairs created = {}".format(self._npairs)
            print(m) if self._logger is None else self._logger.info(m)
        elif self._extract_trimer:

            print("TO BE IMPLEMENTED: TRIMER")
            exit()

        return None

    # # #######################################################################
    def _find_sphere_minatom(self, radius=None,
                             write_gaussian=True, g16_key="#p 6-31g* mp2",
                             g16_mem=4000, g16_nproc=4, pattern="conformers",
                             g16_extra_info="",
                             charge=0, multiplicity=1, onlydifferentmols=False):

        """
        Find molecules in wihch the minimum dustance between atoms is inside a sphere of this radius
        Args:
            radius: in angstroms
            extract_monomers
            extract_pairs
        Returns:

        """

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\t\t    4.a Using method SPHERE_MINATOM with radius {} angstroms ({})".format(self._radius, now)
        print(m) if self._logger is None else self._logger.info(m)

        if self._extract_monomer:
            now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
            m = "\t\t    4.b Extracting monomers ({})".format(now)
            print(m) if self._logger is None else self._logger.info(m)
            self._extract_monomers(write_gaussian=write_gaussian, g16_key=g16_key,
                                   g16_mem=g16_mem, g16_nproc=g16_nproc, pattern=pattern,
                                   charge=charge, multiplicity=multiplicity)

        elif self._extract_pair:
            now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
            m = "\t\t    4.b Extracting pairs ({})".format(now)
            if onlydifferentmols:
                m += "\n\t\t      Only pairs containing different molecules are extracted from the trajectory."
            print(m) if self._logger is None else self._logger.info(m)
            self._extract_pairs_sphere_minatom(self._radius, write_gaussian=write_gaussian, g16_key=g16_key,
                                               g16_mem=g16_mem, g16_nproc=g16_nproc, pattern=pattern,
                                               charge=charge, multiplicity=multiplicity,
                                               onlydifferentmols=onlydifferentmols)
            m = "\t\t        Number of pairs created = {}".format(self._npairs)
            print(m) if self._logger is None else self._logger.info(m)
        elif self._extract_trimer:

            print("TO BE IMPLEMENTED: TRIMER")
            exit()

        return None

    # #######################################################################
    def _extract_monomers(self, write_gaussian=True, g16_key="#p 6-31g* mp2",
                          g16_mem=4000, g16_nproc=4, pattern="conformers", g16_extra_info="",
                          charge=0, multiplicity=1):

        # Write Monomoners
        localdir = os.getcwd()
        path_monomers = os.path.join(localdir, self._pattern+"_g16_monomers")
        if os.path.isdir(path_monomers):
            shutil.rmtree(path_monomers)
        os.mkdir(path_monomers)

        # Write temporal pdbs containing monomers
        tmp_ob_nmols = []
        obconversion = ob.OBConversion()
        informat = 'pdb'
        outformat = 'mol2'
        obconversion.SetInAndOutFormats(informat, outformat)
        for imonomer in self._nmols_array:
            fpdbname = self._write_pdb_cluster(path_monomers, [imonomer], [0, 0, 0])
            imol = ob.OBMol()
            obconversion.ReadFile(imol, fpdbname)
            tmp_ob_nmols.append(imol)
            self._nmonomers += 1

        # WRITE extracted monomers ===================================================
        obconversion = ob.OBConversion()
        outformat = 'mol2'
        obconversion.SetOutFormat(outformat)
        optimized_strings_mols = []
        formula_string_mols = []
        for iconf in range(self._nmonomers):
            optimized_strings_mols.append(obconversion.WriteString(tmp_ob_nmols[iconf]))
            formula_string_mols.append(tmp_ob_nmols[iconf].GetFormula())
        # for iconf in range(self._nmonomers):
        #     optimized_strings_mols.append(obconversion.WriteString(tmp_ob_nmols[iconf]))

        # Write Gaussian
        if write_gaussian:
            self._write_gaussian(localdir, optimized_strings_mols, formula_string_mols, "conformers",
                                 pattern=pattern, g16_key=g16_key, g16_extra_info=g16_extra_info,
                                 g16_mem=g16_mem, g16_nproc=g16_nproc,
                                 charge=charge, multiplicity=multiplicity)
        pass

    # #######################################################################
    def _extract_pairs_sphere_com(self, radius, write_gaussian=True, g16_key="#p 6-31g* mp2",
                                  g16_mem=4000, g16_nproc=4, pattern="conformers",
                                  charge=0, multiplicity=1, onlydifferentmols=False):

        from ext_libc.c_distances_openmp import calc_minimum_image_neighbors_com

        # Write Pairs
        localdir = os.getcwd()
        path_pairs = os.path.join(localdir, self._pattern+"_g16_pairs")
        if os.path.isdir(path_pairs):
            shutil.rmtree(path_pairs)
        os.mkdir(path_pairs)

        # Box
        a = self._universe.dimensions[0]
        b = self._universe.dimensions[1]
        c = self._universe.dimensions[2]
        box = np.array([a, b, c], dtype=np.float32)

        # Check radius
        if radius is None:
            radius = min([0.5*a, 0.5*b, 0.5*c])
        if 2.0*radius > a or 2.0 * radius > b or 2.0 * radius > c:
            m = "Radius {} cannot be greater than half the box ({},{},{})\n".format(radius, a, b, c)
            print(m) if self._logger is None else self._logger.error(m)
            exit()

        # Get pairs
        ncoms = len(self._com_list)
        idx = 0
        com_coords = np.zeros((ncoms, 3), dtype=np.float32)
        for jres, jcom in enumerate(self._com_list):
            com_coords[idx, :] = np.array(jcom)
            idx += 1
        pair_distsq_dict, pair_image_dict = calc_minimum_image_neighbors_com(com_coords, box, radius)

        # Write temporal pdbs containing monomers
        tmp_ob_nmols = []
        obconversion = ob.OBConversion()
        informat = 'pdb'
        outformat = 'mol2'
        obconversion.SetInAndOutFormats(informat, outformat)
        for ipair, _ in pair_image_dict.items():
            if onlydifferentmols:
                isequalmol_pair = self._equal_molecules_in_pair(path_pairs, ipair)
                if isequalmol_pair:
                    continue
            fpdbname = self._write_pdb_cluster(path_pairs, ipair, pair_image_dict[ipair])
            imol = ob.OBMol()
            obconversion.ReadFile(imol, fpdbname)
            tmp_ob_nmols.append(imol)
            self._npairs += 1

        # WRITE extracted monomers ===================================================
        obconversion = ob.OBConversion()
        outformat = 'mol2'
        obconversion.SetOutFormat(outformat)
        optimized_strings_mols = []
        for iconf in range(self._npairs):
            optimized_strings_mols.append(obconversion.WriteString(tmp_ob_nmols[iconf]))

        # Write Gaussian
        if write_gaussian:
            self._write_gaussian(localdir, optimized_strings_mols, "conformers",
                                 pattern=pattern, g16_key=g16_key, g16_extra_info=g16_extra_info,
                                 g16_mem=g16_mem, g16_nproc=g16_nproc,
                                 charge=charge, multiplicity=multiplicity)

        return None

    # #######################################################################
    def _extract_pairs_sphere_minatom(self, radius, write_gaussian=True, g16_key="#p 6-31g* mp2",
                                      g16_mem=4000, g16_nproc=4, pattern="conformers", g16_extra_info="",
                                      charge=0, multiplicity=1, onlydifferentmols=False):

        from MDAnalysis.analysis import distances

        # Write Pairs
        localdir = os.getcwd()
        path_pairs = os.path.join(localdir, self._pattern+"_g16_pairs")
        if os.path.isdir(path_pairs):
            shutil.rmtree(path_pairs)
        os.mkdir(path_pairs)

        # Box
        a = self._universe.dimensions[0]
        b = self._universe.dimensions[1]
        c = self._universe.dimensions[2]
        box = np.array([a, b, c], dtype=np.float32)

        # Check radius
        if radius is None:
            radius = min([0.5*a, 0.5*b, 0.5*c])
        if 2.0*radius > a or 2.0 * radius > b or 2.0 * radius > c:
            m = "Radius {} cannot be greater than half the box ({},{},{})\n".format(radius, a, b, c)
            print(m) if self._logger is None else self._logger.error(m)
            exit()

        # Get pairs
        pair_dist_dict = dict()
        pair_image_dict = dict()
        natoms = len(self._universe.atoms)
        idx = 0
        for imol in range(self._nmols):
            idx1_str = ""
            for i in self._nmols_array[imol]:
                idx1_str += "{} ".format(i)
            idx1_str += "\n"
            imol_select = self._universe.select_atoms('index {}'.format(idx1_str))
            for jmol in range(imol+1, self._nmols):
                idx2_str = ""
                for i in self._nmols_array[jmol]:
                    idx2_str += "{} ".format(i)
                idx2_str += "\n"
                jmol_select = self._universe.select_atoms('index {}'.format(idx2_str))
                dist_arr = distances.distance_array(imol_select.positions,
                                                    jmol_select.positions,
                                                    box=self._universe.dimensions)
                minimum_distance = np.min(dist_arr)
                minimum_distance_indexes = np.unravel_index(np.argmin(dist_arr, axis=None), dist_arr.shape)
                if minimum_distance < radius:
                    pair_dist_dict[(imol, jmol)] = minimum_distance
                    i1 = list(minimum_distance_indexes)[0]
                    i2 = list(minimum_distance_indexes)[1]
                    iatom1 = self._nmols_array[imol][i1]
                    iatom2 = self._nmols_array[jmol][i2]
                    dx = self._universe.atoms[iatom1].position[0] - self._universe.atoms[iatom2].position[0]
                    dy = self._universe.atoms[iatom1].position[1] - self._universe.atoms[iatom2].position[1]
                    dz = self._universe.atoms[iatom1].position[2] - self._universe.atoms[iatom2].position[2]
                    imgx = 0; imgy = 0; imgz = 0
                    if dx > box[0] * 0.5: imgx = -1
                    if dx <= -box[0] * 0.5: imgx = +1
                    if dy > box[1] * 0.5:  imgy = -1
                    if dy <= -box[1] * 0.5: imgy = +1
                    if dz > box[2] * 0.5: imgz = -1
                    if dz <= -box[2] * 0.5: imgz = +1
                    pair_image_dict[(imol, jmol)] = [imgx, imgy, imgz]

        # Write temporal pdbs containing monomers
        tmp_ob_nmols = []
        obconversion = ob.OBConversion()
        informat = 'pdb'
        outformat = 'mol2'
        obconversion.SetInAndOutFormats(informat, outformat)
        for ipair, _ in pair_dist_dict.items():
            if onlydifferentmols:
                isequalmol_pair = self._equal_molecules_in_pair(path_pairs, ipair)
                if isequalmol_pair:
                    continue
            fpdbname = self._write_pdb_cluster(path_pairs, ipair, pair_image_dict[ipair])
            imol = ob.OBMol()
            obconversion.ReadFile(imol, fpdbname)
            tmp_ob_nmols.append(imol)
            formula = imol.GetFormula()
            fpdbold = fpdbname
            tmp = os.path.split(fpdbname)[-1]
            fpdbnew = os.path.join(path_pairs, formula+"_"+tmp)
            os.rename(fpdbold, fpdbnew)
            query = ob.CompileMoleculeQuery(imol)
            self._npairs += 1

        # WRITE extracted monomers ===================================================
        obconversion = ob.OBConversion()
        outformat = 'mol2'
        obconversion.SetOutFormat(outformat)
        optimized_strings_mols = []
        formula_string_mols = []
        for iconf in range(self._npairs):
            optimized_strings_mols.append(obconversion.WriteString(tmp_ob_nmols[iconf]))
            formula_string_mols.append(tmp_ob_nmols[iconf].GetFormula())

        # Write Gaussian
        if write_gaussian:
            self._write_gaussian(localdir, optimized_strings_mols, formula_string_mols, "conformers",
                                 pattern=pattern, g16_key=g16_key, g16_extra_info=g16_extra_info,
                                 g16_mem=g16_mem, g16_nproc=g16_nproc,
                                 charge=charge, multiplicity=multiplicity)
        return None

    # #######################################################################
    def _write_gaussian(self, localdir, optimized_strings_mols, formula_string_mol, typepair, pattern="conformers",
                        g16_key="#p 6-31g* mp2", g16_extra_info="",
                        g16_mem=4000, g16_nproc=4, charge=0, multiplicity=1):

        # Create directory for gaussian inputs
        if self._extract_monomer:
            typepair = "monomers"
        elif self._extract_pair:
            typepair = "pairs"
        elif self._extract_trimer:
            typepair = "trimers"

        parent_dir = os.path.join(localdir, "{}".format("{}_g16_{}/").format(pattern, typepair))
        if not os.path.isdir(parent_dir):
            os.mkdir(parent_dir)

        nconf = len(optimized_strings_mols)

        for item in range(nconf):
            iconf_seed = item
            energy = 0.0

            obconversion = ob.OBConversion()
            obconversion.SetInAndOutFormats('mol2', 'xyz')
            obmol = ob.OBMol()
            xyz_obmol = obconversion.ReadString(obmol, optimized_strings_mols[iconf_seed])
            if not xyz_obmol:
                m = "Something is wrong with Gaussian16 input files\n"
                m += "Gaussian input files are not written"
                print(m) if self._logger is None else self._logger.error(m)
                return None

            xyz_list = obconversion.WriteString(obmol).split("\n")
            legacypattern = os.path.splitext(os.path.split(xyz_list[1])[-1])[0]
            fname = os.path.join(parent_dir, formula_string_mol[item]+"_"+legacypattern+".com")
            fchk_name = formula_string_mol[item]+"_"+legacypattern+".chk"
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

                for item in g16_extra_info:
                    if not item.strip() == '':
                        f.writelines(item)
                        f.writelines("\n")
                f.writelines("\n")

    # #######################################################################
    def _write_pdb_cluster(self, path, listmols, image):

        if len(listmols) == 1:
            tmp = "monomer_{0:03d}.pdb".format(listmols[0])
            filename_pdb = os.path.join(path, tmp)
            natoms_cluster = len(self._nmols_array[listmols[0]])
        elif len(listmols) == 2:
            tmp = "pair_{0:03d}_{1:03d}.pdb".format(listmols[0], listmols[1])
            filename_pdb = os.path.join(path, tmp)
            natoms_cluster = len(self._nmols_array[listmols[0]]) + \
                             len(self._nmols_array[listmols[1]])
        elif len(listmols) == 3:
            tmp = "trimer_{0:03d}_{1:03d}_{2:03d}.pdb".format(listmols[0], listmols[1], listmols[2])
            filename_pdb = os.path.join(path, tmp)
            natoms_cluster = len(self._nmols_array[listmols[0]]) + \
                             len(self._nmols_array[listmols[1]]) + \
                             len(self._nmols_array[listmols[2]])
        else:
            filename_pdb = "cluster.pdb"

        fmt = {
            'ATOM': (
                "ATOM  {serial:5d} {name:<4s}{altLoc:<1s}{resName:<4s}"
                "{chainID:1s}{resSeq:4d}{iCode:1s}"
                "   {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}{occupancy:6.2f}"
                "{tempFactor:6.2f}      {segID:<4s}{element:>2s}\n"),
            'HETATM': (
                "HETATM{serial:5d} {name:<4s}{altLoc:<1s}{resName:<4s}"
                "{chainID:1s}{resSeq:4d}{iCode:1s}"
                "   {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}{occupancy:6.2f}"
                "{tempFactor:6.2f}      {segID:<4s}{element:>2s}\n"),
            'REMARK': "REMARK     {0}\n",
            'COMPND': "COMPND    {0}\n",
            'HEADER': "HEADER    {0}\n",
            'TITLE': "TITLE     {0}\n",
            'MODEL': "MODEL     {0:>4d}\n",
            'NUMMDL': "NUMMDL    {0:5d}\n",
            'ENDMDL': "ENDMDL\n",
            'END': "END\n",
            'CRYST1': ("CRYST1{0:9.3f}{1:9.3f}{2:9.3f}"
                       "{3:7.2f}{4:7.2f}{5:7.2f} "
                       "{6:<11s}{7:4d}\n"),
            'CONECT': "CONECT{0}\n"
        }

        a = self._universe.dimensions[0]
        b = self._universe.dimensions[1]
        c = self._universe.dimensions[2]
        alpha = self._universe.dimensions[3]
        beta = self._universe.dimensions[4]
        gamma = self._universe.dimensions[5]

        atom_kind_molecule_label = None
        spacegroup = "P -1"
        zvalue = 1
        radtodeg = 180 / np.pi

        if len(listmols) == 2:
            listmols, image = self._check_order_mols_in_pairs(path, listmols, image)

        # Minimum Image Convention for all atoms in the lists
        ich_local = 0
        idx_local = 0
        coord_cluster = np.zeros([natoms_cluster, 3], dtype=np.float32)
        for ich in listmols:
            for idx in self._nmols_array[ich]:
                if ich_local == 0:
                    coord_cluster[idx_local, :] = self._universe.coord[idx]
                else:
                    for j in range(0, 3):
                        coord_cluster[idx_local, j] = self._universe.coord[idx][j] - \
                                                      image[j]*self._universe.dimensions[j]

                idx_local += 1
            ich_local += 1

        with open(filename_pdb, 'w') as fpdb:

            fpdb.write(fmt['REMARK'].format('Created with GeCos(J.Ramos)'))
            fpdb.write(fmt['CRYST1'].format(a, b, c, alpha, beta, gamma,
                                            spacegroup, zvalue))

            idx_local = 0
            for ich in listmols:
                for idx in self._nmols_array[ich]:

                    resname = "{0:03d}".format(ich)

                    fpdb.write(fmt['HETATM'].format(
                        serial=idx_local + 1,
                        name=self._universe.atoms[idx].name,
                        altLoc=" ",
                        resName=resname,
                        chainID=" ",
                        resSeq=idx_local + 1,
                        iCode=" ",
                        pos=[i for i in coord_cluster[idx_local, :]],
                        occupancy=1.0,
                        tempFactor=1.0,
                        segID="    ",
                        element=self._universe.atoms[idx].element
                    ))

                    idx_local += 1

            fpdb.write('END\n')

        return filename_pdb

    # #######################################################################
    def _write_tcl_atoms_com(self, filename_tcl="vmd_com.tcl"):

        lines = "proc newRep { sel type color rep imol} {\n"
        lines += "    mol selection $sel\n"
        lines += "    mol representation $type\n"
        lines += "    mol addrep $imol\n"
        lines += "    mol showrep $imol $rep on\n"
        lines += "    mol modcolor $rep $imol $color\n"
        lines += "}\n"
        lines += "\n"
        lines += "set dir \"{}\"\n".format(os.path.split(filename_tcl)[0])
        lines += "\n"
        lines += "display projection orthographic\n"
        lines += "axes location off\n"
        lines += "color Display Background white\n"
        lines += "display depthcue off\n"
        lines += "\n"

        lines += "\n"
        lines += "mol new {} type pdb\n".format(self._filename)
        lines += "set imol1 [molinfo top]\n"
        lines += 'mol delrep 0 $imol1\n'
        lines += 'set rep1 0\n'
        lines += 'newRep "all" "CPK" "Name" $rep1 $imol1\n'
        lines += '\n'
        lines += "mol new {} type pdb\n".format(self._com_filename_pdb)
        lines += "set imol2 [molinfo top]\n"
        lines += 'mol delrep 0 $imol2\n'
        lines += 'set rep2 0\n'
        lines += 'newRep "all" "vdw 0.6" "Index" $rep2 $imol2\n'
        lines += '\n'
        lines += "pbc box\n"
        lines += '\n'

        with open(filename_tcl, 'w') as ftcl:
            ftcl.writelines(lines)

    # #######################################################################
    def _equal_molecules_in_pair(self, path, ipair):

        """
        Using the isomorphism code in OpenBabel find if two molecules in a pair are equals.

        :param path:
        :param ipair:
        :return:
        """

        obconversion = ob.OBConversion()

        # Fragment 1 from the pair
        fpdbname = self._write_pdb_cluster(path, [ipair[0]], [0, 0, 0])
        imol1 = ob.OBMol()
        obconversion.ReadFile(imol1, fpdbname)
        os.remove(fpdbname)

        # Fragment 2 from the pair
        fpdbname = self._write_pdb_cluster(path, [ipair[1]], [0, 0, 0])
        imol2 = ob.OBMol()
        obconversion.ReadFile(imol2, fpdbname)
        os.remove(fpdbname)

        # Isomorphic structure?
        query = ob.CompileMoleculeQuery(imol1)
        mapper = ob.OBIsomorphismMapper.GetInstance(query)
        isomorph = ob.vpairUIntUInt()
        mapper.MapFirst(imol2, isomorph)

        if len(list(isomorph)) == 0:
            return False    # Different molecules
        else:
            return True    # Equal molecules

    # #######################################################################
    def _check_order_mols_in_pairs(self, path, ipair, image):

        """
        Using the isomorphism code in OpenBabel find if two molecules in a pair are equals.

        :param path:
        :param ipair:
        :param image
        :return:
        """

        obconversion = ob.OBConversion()
        # Fragment 1 from the pair
        fpdbname = self._write_pdb_cluster(path, [ipair[0]], [0, 0, 0])
        imol1 = ob.OBMol()
        obconversion.ReadFile(imol1, fpdbname)
        os.remove(fpdbname)
        # Fragment 2 from the pair
        fpdbname = self._write_pdb_cluster(path, [ipair[1]], [0, 0, 0])
        imol2 = ob.OBMol()
        obconversion.ReadFile(imol2, fpdbname)
        os.remove(fpdbname)

        if self._imol1_ref is None and self._imol2_ref is None:

            self._imol1_ref = imol1
            self._imol2_ref = imol2
            self._ipair_ref = ipair
            return [ipair[0], ipair[1]], image
        else:
            query = ob.CompileMoleculeQuery(self._imol1_ref)
            mapper = ob.OBIsomorphismMapper.GetInstance(query)
            isomorph = ob.vpairUIntUInt()
            mapper.MapFirst(imol1, isomorph)
            if len(list(isomorph)) > 0:
                return [ipair[0], ipair[1]], image
            else:
                tup = [ipair[1], ipair[0]]
                image = [-1*i for i in image]
                return tup, image
