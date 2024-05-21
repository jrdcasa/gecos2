import os
import datetime
import glob
from openbabel import openbabel as ob


class GecosPdbFolder(object):

    # #######################################################################
    def __init__(self, pdbfolder, logger=None):

        self._pdbfolder = pdbfolder
        self._logger = logger
        # Number of pdbs
        self._list_pdb_files = sorted(glob.glob(os.path.join(self._pdbfolder, "*.pdb")))
        self._npdb_files = len(self._list_pdb_files)

        filename = self._list_pdb_files[0]
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

        try:
            self._mol_pybabel = mols[0]
        except IndexError:
            self._mol_pybabel = None

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\t\t******** CREATION OF THE GECOS PDB FOLDER OBJECT ********\n\n"
        m += "\t\tConformers will be extracted from a frame of a MD simulation.\n"
        m += "\t\t Find number of molecules in the frame.({})\n".format(now)
        m += "\n"
        m += "\t\tFormula             = {}\n".format(self._mol_pybabel.GetFormula())
        m += "\t\tNumber of atoms     = {}\n".format(self._mol_pybabel.NumAtoms())
        m += "\t\tMolecular weight    = {0:.2f} g/mol\n".format(self._mol_pybabel.GetMolWt())
        m += "\t\tNumber of rings     = {}\n".format(len(self._mol_pybabel.GetSSSR()))
        m += "\t\tNumber of molecules = {}\n".format(self._npdb_files)
        m += "\t\t******** CREATION OF THE GECOS PDB FOLDER OBJECT ********\n\n"
        print(m) if self._logger is None else self._logger.info(m)

    # #######################################################################
    def generate_pdbfolder(self, write_gaussian=True, g16_key="#p 6-31g* mp2",
                           g16_mem=4000, g16_extra_info="",
                           g16_nproc=4, pattern="conformers", charge=0, multiplicity=1):

        # Write Monomoners
        localdir = os.getcwd()

        # Write temporal pdbs containing monomers
        tmp_ob_nmols = []
        obconversion = ob.OBConversion()
        informat = 'pdb'
        outformat = 'mol2'
        obconversion.SetInAndOutFormats(informat, outformat)
        for ipdb in range(0, self._npdb_files):
            fpdbname = self._list_pdb_files[ipdb]
            imol = ob.OBMol()
            obconversion.ReadFile(imol, fpdbname)
            tmp_ob_nmols.append(imol)

        # WRITE extracted conformers ===================================================
        obconversion = ob.OBConversion()
        outformat = 'mol2'
        obconversion.SetOutFormat(outformat)
        optimized_strings_mols = []
        formula_string_mols = []
        for iconf in range(self._npdb_files):
            optimized_strings_mols.append(obconversion.WriteString(tmp_ob_nmols[iconf]))
            formula_string_mols.append(tmp_ob_nmols[iconf].GetFormula())

        # Write Gaussian
        if write_gaussian:
            self._write_gaussian(localdir, optimized_strings_mols, formula_string_mols, "conformers",
                                 pattern=pattern, g16_key=g16_key, g16_extra_info=g16_extra_info,
                                 g16_mem=g16_mem, g16_nproc=g16_nproc,
                                 charge=charge, multiplicity=multiplicity)

    # #######################################################################
    def _write_gaussian(self, localdir, optimized_strings_mols, formula_string_mol, typepair, pattern="conformers",
                        g16_key="#p 6-31g* mp2", g16_extra_info="",
                        g16_mem=4000, g16_nproc=4, charge=0, multiplicity=1):

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

                for jtem in g16_extra_info:
                    if not jtem.strip() == '':
                        f.writelines(jtem)
                        f.writelines("\n")
                f.writelines("\n")
