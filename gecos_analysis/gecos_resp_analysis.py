import glob
import os.path
from collections import defaultdict
from openbabel import openbabel as ob
import numpy as np
from utils.resp_files_pyresp_server import system_config, project_config


class GecosRespAnalysis:

    """
    This class is use to analyze RESP calculations from PYRED server
    (https://upjv.q4md-forcefieldtools.org/REDServer-Development/).
    """

    # =========================================================================
    def __init__(self, workdir, gaussian_gecos_obj, is_prepareinputs, resp_software,
                 deltaenergythreshold=None, nconfsthreshold=None, logger=None):

        """
        Initialize the instance
        """

        self._workdir = workdir
        self._logger = logger

        if nconfsthreshold is not None and deltaenergythreshold is not None:
            self._nconfsthreshold = nconfsthreshold
            self._deltaenergythreshold = None
            m = "\n\t\t WARNING: Energy and configuration thresholds are activated.\n"
            m += "\t\t            These parameters cannot be used at the same time.\n"
            m += "\t\t Configuration threshold activated by default.\n"
            m += "\t\t Configuration threshold = {}.\n".format(nconfsthreshold)
        elif nconfsthreshold is not None:
            self._nconfsthreshold = nconfsthreshold
            self._deltaenergythreshold = None
            m = "\t\t Configuration threshold activated by default.\n"
            m += "\t\t Configuration threshold = {}.\n".format(nconfsthreshold)
        elif deltaenergythreshold is not None:
            self._nconfsthreshold = None
            self._deltaenergythreshold = deltaenergythreshold
            m = "\t\t Energy threshold activated by default.\n"
            m += "\t\t Energy threshold = {} kcal/mol.\n".format(deltaenergythreshold)
        print(m) if self._logger is None else self._logger.info(m)

        self._isprepareinputs = is_prepareinputs
        self._gaussian_gecos_obj = gaussian_gecos_obj

        if self._isprepareinputs and resp_software.upper() == "FORPYRED":
            self._writepdb_forpyred_server()
        elif self._isprepareinputs and resp_software.upper() == "MULTIWFN":
            self._write_script_for_multiwfn_resp()

    # =========================================================================
    def _writepdb_forpyred_server(self):

        with open(os.path.join("./", "Mol_red1.pdb"), 'w') as fpdb:
            for idx, item in enumerate(self._gaussian_gecos_obj._xyzlist):
                if self._deltaenergythreshold is not None:
                    if self._gaussian_gecos_obj._df['DeltaEnergy(kcal/mol)'][idx] > self._deltaenergythreshold:
                        continue
                if self._nconfsthreshold is not None:
                    if idx >= self._nconfsthreshold :
                        continue
                igeomxyz = item
                molref = ob.OBMol()
                obref = ob.OBConversion()
                obref.SetInAndOutFormats('xyz', 'pdb')
                obref.ReadString(molref, igeomxyz)
                namefileref = os.path.join("./", "tmp_resp.pdb")
                obref.WriteFile(molref, namefileref)

                line = "MODEL {}\n".format(idx+1)

                with open(namefileref, 'r') as fref:
                    linestmp = fref.readlines()
                    for iline in linestmp:
                        if iline.count("HETATM") != 0 or iline.count("ATOM") != 0:
                            line += iline
                line += "ENDMDL\n".format(idx+1)

                fpdb.writelines(line)

        with open("System.config", 'w') as fsys:
            fsys.writelines(system_config)

        with open("Project.config", 'w') as fsys:
            fsys.writelines(project_config)

    # =========================================================================
    def _write_script_for_multiwfn_resp(self, temp=450):

        fconflist = open("conflist.txt", 'w')
        fmultiwfn_script = open("input_multiwfn.txt", 'w')
        population_dict = defaultdict()

        os.path.join(self._workdir, "*.wfn")
        lwfn = glob.glob(os.path.join(self._workdir, "*.wfn"))
        lfchk = glob.glob(os.path.join(self._workdir, "*.fchk"))

        if len(lwfn) > 0:
            ext = ".wfn"
        elif len(lfchk) > 0:
            ext = ".fchk"
        else:
            m = "\n\t\t ERROR: No wfn or fchk files in the working directory.\n"
            print(m) if self._logger is None else self._logger.error(m)
            exit()

        rt = temp * 1.987e-03  # kcal/mol
        sum = 0.0
        for idx, item in enumerate(self._gaussian_gecos_obj._xyzlist):
            if self._deltaenergythreshold is not None:
                if self._gaussian_gecos_obj._df['DeltaEnergy(kcal/mol)'][idx] > self._deltaenergythreshold:
                    continue
            if self._nconfsthreshold is not None:
                if idx >= self._nconfsthreshold:
                    continue

            filepattern = os.path.splitext(os.path.split(self._gaussian_gecos_obj._logfiles[idx])[-1])[0]

            if idx == 0:
                linemultiwfn = filepattern + ext + "\n"
                linemultiwfn += "7\n"
                linemultiwfn += "18\n"
                linemultiwfn += "5\n"
                linemultiwfn += "1\n"
                linemultiwfn += "eqvcons.txt\n"
                linemultiwfn += "6\n"
                linemultiwfn += "1\n"
                linemultiwfn += "chgcons.txt\n"
                linemultiwfn += "-1\n"
                linemultiwfn += "conflist.txt\n"
                linemultiwfn += "2\n"
                linemultiwfn += "EOF\n"
                fmultiwfn_script.writelines(linemultiwfn)


            try:
                if np.isnan(self._gaussian_gecos_obj._df['Delta_G(kcal/mol)'][idx]):
                    population_dict[filepattern] = self._gaussian_gecos_obj._df['DeltaEnergy(kcal/mol)'][idx]
                else:
                    population_dict[filepattern] = self._gaussian_gecos_obj._df['Delta_G(kcal/mol)'][idx]
            except AttributeError:
                population_dict[filepattern] = self._gaussian_gecos_obj._df['DeltaEnergy(kcal/mol)'][idx]

            sum += np.exp(-population_dict[filepattern] / rt)

        line = ""
        for key, value in population_dict.items():
            term = np.exp(-value / rt)
            population_dict[key] = term/sum
            line += "./{0:s} {1:10.4f}\n".format(key+ext, population_dict[key])
        fconflist.writelines(line)

        m = "\t\t Input file for MultiWfn software has been saved in {}\n".format("input_multiwfn.txt")
        m += "\t\t and conformer populations have been saved in {}\n".format("conflist.txt")
        m += "\t\t You need to create the followring files by yourself:\n"
        m += "\t\t\t * eqvcons.txt --> Equivalence constraint (see 4.7.7.3 section in Multiwfn manual).\n"
        m += "\t\t\t * chgcons.txt --> Charge constraint (see 4.7.7.4 section in Multiwfn manual).\n"
        print(m) if self._logger is None else self._logger.info(m)

        fconflist.close()
        fmultiwfn_script.close()