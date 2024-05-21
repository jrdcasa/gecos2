import cclib
import glob
import os
from collections import defaultdict
from openbabel import openbabel as ob
import pandas as pd
import rmsd
import copy
import numpy as np
import datetime
import MDAnalysis
from utils.internal_coordinates import dihedral_py
from utils.atomic_data import atomic_number_to_element
try:
    from gecos_analysis.gecos_closecontacts import CloseContacts
except ModuleNotFoundError:
    from gecos_closecontacts import CloseContacts


class GaussianGecos:

    # =========================================================================
    def __init__(self, pathlogfiles, ext="log", logger=None):

        """
        Initialize the instance

        :param pathlogfiles: Path to the output files from Gaussian16
        :param ext: Extension of the output files
        :param logger: Logger to write results
        """

        # Logger to write results
        self._logger = logger
        # Path where the log iles are located
        self._logfiles = sorted(glob.glob(os.path.join(pathlogfiles, "*." + ext)))
        self._optscfenergies = defaultdict()
        self._freenergies = defaultdict()
        self._optgeometry = defaultdict(list)
        self._rotconsts = defaultdict(list)
        self._vibfreqs = defaultdict()
        self._vibirs = defaultdict()
        self._temperature = defaultdict()
        self._xyzlist = []
        self._internalcoords = defaultdict(dict)
        self._minscfenergy_ha = None
        self._2dgrid = None
        self._nameID_to_idx = defaultdict()
        self._idx_to_nameID = defaultdict()
        self._cre_dict = dict()

        if len(self._logfiles) == 0:
            m = "\n\t\t ERROR. No G16 log files in the folder\n"
            m += "\t\t {}".format(pathlogfiles)
            print(m) if self._logger is None else self._logger.error(m)
            exit()

        # DataFrame
        self._df = pd.DataFrame({'ID': [],
                                 'IDCluster': [],
                                 'Type': [],
                                 'OptEnergy(Ha)': [],
                                 'DeltaEnergy(kcal/mol)': [],
                                 'Path': [],
                                 'RMSDwithHs': [],
                                 'RMSDwithoutHs': [],
                                 'RotConstA(Ghz)': [],
                                 'RotConstB(Ghz)': [],
                                 'RotConstC(Ghz)': [],
                                 'AngleMainPrincipalAxes(deg)': [],
                                 'RadiusGyration(ang)': [],
                                 'NHbonds': [],
                                 'NvdwContacts_Intermol': [],
                                 'NvdwContacts_Intramol': [],
                                 'Delta_G(kcal/mol)': [],
                                 'Delta_-TS(kcal/mol)': [],
                                 'Delta_H(kcal/mol)': []})

    # =========================================================================
    def getlogfiles(self):

        return self._logfiles

    # =========================================================================
    def getvibfreqs(self):

        return self._vibfreqs

    # =========================================================================
    def getvibirs(self):

        return self._vibirs

    # =========================================================================
    def gettemperature(self):

        return self._temperature

    # =========================================================================
    def getdeltag(self):

        dd = defaultdict()
        ll = list(self._df['Delta_G(kcal/mol)'])
        for i, item in enumerate(list(self._df['ID'])):
            dd[item] = ll[i]

        return dd

    # =========================================================================
    def getoptxyz_coordinates(self):

        return self._xyzlist

    # =========================================================================
    def extract_vibrational_ir(self):

        # ==================== THERMOCHEMISTRY ============================
        # Parse data from gaussian logs
        tmp_optgeometry = defaultdict(list)
        entropy_list = defaultdict()
        temperature_list = defaultdict()
        enthalpy_list = defaultdict()
        scf_list = defaultdict()
        vibfreqs_list = defaultdict()
        vibirs_list = defaultdict()
        nl = len(self._logfiles)
        for idx, ipath in enumerate(self._logfiles):

            # TODO: Error if there is only one log file. Check if this is correct
            if nl >= 10 and idx % int(nl * 0.1) == 0:
                now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
                m = "\t\t\tParsing {} log files out of {} ({})".format(idx, nl, now)
                print(m) if self._logger is None else self._logger.info(m)
            elif nl <= 10:
                now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
                m = "\t\t\tParsing {} log files out of {} ({})".format(idx, nl, now)
                print(m) if self._logger is None else self._logger.info(m)
            ifile = os.path.splitext(os.path.split(ipath)[-1])[0]
            parser = cclib.io.ccopen(ipath)
            data = parser.parse()
            self._freenergies[ifile] = data.freeenergy  # hartrees
            entropy_list[ifile] = data.entropy   # hartree/K
            enthalpy_list[ifile] = data.enthalpy  # hartrees
            temperature_list[ifile] = data.temperature  # K
            scf_list[ifile] = data.scfenergies[-1]  # in eV
            vibfreqs_list[ifile] = data.vibfreqs  # cm-1
            vibirs_list[ifile] = data.vibirs

            tmp_optgeometry[ifile].append(data.atomnos)
            tmp_optgeometry[ifile].append(data.atomcoords[-1])
            self._rotconsts[ifile].append(data.rotconsts[-1])

        # Ordering the structures by free energy
        self._freenergies = dict(sorted(self._freenergies.items(), key=lambda v: v[1]))
        pathlist = [i for i in self._freenergies.keys()]
        idlist = [os.path.splitext(os.path.split(i)[-1])[0] for i in self._freenergies.keys()]
        freeenergylist = [i for key, i in self._freenergies.items()]
        minfreeenergy = freeenergylist[0]

        deltafreeenergylist = [cclib.parser.utils.convertor(i - minfreeenergy, 'hartree', 'kcal/mol')
                               for i in freeenergylist]

        # Get all values of thermochemistry
        iref = list(self._freenergies.items())[0][0]
        self._df['ID'] = idlist
        self._df['Delta_G(kcal/mol)'] = deltafreeenergylist
        self._df['Path'] = pathlist
        delta_tentropy_list = []
        delta_enthalpy_list = []
        delta_scf_list = []
        rotconstlist_a = []
        rotconstlist_b = []
        rotconstlist_c = []

        for ikey, ivalue in self._freenergies.items():
            delta_tentropy = entropy_list[ikey] - entropy_list[iref]
            delta_tentropy = cclib.parser.utils.convertor(delta_tentropy, 'hartree', 'kcal/mol')\
                             * temperature_list[iref] * (-1)
            delta_tentropy_list.append(delta_tentropy)

            delta_enthalpy = enthalpy_list[ikey] - enthalpy_list[iref]
            delta_enthalpy = cclib.parser.utils.convertor(delta_enthalpy, 'hartree', 'kcal/mol')
            delta_enthalpy_list.append(delta_enthalpy)

            delta_scf = scf_list[ikey] - scf_list[iref]
            delta_scf = cclib.parser.utils.convertor(delta_scf, 'eV', 'kcal/mol')
            delta_scf_list.append(delta_scf)

            self._vibfreqs[ikey] = vibfreqs_list[ikey]
            self._vibirs[ikey] = vibirs_list[ikey]
            self._temperature[ikey] = temperature_list[ikey]

            rotconstlist_a.append(self._rotconsts[ikey][0][0])
            rotconstlist_b.append(self._rotconsts[ikey][0][1])
            rotconstlist_c.append(self._rotconsts[ikey][0][2])
            self._optgeometry[ikey] = tmp_optgeometry[ikey]

        self._df['Delta_-TS(kcal/mol)'] = delta_tentropy_list   # kcal/mol
        self._df['Delta_H(kcal/mol)'] = delta_enthalpy_list  # kcal/mol
        self._df['Delta_E(kcal/mol)'] = delta_scf_list  # kcal/mol
        self._df['RotConstA(Ghz)'] = rotconstlist_a
        self._df['RotConstB(Ghz)'] = rotconstlist_b
        self._df['RotConstC(Ghz)'] = rotconstlist_c
        # Hartrees, In this case, the free energy instead of the scf one is used in this case
        self._df['OptEnergy(Ha)'] = freeenergylist

        # Write all conformers in a xyx file
        with open("Optimized_Conformers_sorted.xyz", 'w') as fxyz:
            for key, value in self._optgeometry.items():
                natoms = len(value[0])
                line = "{0:<10d}\n".format(natoms)
                line += "{0:<s}\n".format(key)
                for idx in range(0, natoms):
                    atomic_number = atomic_number_to_element[value[0][idx]-1]
                    x, y, z = value[1][idx]
                    line += "{0:<3s} {1:<12.6f} {2:<12.6f} {3:<12.6f}\n".\
                        format(atomic_number, float(x), float(y), float(z))
                fxyz.writelines(line)

    # =========================================================================
    def extract_energy(self):

        # Extracting energy of the optimized structure from the logfiles
        tmp_optgeometry = defaultdict(list)
        nl = len(self._logfiles)
        for idx, ifile in enumerate(self._logfiles):
            # TODO: Error if there is only one log file. Check if this is correct
            if nl >= 10 and idx % int(nl*0.1) == 0:
                now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
                m = "\t\t\tParsing {} log files out of {} ({})".format(idx, nl, now)
                print(m) if self._logger is None else self._logger.info(m)
            elif nl <= 10:
                now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
                m = "\t\t\tParsing {} log files out of {} ({})".format(idx, nl, now)
                print(m) if self._logger is None else self._logger.info(m)

            parser = cclib.io.ccopen(ifile)
            data = parser.parse()
            self._optscfenergies[ifile] = data.scfenergies[-1]  # in eV,
            tmp_optgeometry[ifile].append(data.atomnos)
            tmp_optgeometry[ifile].append(data.atomcoords[-1])
            self._rotconsts[ifile].append(data.rotconsts[-1])

        # Creating the dataframe
        self._optscfenergies = dict(sorted(self._optscfenergies.items(), key=lambda v: v[1]))
        pathlist = [i for i in self._optscfenergies.keys()]
        idlist = [os.path.splitext(os.path.split(i)[-1])[0] for i in self._optscfenergies.keys()]
        energylist = [cclib.parser.utils.convertor(i, 'eV', 'hartree') for key, i in self._optscfenergies.items()]
        minenergy = energylist[0]
        self._minscfenergy_ha = minenergy
        deltaenergylist = [cclib.parser.utils.convertor(i - minenergy, 'hartree', 'kcal/mol') for i in energylist]

        rotconstlist_a = []
        rotconstlist_b = []
        rotconstlist_c = []

        for ikey, ivalue in self._optscfenergies.items():
            rotconstlist_a.append(self._rotconsts[ikey][0][0])
            rotconstlist_b.append(self._rotconsts[ikey][0][1])
            rotconstlist_c.append(self._rotconsts[ikey][0][2])
            self._optgeometry[ikey] = tmp_optgeometry[ikey]

        self._df['ID'] = idlist
        self._df['OptEnergy(Ha)'] = energylist
        self._df['DeltaEnergy(kcal/mol)'] = deltaenergylist
        self._df['Path'] = pathlist
        self._df['RotConstA(Ghz)'] = rotconstlist_a
        self._df['RotConstB(Ghz)'] = rotconstlist_b
        self._df['RotConstC(Ghz)'] = rotconstlist_c

    # =========================================================================
    def extract_rmsd(self):

        """
        Calculate the rmsd between all molecules and the lowest energy molecule
        as reference

        """

        rmsd_all_list = []
        rmsd_noh_list = []
        idxfile = 0
        for ifile, values in self._optgeometry.items():
            localdir, basename = os.path.split(ifile)
            basename = os.path.splitext(basename)[0]
            natoms = len(values[0])
            xyzstring = "{0:<s}\n{1:<s}\n".format(str(natoms), basename)
            for idx in range(len(values[0])):
                element = ob.GetSymbol(int(values[0][idx]))
                xyzstring += "{0:<s}  {1:10.4f}  {2:10.4f}  {3:10.4f}\n".format \
                    (element, values[1][idx][0], values[1][idx][1], values[1][idx][2])
            self._xyzlist.append(xyzstring)
            # xyzlist[0] = "3\nCarb\nC 1.54 0.0 0.0\nC 0.0 1.54 0.0\nC 0.0 0.0 1.54"
            if idxfile == 0:
                molref = ob.OBMol()
                obref = ob.OBConversion()
                obref.SetInAndOutFormats('xyz', 'xyz')
                obref.ReadString(molref, self._xyzlist[0])
                namefileref = os.path.join(localdir, "reference_allign.xyz")
                obref.WriteFile(molref, namefileref)
                atoms_ref, coord_ref = rmsd.get_coordinates_xyz(namefileref)
                atoms_ref_noh = np.where(atoms_ref != 'H')
                coord_ref_noh = copy.deepcopy(coord_ref[atoms_ref_noh])
                rmsd_all_list.append(0.0)
                rmsd_noh_list.append(0.0)
                os.remove(namefileref)
            else:
                obtarget = ob.OBConversion()
                obtarget.SetInAndOutFormats('xyz', 'xyz')
                moltarget = ob.OBMol()
                obtarget.ReadString(moltarget, self._xyzlist[idxfile])
                a = ob.OBAlign(False, False)
                a.SetRefMol(molref)
                a.SetTargetMol(moltarget)
                a.Align()
                a.UpdateCoords(moltarget)
                namefile = os.path.join(localdir, "target_{}_allign.xyz".format(idxfile))
                obref.WriteFile(moltarget, namefile)
                atoms_target, coord_target = rmsd.get_coordinates_xyz(namefile)
                atoms_target_noh = np.where(atoms_target != 'H')
                coord_target_noh = copy.deepcopy(coord_target[atoms_target_noh])
                rmsd_noh_list.append(rmsd.rmsd(coord_target_noh, coord_ref_noh))
                rmsd_all_list.append(rmsd.rmsd(coord_target, coord_ref))
                os.remove(namefile)
            idxfile += 1

        self._df["RMSDwithHs"] = rmsd_all_list
        self._df["RMSDwithoutHs"] = rmsd_noh_list

    # =========================================================================
    def extract_internalcoords(self, args):

        """
         The ndx file must contiain labels [ distances ], [ angles ] and/or [ dihedral ].
         Following the label a line for each the internal coordinate must be listed.
         The internal coodinate must to be present in the log file otherwise is not considered.
         Example:
              [ distances ]
              1 2
              2 3
              [ dihedrals ]
              4 3 2 1
              5 4 3 1
         This file extracts the values of two distances (1-2 and 2-3) and two dihedrals (4-3-2-1 and 5-4-3-1)
         The numbering starts at 1 in th file internally the numbering starts at 0.
        """

        distances_list = []
        angles_list = []
        dihedrals_list = []

        with open(args.indxfile, 'r') as fout:
            lines = fout.readlines()
            for iline in lines:
                if iline.count("[") != 0:
                    continue
                else:
                    tokens = iline.split()
                    if len(tokens) == 2:
                        # Distance
                        try:
                            distances_list.append([int(i) for i in tokens])
                        except (Exception, ):
                            m = "\t\t ERROR in {} file. Distances".format(args.indxfile)
                            print(m) if self._logger is None else self._logger.info(m)
                            exit()
                    elif len(tokens) == 3:
                        # Angles
                        try:
                            angles_list.append([int(i) for i in tokens])
                        except (Exception, ):
                            m = "\t\t ERROR in {} file. Angles".format(args.indxfile)
                            print(m) if self._logger is None else self._logger.info(m)
                            exit()
                    elif len(tokens) == 4:
                        # Dihedrals
                        try:
                            dihedrals_list.append([int(i) for i in tokens])
                        except (Exception, ):
                            m = "\t\t ERROR in {} file. Distances".format(args.indxfile)
                            print(m) if self._logger is None else self._logger.info(m)
                            exit()

        dict_2d_dih_angles = defaultdict(list)
        for key, values in self._optgeometry.items():
            for idist in distances_list:
                m = "\t\t ERROR in {} file. " \
                    "Distances are not yet implemented in gecos_outputgaussian.py".format(args.indxfile)
                print(m) if self._logger is None else self._logger.info(m)
                exit()
            for iangle in angles_list:
                m = "\t\t ERROR in {} file. " \
                    "Angles are not yet implemented in gecos_outputgaussian.py".format(args.indxfile)
                print(m) if self._logger is None else self._logger.info(m)
                exit()
            for idih in dihedrals_list:
                at1, at2, at3, at4 = idih[0:4]
                c1 = values[1][at1-1]
                c2 = values[1][at2-1]
                c3 = values[1][at3-1]
                c4 = values[1][at4-1]
                dih_angle = dihedral_py(c1, c2, c3, c4, units="degree")
                if at4 < at1:
                    label = "{}-{}-{}-{}".format(at4, at3, at2, at1)
                else:
                    label = "{}-{}-{}-{}".format(at1, at2, at3, at4)
                self._internalcoords[key][label] = dih_angle
                dict_2d_dih_angles[label].append(dih_angle)

        m = "\n\t\tInternal coordinates have been written in distances.dat, angles.dat and/or dihedral.dat files.\n"
        print(m) if self._logger is None else self._logger.info(m)

        de_energy_array = np.zeros(len(self._internalcoords), dtype=np.float32)
        idx = 0
        with open("dihedral.dat", 'w') as fdih:
            for ikey, ivalues in self._internalcoords.items():
                line = ikey
                for jkey, jvalues in ivalues.items():
                    line += " {0:.1f} ".format(jvalues)
                ecurr = cclib.parser.utils.convertor(self._optscfenergies[ikey], 'eV', 'hartree')
                de = ecurr - self._minscfenergy_ha
                de_kcalmol = cclib.parser.utils.convertor(de, 'hartree', 'kcal/mol')
                de_energy_array[idx] = de_kcalmol
                line += "{0:.2f} ".format(de_kcalmol)
                fdih.writelines(line+"\n")
                idx += 1

        # Output data in grid format for GNUPLOT
        try:
            res_grid = 1
            half_width = 0.5*res_grid
            npoints = int(360 / res_grid) + 1
            mask = np.linspace(-180, 180, npoints)
            self._2dgrid = np.full((npoints, npoints, 2), 999999, dtype=np.float32)

            pair_of_angles = np.full([len(self._internalcoords), 2], -1000, dtype=np.float32)
            icol = 0
            for key, values in dict_2d_dih_angles.items():
                pair_of_angles[:, icol] = np.asarray(values)
                icol += 1

            idx = 0
            for ipair in pair_of_angles:
                xreal = ipair[0]
                yreal = ipair[1]
                for ix in range(0, len(mask)):
                    if mask[ix] < xreal < mask[ix+1]:
                        if mask[ix]+half_width > xreal:
                            xgrid = int(mask[ix])
                        else:
                            xgrid = int(mask[ix+1])
                        break
                for iy in range(0, len(mask)):
                    if mask[iy] < yreal < mask[iy+1]:
                        if mask[iy]+half_width > yreal:
                            ygrid = int(mask[iy])
                        else:
                            ygrid = int(mask[iy+1])
                        break

                de_kcalmol = de_energy_array[idx]

                if np.abs(self._2dgrid[xgrid, ygrid, 0] - 999999.0) < 1e-08:
                    self._2dgrid[xgrid, ygrid, 0] = de_kcalmol
                    self._2dgrid[xgrid, ygrid, 1] = 1
                else:
                    self._2dgrid[xgrid, ygrid, 0] += de_kcalmol
                    self._2dgrid[xgrid, ygrid, 1] += 1
                print(xreal, yreal, xgrid, ygrid, de_kcalmol)
                idx += 1

            nx, ny, _ = self._2dgrid.shape
            line = ""
            for ix in range(0, nx):
                for iy in range(0, ny):
                    if np.abs(self._2dgrid[ix, iy, 0] - 999999.0) < 1e-08:
                        line += "{} {} {}\n".format(mask[ix], mask[iy], 1000.)
                    else:
                        line += "{} {} {}\n".format(mask[ix], mask[iy], self._2dgrid[ix, iy, 0] / self._2dgrid[ix, iy, 1])
                line += "\n"

            with open("dihedral_grid.dat", 'w') as fdihgrid:
                fdihgrid.writelines(line)
        except UnboundLocalError:
            pass
        except IndexError:
            pass

    # =========================================================================
    def write_to_log(self, logfolder, generate_data_gnuplot=True):

        df = self._df

        # Files
        m = "\t\tLocaldir                      : {}\n".format(os.getcwd())
        m += "\t\tDirectory with Log files      : {}\n".format(logfolder)
        print(m) if self._logger is None else self._logger.info(m)

        # Write table
        m = "\n\t\t{0:1s} {1:^48s} {2:^14s} {3:^19s} {4:^10s} {5:^17s} " \
            "{6:^16s} {7:^16s} {8:^16s} {9:^28s} {10:^8s} " \
            "{11:^17s} {12:^17s}\n".format('#', 'ID', 'Type', 'Rel. Energy(kcal/mol)',
                                           'RMSDwithHs(A)', 'RMSDwithoutHs(A)',
                                           'RotConstA(Ghz)', 'RotConstB(Ghz)',
                                           'RotConstC(Ghz)', 'AngleMainPrincipalAxes(deg)',
                                           'NHbonds', 'NvdwCont_Intermol', 'NvdwCont_Intramol')

        lenm = len(m)
        m += "\t\t# " + len(m) * "=" + "\n"
        for ind in df.index:
            line = "\t\t{0:6d} {1:^42s} {2:<14s} {3:>14.2f} {4:>20.3f} {5:>14.3f} {6:>17.6f}  " \
                   "{7:>14.6f}   {8:>14.6f}  {9:>20.1f}  {10:14d}  {11:14d}  {12:14d}\n" \
                .format(ind, df['ID'][ind],
                        str(df['Type'][ind]),
                        df['DeltaEnergy(kcal/mol)'][ind],
                        df['RMSDwithHs'][ind],
                        df['RMSDwithoutHs'][ind],
                        df['RotConstA(Ghz)'][ind],
                        df['RotConstB(Ghz)'][ind],
                        df['RotConstC(Ghz)'][ind],
                        df['AngleMainPrincipalAxes(deg)'][ind],
                        df['NHbonds'][ind],
                        df['NvdwContacts_Intermol'][ind],
                        df['NvdwContacts_Intramol'][ind], )
            m += line
            self._nameID_to_idx[df['ID'][ind]] = ind
            self._idx_to_nameID[ind] = df['ID'][ind]

        # Skip last carrier return
        m = m[0:-1]
        print(m) if self._logger is None else self._logger.info(m)

        # Job Time
        m1 = "\t\t# " + lenm * "=" + "\n"
        print(m1) if self._logger is None else self._logger.info(m1)
        end = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m1 = "\t\t Job: {}\n".format(end)
        print(m1) if self._logger is None else self._logger.info(m1)

        # Generate data to draw in gnuplot
        if generate_data_gnuplot:
            if self._logger is not None:
                folder, basename = os.path.split(self._logger.handlers[0].baseFilename)
                basename = os.path.splitext(basename)[0]
            else:
                basename = "gecos_energy_analysis.dat"
            with open(os.path.join(folder, basename + ".dat"), 'w') as fdata:
                lines = m.split("\n")
                for iline in lines:
                    if iline.count("Rotamer") != 0 or iline.count("Identical") != 0:
                        continue
                    fdata.writelines(iline+"\n")

            fname = os.path.join(folder, basename + ".dat")
            self._gnuplot_template(fname)

    # =========================================================================
    def write_vib_to_log(self, logfolder, generate_data_gnuplot=True):

        df = self._df

        # Files
        m = "\t\tLocaldir                      : {}\n".format(os.getcwd())
        m += "\t\tDirectory with Log files      : {}\n".format(logfolder)
        print(m) if self._logger is None else self._logger.info(m)

        # Write table
        m = "\n\t\t{0:1s} {1:^48s} {2:10s} {3:^19s} {4:^19s} {5:^17s} " \
            "{5:^16s} \n".format('#', 'ID', 'Type', 'DeltaG(kcal/mol)',
                                 'Delta_-TS(kcal/mol)', 'DeltaH(kcal/mol)',
                                 'DeltaEscf(kcal/mol)')

        lenm = len(m)
        m += "\t\t# " + len(m) * "=" + "\n"
        for ind in df.index:
            try:
                line = "\t\t{0:6d} {1:^40s} {2:<14s} {3:^18.2f} {4:>18.2f} {5:>18.2f} {6:>18.2f}\n" \
                    .format(ind, df['ID'][ind],
                            df['Type'][ind],
                            df['Delta_G(kcal/mol)'][ind],
                            df['Delta_-TS(kcal/mol)'][ind],
                            df['Delta_H(kcal/mol)'][ind],
                            df['Delta_E(kcal/mol)'][ind])
            except ValueError:
                line = "\t\t{0:6d} {1:^40s} {2:<12s} {3:^18.2f} {4:>18.2f} {5:>18.2f} {6:>18.2f}\n" \
                    .format(ind, df['ID'][ind],
                            'N/A',
                            df['Delta_G(kcal/mol)'][ind],
                            df['Delta_-TS(kcal/mol)'][ind],
                            df['Delta_H(kcal/mol)'][ind],
                            df['Delta_E(kcal/mol)'][ind])

            m += line
        m = m[:-1]   # Remove last \n
        print(m) if self._logger is None else self._logger.info(m)

        # Job Time
        m1 = "\t\t# " + lenm * "=" + "\n"
        print(m1) if self._logger is None else self._logger.info(m1)
        end = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m1 = "\t\t Job: {}\n".format(end)
        print(m1) if self._logger is None else self._logger.info(m1)

        # Generate data to draw in gnuplot
        if generate_data_gnuplot:
            if self._logger is not None:
                folder, basename = os.path.split(self._logger.handlers[0].baseFilename)
                basename = os.path.splitext(basename)[0]
            else:
                basename = "gecos_thermo_analysis.dat"
            with open(os.path.join(folder, basename + ".dat"), 'w') as fdata:
                lines = m.split("\n")
                for iline in lines:
                    if iline.count("Rotamer") != 0 or iline.count("Identical") != 0:
                        continue
                    fdata.writelines(iline+"\n")

    # =========================================================================
    def write_resp_to_log(self, logfolder, generate_data_gnuplot=True):

        df = self._df

        # Files
        m = "\t\tLocaldir                      : {}\n".format(os.getcwd())
        m += "\t\tDirectory with Log files      : {}\n".format(logfolder)
        print(m) if self._logger is None else self._logger.info(m)

        # Write table
        m = "\n\t\t{0:1s} {1:^30s} {2:^5s} {3:^19s} {4:^8s} {5:^17s} " \
            "{6:^16s} {7:^16s} {8:^16s} \n".format('#', 'ID', 'Type', 'DeltaEnergy(kcal/mol)',
                                                   'RMSDwithHs(A)', 'RMSDwithoutHs(A)',
                                                   'RotConstA(Ghz)', 'RotConstB(Ghz)',
                                                   'RotConstC(Ghz)')

        lenm = len(m)
        m += "\t\t# " + len(m) * "=" + "\n"
        for ind in df.index:
            line = "\t\t{0:^42s} {1:<9s} {2:>14.2f} {3:>20.3f} {4:>14.3f} {5:>17.6f}  " \
                   "{6:>14.6f}   {7:>14.6f}\n" \
                .format(df['ID'][ind],
                        str(df['Type'][ind]),
                        df['DeltaEnergy(kcal/mol)'][ind],
                        df['RMSDwithHs'][ind],
                        df['RMSDwithoutHs'][ind],
                        df['RotConstA(Ghz)'][ind],
                        df['RotConstB(Ghz)'][ind],
                        df['RotConstC(Ghz)'][ind])
            m += line
        print(m) if self._logger is None else self._logger.info(m)

        # Job Time
        m1 = "\t\t# " + lenm * "=" + "\n"
        print(m1) if self._logger is None else self._logger.info(m1)
        end = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m1 = "\t\t Job: {}\n".format(end)
        print(m1) if self._logger is None else self._logger.info(m1)

        # Generate data to draw in gnuplot
        if generate_data_gnuplot:
            if self._logger is not None:
                folder, basename = os.path.split(self._logger.handlers[0].baseFilename)
                basename = os.path.splitext(basename)[0]
            else:
                basename = "gecos_resp_prep_analysis.dat"
            with open(os.path.join(folder, basename + ".dat"), 'w') as fdata:
                fdata.writelines(m)

    # =========================================================================
    def write_clustering_to_log(self, logfolder, generate_data_gnuplot=True):

        df = self._df

        # Files
        m = "\t\tLocaldir                      : {}\n".format(os.getcwd())
        m += "\t\tDirectory with Log files      : {}\n".format(logfolder)
        print(m) if self._logger is None else self._logger.info(m)

        # Write table
        m = "\n\t\t{0:1s} {1:^30s} {2:^5s} {3:^19s} {4:^8s} {5:^17s} " \
            "{6:^16s} {7:^16s} {8:^16s} \n".format('#', 'ID', 'Type', 'DeltaEnergy(kcal/mol)',
                                                   'RMSDwithHs(A)', 'RMSDwithoutHs(A)',
                                                   'RotConstA(Ghz)', 'RotConstB(Ghz)',
                                                   'RotConstC(Ghz)')

        lenm = len(m)
        m += "\t\t# " + len(m) * "=" + "\n"
        for ind in df.index:
            line = "\t\t{0:^42s} {1:<9s} {2:>14.2f} {3:>20.3f} {4:>14.3f} {5:>17.6f}  " \
                   "{6:>14.6f}   {7:>14.6f}\n" \
                .format(df['ID'][ind],
                        str(df['Type'][ind]),
                        df['DeltaEnergy(kcal/mol)'][ind],
                        df['RMSDwithHs'][ind],
                        df['RMSDwithoutHs'][ind],
                        df['RotConstA(Ghz)'][ind],
                        df['RotConstB(Ghz)'][ind],
                        df['RotConstC(Ghz)'][ind])
            m += line
        print(m) if self._logger is None else self._logger.info(m)

        # Job Time
        m1 = "\t\t# " + lenm * "=" + "\n"
        print(m1) if self._logger is None else self._logger.info(m1)
        end = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m1 = "\t\t Job: {}\n".format(end)
        print(m1) if self._logger is None else self._logger.info(m1)

        # Generate data to draw in gnuplot
        if generate_data_gnuplot:
            if self._logger is not None:
                folder, basename = os.path.split(self._logger.handlers[0].baseFilename)
                basename = os.path.splitext(basename)[0]
            else:
                basename = "gecos_resp_prep_analysis.dat"
            with open(os.path.join(folder, basename + ".dat"), 'w') as fdata:
                fdata.writelines(m)

    # =========================================================================
    def close_contacts(self, args):

        if len(self._xyzlist) == 0:
            m1 = "\t\t No XYZ coordinates are found.!!!!!!\n"
            m2 = "\t\t Close Contacts cannot be calculated.!!!!!!"
            m3 = "\n\t\t" + len(m2) * "*" + "\n"
            print(m3 + m1 + m2 + m3) if self._logger is None else self._logger.warn(m3 + m1 + m2 + m3)
            return False

        # Initialize instance Contacts
        contacts = CloseContacts(self._logger)

        # For each molecule in the class search for close_contacts
        nhbonds_list = []
        nvdwcontacts_inter = []
        nvdwcontacts_intra = []
        for idx, ixyzmol in enumerate(self._xyzlist):
            molob = ob.OBMol()
            obmol = ob.OBConversion()
            obmol.SetInAndOutFormats('xyz', 'xyz')
            obmol.ReadString(molob, self._xyzlist[idx])

            # Number of fragments. Getting info from fragments
            numfrag = molob.Separate()
            if len(numfrag) != 2 and len(numfrag) != 1:
                m = "\t\tTwo fragments are expected ({} fragments found.\n)".format(len(numfrag))
                m += "\t\t\t Close contacts are not calculated.)".format(len(numfrag))
                print(m) if self._logger is None else self._logger.warning(m)
                self._df['NHbonds'] = len(self._xyzlist)*[0]
                self._df['NvdwContacts_Intermol'] = len(self._xyzlist)*[0]
                self._df['NvdwContacts_Intramol'] = len(self._xyzlist)*[0]
                return False

            # For each fragment get the necessary info in idatoms_info_frag dictionary
            idatoms_info_frag = defaultdict(list)
            list_nb_pairs = defaultdict(list)
            for idx_ifrag, ifrag_ob in enumerate(numfrag):
                for pair in ob.OBMolPairIter(ifrag_ob):
                    (first, second) = pair
                    begin = ifrag_ob.GetAtom(first).GetId()
                    end = ifrag_ob.GetAtom(second).GetId()
                    list_nb_pairs[idx_ifrag].append({begin, end})
                for iat in ob.OBMolAtomIter(ifrag_ob):
                    symbol = ob.GetSymbol(iat.GetAtomicNum())
                    idatoms_info_frag[idx_ifrag].append([iat.GetId(),
                                                         symbol,
                                                         [iat.GetX(), iat.GetY(), iat.GetZ()]])
            if hasattr(args, "vdw") and args.vdw is not None:
                contacts.vanderwaals(self._df['ID'][idx],
                                     idatoms_info_frag,
                                     list_nb_pairs,
                                     delta=args.vdw,
                                     noignore_hs=args.noignoreh)

            if hasattr(args, "hbonds") and args.hbonds is not None:

                # HD means H atom bonded to Donor atom (D), A=acceptor
                hd_atom_idx_list = []
                acc_atom_idx_list = []
                hd_d_idx_dict = {}
                coords = []
                for iatom in ob.OBMolAtomIter(molob):
                    coords.append([iatom.GetX(), iatom.GetY(), iatom.GetZ()])
                coords = np.array(coords)

                for iatom in ob.OBMolAtomIter(molob):
                    idx_atom = iatom.GetId()
                    if iatom.GetAtomicNum() == 1:
                        for neighbor in ob.OBAtomAtomIter(iatom):
                            if neighbor.IsHbondDonor():
                                hd_atom_idx_list.append(idx_atom)  # list of H atoms
                                hd_d_idx_dict[idx_atom] = neighbor.GetId()  # dictionary of H
                    if iatom.IsHbondAcceptor():
                        acc_atom_idx_list.append(idx_atom)

                contacts.hbonds(self._df['ID'][idx],
                                idatoms_info_frag,
                                coords,
                                hd_atom_idx_list,
                                hd_d_idx_dict,
                                acc_atom_idx_list,
                                args.hbonds[0],
                                args.hbonds[1])

            nhbonds_list.append(len(contacts._hbond_list[self._df['ID'][idx]]))
            nvdwcontacts_inter.append(len(contacts._vdwcontacts_intermol[self._df['ID'][idx]]))
            nvdwcontacts_intra.append(len(contacts._vdwcontacts_intramol[self._df['ID'][idx]]))

        self._df['NHbonds'] = nhbonds_list
        self._df['NvdwContacts_Intermol'] = nvdwcontacts_inter
        self._df['NvdwContacts_Intramol'] = nvdwcontacts_intra

        if hasattr(args, "hbonds") and args.hbonds is not None:
            with open("hydrogen_bonds_info.dat", 'w') as fhbond:
                line = "#        Name        Donor     H     Acceptor    Distance(A)    Angle(degree)\n"
                fhbond.writelines(line)
                for idmol in self._df["ID"]:
                    line = "#{0:30s}\n".format(idmol)
                    fhbond.writelines(line)
                    for ihbond in contacts._hbond_list[idmol]:
                        line = "{0:10d} {1:10d} {2:10d} {3:10.3f} {4:10.1f}\n".format(ihbond[0], ihbond[1],
                                                                                      ihbond[2], ihbond[3], ihbond[4])
                        fhbond.writelines(line)

        return True

    # =========================================================================
    def moment_of_inertia(self):

        """
        Calculate the angle between the principal axis of the moment of inertia

        :return:
        """

        # For each molecule in the class search for close_contacts
        localdir = os.getcwd()
        angle_list = []
        for idxmol, ixyzmol in enumerate(self._xyzlist):
            nmols_array = defaultdict(list)
            molob = ob.OBMol()
            obmol = ob.OBConversion()
            obmol.SetInAndOutFormats('xyz', 'pdb')
            obmol.ReadString(molob, self._xyzlist[idxmol])
            namefileref = os.path.join(localdir, "tmp.pdb")
            obmol.WriteFile(molob, namefileref)

            universe = MDAnalysis.Universe("tmp.pdb")
            vects = defaultdict()
            for idxfrg, frg in enumerate(universe.atoms.fragments):
                nmols_array[idxfrg] = list(frg.indices)
                seltxt = "index "
                for iat in nmols_array[idxfrg]:
                    seltxt += "{0:d} ".format(iat)
                g = universe.select_atoms(seltxt)
                utransp = g.principal_axes()
                u = utransp.T
                vects[idxfrg] = u[:, 0]  # main axis of the principal axes

            if len(vects) == 2:
                angle_list.append(np.arccos(np.dot(vects[0], vects[1])) * 180. / np.pi)
        try:
            self._df["AngleMainPrincipalAxes(deg)"] = angle_list
        except ValueError:
            pass

    # =========================================================================
    def radius_of_gyration(self):

        """
        Calculate the radius of gyration of the molecule.

        :return:
        """

        # For each molecule in the class search for close_contacts
        localdir = os.getcwd()
        rg_list = []
        for idxmol, ixyzmol in enumerate(self._xyzlist):
            nmols_array = defaultdict(list)
            molob = ob.OBMol()
            obmol = ob.OBConversion()
            obmol.SetInAndOutFormats('xyz', 'pdb')
            obmol.ReadString(molob, self._xyzlist[idxmol])
            namefileref = os.path.join(localdir, "tmp.pdb")
            obmol.WriteFile(molob, namefileref)

            universe = MDAnalysis.Universe("tmp.pdb")

            coordinates = universe.atoms.positions
            center_of_mass = universe.atoms.center_of_mass()
            masses = universe.atoms.masses
            total_mass = sum(masses)

            # get squared distance from center
            ri_sq = (coordinates-center_of_mass)**2
            # sum the unweighted positions
            sq = np.sum(ri_sq, axis=1)
            sq_x = np.sum(ri_sq[:, [1, 2]], axis=1)  # sum over y and z
            sq_y = np.sum(ri_sq[:, [0, 2]], axis=1)  # sum over x and z
            sq_z = np.sum(ri_sq[:, [0, 1]], axis=1)  # sum over x and y

            # make into array
            sq_rs = np.array([sq, sq_x, sq_y, sq_z])

            # weight positions
            rog_sq = np.sum(masses*sq_rs, axis=1)/total_mass
            rg_list.append(np.sqrt(rog_sq)[0])

        try:
            self._df['RadiusGyration(ang)'] = rg_list
        except ValueError:
            pass

        with open("radius_gyration.dat", 'w') as frg:
            line = "#ID Name Delta_Energy(kcal/mol) Rg(angstroms)\n"
            frg.writelines(line)
            for ind in self._df.index:
                line = "\t\t{0:6d} {1:^42s} {2:>14.2f} {3:>14.3f} \n" \
                    .format(ind, self._df['ID'][ind],
                            self._df['DeltaEnergy(kcal/mol)'][ind],
                            self._df['RadiusGyration(ang)'][ind])
                frg.writelines(line)

    # =========================================================================
    def cluster_conformers(self, rmsd_only_heavy=True, energy_thr=0.1, rot_constant_thr=0.0005,
                           rmsd_thr=1.0, window_energy=1000.0):

        hartrees_to_kcalmol = 627.509391

        # Get heavy atoms indexes
        list_indices_rmsd_atoms = []
        all_log_atoms = []
        for key, values in self._optgeometry.items():
            natoms = len(values[0])
            all_log_atoms.append(natoms)

        # Get heavy atoms indexes
        list_indices_rmsd_atoms = []
        all_log_atoms = []
        for key, values in self._optgeometry.items():
            natoms = len(values[0])
            atomic_numbers_log = values[0]
            all_log_atoms.append(natoms)

        if len(set(all_log_atoms)) != 1:
            now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
            m = "\t\t\t Logs must have the same number ot atoms to clusterize ({})".format(now)
            print(m) if self._logger is None else self._logger.error(m)
            exit()

        for iatom in range(0, natoms):
            if rmsd_only_heavy:
                if atomic_numbers_log[iatom] != 1:
                    list_indices_rmsd_atoms.append(iatom)
            elif not rmsd_only_heavy:
                list_indices_rmsd_atoms.append(iatom)

        idx = 0
        cluster = defaultdict(dict)
        icluster = 0
        n_conformers = self._df.shape[0]
        self._cre_dict['Conformers'] = 0
        self._cre_dict['Rotamers'] = 0
        self._cre_dict['Identical'] = 0
        self._cre_dict['Others'] = 0

        dict_conftype = defaultdict()
        for index in self._df.index:
            idx += 1
            try:
                if index % int(n_conformers * 0.10) == 0:
                    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
                    m = "\t\t\t Clustering Ensemble {} of {} ({})".format(idx, n_conformers, now)
                    print(m) if self._logger is None else self._logger.info(m)
            except ZeroDivisionError:
                now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
                m = "\t\t\t Clustering Ensemble {} of {} ({})".format(idx, n_conformers, now)
                print(m) if self._logger is None else self._logger.info(m)

            energy = self._df['OptEnergy(Ha)'][index] * hartrees_to_kcalmol
            rot_constant = [self._df['RotConstA(Ghz)'][index],
                            self._df['RotConstB(Ghz)'][index],
                            self._df['RotConstC(Ghz)'][index]]

            if icluster == 0:
                icluster += 1
                cluster[icluster] = {"seed": index, "lowest_energy": energy, "highest_energy": energy,
                                     "nelements": 0, "elements": [], "pairs": [], "files": [], "rot_constant": []}
                cluster[icluster]["pairs"].append([0.000, energy])
                cluster[icluster]["files"].append("")
                cluster[icluster]["nelements"] += 1
                cluster[icluster]["elements"].append(index)
                cluster[icluster]["rot_constant"] = rot_constant
                self._df.at[index, 'IDCluster'] = icluster
                self._df.at[index, 'Type'] = "Conformer"
                iconf_min_energy = index
                self._cre_dict['Conformers'] = 1
                dict_conftype[self._df.at[index, 'ID']] = "Conformer"
            else:
                found = False
                iconf_target = index
                for idx_cluster in range(1, icluster + 1):
                    iconf_ref = cluster[idx_cluster]["seed"]
                    iener_ref = cluster[idx_cluster]["lowest_energy"]
                    irot_constant_ref = cluster[idx_cluster]["rot_constant"]

                    # Check rotamers in accordance of https://doi.org/10.1039/C9CP06869D
                    #   Conditions:
                    #      1. DE > E_thr AND DRMSD > RMSD_thr AND DB > B_thr : Conformers
                    #      2. DE < E_thr AND DRMSD > RMSD_thr AND DB > B_thr : Rotamers
                    #      3. DE < E_thr AND DRMSD < RMSD_thr AND DB < B_thr : Identical
                    delta_energy = np.abs(energy - iener_ref)

                    # Remove conformers with high energy difference from the lowest one
                    if delta_energy > window_energy:
                        found = True
                        break

                    if delta_energy > energy_thr:
                        pass
                    else:
                        rmsd_noh = self.getconformeropenbabelrms(iconf_ref, iconf_target,
                                                                 atomids=list_indices_rmsd_atoms, align=True)
                        if rmsd_noh > rmsd_thr:
                            pass
                        else:
                            db = [np.abs(irot_constant_ref[i] - rot_constant[i]) for i in range(0, 3)]
                            if all(i < rot_constant_thr for i in db):
                                self._cre_dict['Identical'] += 1
                                self._df.at[index, 'Type'] = "Identical_"+str(iconf_ref)
                                dict_conftype[self._df.at[index, 'ID']] = "Identical_"+str(iconf_ref)
                            else:
                                self._cre_dict['Rotamers'] += 1
                                self._df.at[index, 'Type'] = "Rotamer_"+str(iconf_ref)
                                dict_conftype[self._df.at[index, 'ID']] = "Rotamer_" + str(iconf_ref)
                            cluster[idx_cluster]["pairs"].append([rmsd_noh, energy])
                            cluster[idx_cluster]["files"].append("")
                            cluster[idx_cluster]["highest_energy"] = energy
                            cluster[idx_cluster]["nelements"] += 1
                            cluster[idx_cluster]["elements"].append(index)
                            self._df.at[index, 'IDCluster'] = idx_cluster
                            found = True
                            break

                if not found:
                    icluster += 1
                    self._cre_dict['Conformers'] += 1
                    self._df.at[index, 'Type'] = "Conformer"
                    dict_conftype[self._df.at[index, 'ID']] = "Conformer"
                    cluster[icluster] = {"seed": index, "lowest_energy": energy, "highest_energy": energy,
                                         "nelements": 0, "elements": [], "pairs": [], "files": [], "rot_constant": []}
                    cluster[icluster]["files"].append(index)
                    cluster[icluster]["nelements"] += 1
                    cluster[icluster]["elements"].append(index)
                    rmsd_noh = self.getconformeropenbabelrms(iconf_min_energy, iconf_target,
                                                             atomids=list_indices_rmsd_atoms, align=True)
                    cluster[icluster]["pairs"].append([rmsd_noh, energy])
                    cluster[icluster]["rot_constant"] = rot_constant
                    self._df.at[index, 'IDCluster'] = icluster

        # Write Conformer/Rotamers/Identical info
        m = "\n"
        s = 0
        for key, values in self._cre_dict.items():
            m += "\t\t\t {0:<20s}: {1:d}\n".format(key, values)
            s += values
        key = "Total"
        m += "\t\t\t {0:<20s}: {1:d}\n".format(key, s)
        print(m) if self._logger is None else self._logger.info(m)

        m = "\t\tThresholds used to classify conformers are:\n" \
            "\t\t energy_thr={} kcal/mol, rot_constant_thr={} GHz,\n" \
            "\t\t rmsd_thr={} Angs window_energy={} kcal/mol rms_only_heavy={}\n".\
            format(energy_thr, rot_constant_thr, rmsd_thr, window_energy, rmsd_only_heavy)
        print(m) if self._logger is None else self._logger.info(m)

        return cluster, dict_conftype

    # =========================================================================
    def getconformeropenbabelrms(self, iconf, jconf, atomids=None, align=True):

        """
        Calculate the rmsd between all molecules and the lowest energy molecule
        as reference

        """

        rmsd_value = 0.0

        igeomxyz = self._xyzlist[iconf]
        jgeomxyz = self._xyzlist[jconf]

        molref = ob.OBMol()
        obref = ob.OBConversion()
        obref.SetInAndOutFormats('xyz', 'xyz')
        obref.ReadString(molref, igeomxyz)
        namefileref = os.path.join("./", "reference_allign.xyz")
        obref.WriteFile(molref, namefileref)
        _, coord_ref = rmsd.get_coordinates_xyz(namefileref)
        coord_ref = copy.deepcopy(coord_ref[atomids])
        os.remove(namefileref)

        obtarget = ob.OBConversion()
        obtarget.SetInAndOutFormats('xyz', 'xyz')
        moltarget = ob.OBMol()
        obtarget.ReadString(moltarget, jgeomxyz)
        a = ob.OBAlign(False, False)
        a.SetRefMol(molref)
        a.SetTargetMol(moltarget)
        a.Align()
        a.UpdateCoords(moltarget)
        namefile = os.path.join("./", "target_allign.xyz")
        obref.WriteFile(moltarget, namefile)
        _, coord_target = rmsd.get_coordinates_xyz(namefile)
        coord_target = copy.deepcopy(coord_target[atomids])
        rmsd_value = rmsd.rmsd(coord_target, coord_ref)
        os.remove(namefile)

        return rmsd_value

    # =========================================================================
    def _gnuplot_template(self, fname, nx=3, ny=2):

        try:
            nconf = self._cre_dict['Conformers']
        except KeyError:
            nconf = len(self._idx_to_nameID)
        startconf = -1
        endconf = nconf + 1
        ticx_conf = int(np.floor((nconf + 1)/10))
        if ticx_conf == 0:
            ticx_conf = 1
        max_nhbonds = max(self._df.NHbonds)
        min_nhbonds = min(self._df.NHbonds)
        max_intervdw = max(self._df.NvdwContacts_Intermol)
        min_intervdw = min(self._df.NvdwContacts_Intermol)
        max_intravdw = max(self._df.NvdwContacts_Intramol)
        min_intravdw = min(self._df.NvdwContacts_Intramol)


        linegnuplot = "reset\n"
        linegnuplot += 'f1="{}"\n'.format(fname)
        linegnuplot += "\n"
        linegnuplot += "# ============== COMMON FORMATS ====================\n"
        linegnuplot += 'set style line 1 lt 1 ps 1.0 lc rgb "black"  pt 6 lw 2.0\n'
        linegnuplot += 'set style line 2 lt 1 ps 1.0 lc rgb "black"  pt 4 lw 2.0\n'
        linegnuplot += 'set style line 3 lt 1 ps 1.0 lc rgb "black"  pt 6 lw 2.0\n'
        linegnuplot += "\n"
        linegnuplot += "# ============== 2D PLOT =====================\n"
        linegnuplot += 'set term wxt 1 enhanced dashed size 900,800 font "Arial,12"\n'
        linegnuplot += "set multiplot layout {},{}\n".format(nx, ny)
        linegnuplot += "\n"

        linegnuplot += 'set title "Relative energy of conformers"\n'
        linegnuplot += 'set xlabel "Conformer ID"\n'
        linegnuplot += 'set ylabel "Rel. Energy (kcal/mol)\n'
        linegnuplot += "set xtics {}\n".format(ticx_conf)
        linegnuplot += "set xrange[{}:{}]\n".format(startconf, endconf)
        linegnuplot += 'set format y "%.1f"\n'
        linegnuplot += "set grid\n"
        linegnuplot += "p f1 u 1:4 ls 1 notitle\n"
        linegnuplot += "unset xrange\n"
        linegnuplot += "\n"

        linegnuplot += 'set title "Relative energy vs RMSD"\n'
        linegnuplot += 'set xlabel "Rel. Energy (kcal/mol)\n'
        linegnuplot += 'set ylabel "RMSD (Angstroms)"\n'
        linegnuplot += 'set format x "%.1f"\n'
        linegnuplot += 'set format y "%.1f"\n'
        linegnuplot += "set xtics auto\n"
        linegnuplot += "p f1 u 4:6 ls 2 notitle\n"
        linegnuplot += "\n"

        linegnuplot += 'set title "Relative energy vs Number H bonds"\n'
        linegnuplot += 'set xlabel "Rel. Energy (kcal/mol)"\n'
        linegnuplot += 'set ylabel "# H-bonds"\n'
        linegnuplot += 'set format x "%.1f"\n'
        linegnuplot += 'unset format y\n'
        linegnuplot += 'set ytics 1\n'
        linegnuplot += "set xtics auto\n"
        linegnuplot += "set yrange [{}:{}]\n".format(min_nhbonds-1, max_nhbonds+1)
        linegnuplot += "p f1 u 4:11 ls 3 notitle\n"
        linegnuplot += "unset yrange\n"
        linegnuplot += "\n"

        linegnuplot += 'set title "Relative energy vs Intermolecular Contacts"\n'
        linegnuplot += 'set xlabel "Rel. Energy (kcal/mol)"\n'
        linegnuplot += 'set ylabel "# Intermolecular Contacts"\n'
        linegnuplot += 'set format x "%.1f"\n'
        linegnuplot += 'unset format y\n'
        linegnuplot += 'set ytics 1\n'
        linegnuplot += "set xtics auto\n"
        linegnuplot += "set yrange [{}:{}]\n".format(min_intervdw - 1, max_intervdw + 1)
        linegnuplot += "p f1 u 4:12 ls 3 notitle\n"
        linegnuplot += "unset yrange\n"
        linegnuplot += "\n"

        linegnuplot += 'set title "Relative energy vs Intramolecular Contacts"\n'
        linegnuplot += 'set xlabel "Rel. Energy (kcal/mol)"\n'
        linegnuplot += 'set ylabel "# Intramolecular Contacts"\n'
        linegnuplot += 'set format x "%.1f"\n'
        linegnuplot += 'set format y "%.0f"\n'
        linegnuplot += 'unset format y\n'
        linegnuplot += 'set ytics 1\n'
        linegnuplot += "set xtics auto\n"
        linegnuplot += "set yrange [{}:{}]\n".format(min_intravdw - 1, max_intravdw + 1)
        linegnuplot += "p f1 u 4:13 ls 3 notitle\n"
        linegnuplot += "unset yrange\n"
        linegnuplot += "\n"

        linegnuplot += 'set title "Number H bonds vs RMDS (with energy)"\n'
        linegnuplot += 'set xlabel "RMSD (Angstroms)"\n'
        linegnuplot += 'set ylabel "# H-bonds"\n'
        linegnuplot += 'set format x "%.1f"\n'
        linegnuplot += 'set format y "%.0f"\n'
        linegnuplot += "set xtics auto\n"
        linegnuplot += "set yrange [{}:{}]\n".format(min_nhbonds-1, max_nhbonds+1)
        linegnuplot += "set grid\n"
        linegnuplot += "set palette rgb 33,13,10\n"
        linegnuplot += "plot f1 u 6:11:4 palette pt 7 notitle\n"
        linegnuplot += "unset yrange\n"
        linegnuplot += "\n"

        linegnuplot += "unset multiplot\n"
        linegnuplot += "\n"
        linegnuplot += '# ============== 3D-POINT PLOT ===============\n'
        linegnuplot += 'set term wxt 2 enhanced dashed size 500,400 font "Arial,12"\n'
        linegnuplot += 'set multiplot layout 1,1\n'
        linegnuplot += '\n'
        linegnuplot += 'unset view\n'
        linegnuplot += 'set zlabel "Rel. Energy (kcal/mol)\n'
        linegnuplot += 'set xlabel "RMSD (Angstroms)"\n'
        linegnuplot += 'set ylabel "# H-bonds"\n'
        linegnuplot += "set xtics auto\n"
        linegnuplot += "set ytics auto\n"
        linegnuplot += "set ztics auto\n"
        linegnuplot += 'splot f1 u 6:11:4 ls 1 notitle w p\n'
        linegnuplot += 'unset multiplot\n'

        if os.path.exists("radius_gyration.dat"):
            linegnuplot += "\n"
            linegnuplot += '# ============== Radius of gyration ===============\n'
            linegnuplot += 'f2="./radius_gyration.dat"\n'
            linegnuplot += 'set term wxt 3 enhanced dashed size 500,400 font "Arial,12"\n'
            linegnuplot += 'set multiplot layout 1,1\n'
            linegnuplot += 'set title "Relative energy vs Radius of gyration"\n'
            linegnuplot += 'set xlabel "Rel. Energy (kcal/mol)"\n'
            linegnuplot += 'set ylabel "# Rg (Angstoms)"\n'
            linegnuplot += 'set format x "%.1f"\n'
            linegnuplot += 'set format y "%.1f"\n'
            linegnuplot += 'unset format y\n'
            linegnuplot += 'set ytics 1\n'
            linegnuplot += "set xtics auto\n"
            linegnuplot += "unset yrange\n"
            linegnuplot += "p f2 u 3:4 ls 3 notitle\n"
            linegnuplot += "unset yrange\n"
            linegnuplot += 'unset multiplot\n'
            linegnuplot += "\n"

        with open("template_gnuplot_energy_analysis.gnu", "w") as fgnuplot:
            fgnuplot.writelines(linegnuplot)
