import argparse
import numpy as np

import utils
import os
import datetime
try:
    from gecos_analysis.gecos_outputgaussian import GaussianGecos
    from gecos_analysis.gecos_irspectra import GecosIRSpectra
    from gecos_analysis.gecos_exp_irspectra import GecosExpIRSpectra
    from gecos_analysis.gecos_resp_analysis import GecosRespAnalysis
    from gecos_analysis.gecos_CnnClusteringWrapper import CnnClusteringWrapper
    from gecos_analysis.gecos_compare_spectrum import GecosCompareSpectrum
except ModuleNotFoundError:
    from gecos_outputgaussian import GaussianGecos
    from gecos_irspectra import GecosIRSpectra
    from gecos_exp_irspectra import GecosExpIRSpectra
    from gecos_resp_analysis import GecosRespAnalysis
    from gecos_CnnClusteringWrapper import CnnClusteringWrapper
    from gecos_compare_spectrum import GecosCompareSpectrum


class SmartFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)


# =============================================================================
def parse_arguments():

    desc = """Energy QM analysis.
    This is part of the gecos library"""

    parser = argparse.ArgumentParser(description=desc, formatter_class=SmartFormatter)

    # Subparsers
    subparser = parser.add_subparsers(dest='command', required=True)
    energy = subparser.add_parser('energy', formatter_class=SmartFormatter)
    spectrum = subparser.add_parser('spectrum', formatter_class=SmartFormatter)
    exp_spectrum = subparser.add_parser('exp_spectrum', formatter_class=SmartFormatter)
    resp_prep = subparser.add_parser('resp_prep', formatter_class=SmartFormatter)
    clustering = subparser.add_parser("clustering", formatter_class=SmartFormatter)
    similarity = subparser.add_parser("similarity", formatter_class=SmartFormatter)

    # Arguments for energy command
    energy.add_argument("-d", "--logfolder", dest="logfolder",
                        help="Folder containing the QM outputs.\n ",
                        action="store", required=True, default=None)

    energy.add_argument("-p", "--package", dest="qmpackage", choices=['g16'],
                        help="QM package to analyse (Default: Gaussian16).\n ",
                        action="store", required=False, default="g16")

    energy.add_argument("--log", dest="log",
                        help="Name of the file to write logs from this command",
                        action="store", required=False, default="gecos_energy_analysis.log")

    energy.add_argument("--vdw", dest="vdw", nargs="?",
                        help="Calculate close contacts of type Van der Waals\n"
                             "with delta (default = 0.25A) as:\n"
                             " d(at1-at2) < r_vdw(at1)+r_vdw(at2)+delta",
                        action="store", required=False, const=0.25, type=float)

    energy.add_argument("--hbonds", dest="hbonds", nargs='*',
                        help="Calculate hbonds in the molecule\n"
                             "using a distance and an angle.\n"
                             "Default arguments: 3.0 A and 140º\n",
                        metavar="hbond_dist_angstroms, hbond_angle_degrees",
                        action="store", required=False)

    energy.add_argument("--noignoreh", dest="noignoreh",
                        help="No ignore hydrogens for VdW Contacts",
                        action="store_true", required=False)

    energy.add_argument("--indxfile", dest="indxfile",
                        help="A file using the ndx format from "
                             "GROMACS to extract internal coordinates from the log"
                             "files.",
                        action="store", required=False, default=None)

    energy.add_argument("--clusterconf", dest="clusterconf",
                        help="Clusterize conformers",
                        action="store_true", required=False)

    energy.add_argument("--clusterparam", dest="clusterparam", nargs=5,
                        help="R|Parameters:\n"
                             " - Energy threshold (kcal/mol) \n"
                             " - Difference in Rotation constants (GHz)\n"
                             " - RMSD threshold (Angstroms).\n"
                             " - Window energy (kcal/mol)\n"
                             " - Only heavy atoms to calculate RMSD (True/False)\n"
                             "For more info: DOI	https://doi.org/10.1039/C9CP06869D",
                        action="store", required=False, default=[0.1, 0.0005, 1.0, 1000.0, True],
                        metavar=("Energy_Thresh", "RotConst_Thresh", "RMSD_Thresh",
                                 "Window_energy", "Only_Backbone_RMSD"))

    energy.add_argument("--radius_gyration", dest="rg",
                        help="Calculate the radius of gyration",
                        action="store_true", required=False)

    # Arguments for spectrum command
    spectrum.add_argument("-d", "--logfolder", dest="logfolder",
                          help="Folder containing the QM outputs.\n ",
                          action="store", required=True, default=None)

    spectrum.add_argument("-p", "--package", dest="qmpackage", choices=['g16'],
                          help="QM package to analyse (Default: Gaussian16).\n ",
                          action="store", required=False, default="g16")

    spectrum.add_argument("--log", dest="log",
                          help="Name of the file to write logs from this command",
                          action="store", required=False, default="gecos_spectrum_analysis.log")

    spectrum.add_argument("--scale", dest="scale",
                          help="Scale factor to apply to frequencies",
                          action="store", required=False, default=None, type=float)

    spectrum.add_argument("--start", dest="start",
                          help="Start frequency to fit the spectrum (cm^-1)",
                          action="store", required=False, default=0.0, type=float)

    spectrum.add_argument("--end", dest="end",
                          help="End frequency to fit the spectrum (cm^-1)",
                          action="store", required=False, default=4000.0, type=float)

    spectrum.add_argument("--npoints", dest="npoints",
                          help="Number of points to calculate the spectrum",
                          action="store", required=False, default=500, type=int)

    spectrum.add_argument("--width", dest="width",
                          help="Width used in the fit of the spectrum (cm^-1).",
                          action="store", required=False, default=10.0, type=float)

    spectrum.add_argument("--function", dest="function",
                          help="Function to plot the IR peaks", choices=["lorentzian", "lorentzian_gaussian"],
                          action="store", required=False, default="lorentzian")

    spectrum.add_argument("--avg", dest="averaged",
                          help="Calculate average spestrum",
                          choices=["boltzmann", "simple", "all"],
                          action="store", required=False, default=None)

    spectrum.add_argument("--clusterconf", dest="clusterconf",
                          help="Clusterize conformers",
                          action="store_true", required=False)

    spectrum.add_argument("--clusterparam", dest="clusterparam", nargs=5,
                          help="R|Parameters:\n"
                               " - Energy threshold (kcal/mol) \n"
                               " - Difference in Rotation constants (GHz)\n"
                               " - RMSD threshold (Angstroms).\n"
                               " - Window energy (kcal/mol)\n"
                               " - Only heavy atoms to calculate RMSD (True/False)\n"
                               "For more info: DOI	https://doi.org/10.1039/C9CP06869D",
                          action="store", required=False, default=[0.1, 0.0005, 1.0, 1000.0, True],
                          metavar=("Energy_Thresh", "RotConst_Thresh", "RMSD_Thresh",
                                   "Window_energy", "Only_Backbone_RMSD"))

    spectrum.add_argument("--similarity", dest="similarity",
                          help="R|The similarity of the calculated spectra with\n "
                               "the reference spectra (lowest energy structure)\n"
                               "using the Pearson’s product moment correlation coefficient\n"
                               "and the Spearman's rank correlation coefficient"
                               "For more info: DOI	https://dx.doi.org/10.1021/acs.jctc.0c00126",
                          action="store_true", required=False)

    spectrum.add_argument("-n", "--normalize", dest="normalize", choices=['minmax', 'vectornorm', '1-norm', 'snv'],
                          help="R|Method to normalize the spectra. Methods are explained in:\n"
                               "  \"A preliminary study on the importance of normalization method\n"
                               "  in Infrared MicroSpectroscopy for biomedical applications\".\n"
                               "  Andrea Zancla et al.\n"
                               "  24th IMEKO TC4 International Symposium\n"
                               "  22nd International Workshop on ADC and DAC Modelling and Testing\n"
                               "  IMEKO TC-4 2020\n"
                               "  Palermo, Italy, September 14-16, 2020\n",
                          action="store", required=False, default=None)

    # Arguments for preprocesing experimental spectrum command.
    exp_spectrum.add_argument("-f", "--file", dest="fileexpspectrum", type=str,
                              help="Name of the file containing the experimental spectrum",
                              action="store", required=True)

    exp_spectrum.add_argument("-r", "--resample", dest="resample_list", nargs=3, type=float,
                              help="R|Resample data from start to end with npoints.\n"
                                   "Paramters:\n"
                                   "  - start: First frecuency (cm^-1)\n"
                                   "  - end: Last frecuency (cm^-1)\n"
                                   "  - npoints: Number of points to resanmple",
                              action="store", metavar=["START", "END", "NPOINTS"], required=False, default=None)

    exp_spectrum.add_argument("--log", dest="log", type=str,
                              help="Name of the file to write logs from this command",
                              action="store", required=False, default="gecos_experimental_spectrum.log")

    exp_spectrum.add_argument("-n", "--normalize", dest="normalize", choices=['minmax', 'vectornorm', '1-norm', 'snv'],
                              help="R|Method to normalize the spectra. Methods are explained in:\n"
                                   "  \"A preliminary study on the importance of normalization method\n"
                                   "  in Infrared MicroSpectroscopy for biomedical applications\".\n"
                                   "  Andrea Zancla et al.\n"
                                   "  24th IMEKO TC4 International Symposium\n"
                                   "  22nd International Workshop on ADC and DAC Modelling and Testing\n"
                                   "  IMEKO TC-4 2020\n"
                                   "  Palermo, Italy, September 14-16, 2020\n",
                              action="store", required=False, default=None)

    exp_spectrum.add_argument("-b", "--baseline", dest="baseline",
                              choices=['ASLS', 'Mixture_Model', 'Optimize_extended_range', '2021Munich'],
                              help="R|Baseline correction of the spectrum.\n"
                                    " - ASLS (Asymmetric Least Squares) as implemented in pybaselines library " 
                                    "(Analytical Chemistry, 2003, 75(14), 3631-3636) \n"
                                    " - Mixture_Model as implemented in pybaselines library " 
                                    "(Chemometric and Intelligent Laboratory Systems, 2012, 117, 56-60) \n"
                                    " - Optimize_extended_range as implemented in pybaselines library " 
                                    "(Journal of Raman Spectroscopy, 2012, 43(12), 1884-1894) \n"
                                    " - 2021Muchich (J.Chem.TheoryComp.2021,17,985)",
                                    action="store",
                              required=False, default=None)

    exp_spectrum.add_argument("--baselineparams", dest="baselineparams", nargs='+',
                              help="R|Baseline correction parameters:\n"
                                   " - ASLS. Smoothing parameter (lambda) and penalizing weigth (p)"
                                   " factor. Recommended values: 1e6, 1e-2\n"
                                   "   https://pybaselines.readthedocs.io/en/latest/api/pybaselines"
                                   "/api/index.html#pybaselines.api.Baseline.asls\n"
                                   " - Mixture_Model. Smoothing parameter (lambda) and penalizing weigth (p)"
                                   " factor. Recommended values: 1e5, 1e-2\n"
                                   "   https://pybaselines.readthedocs.io/en/latest/api/pybaselines"
                                   "/api/index.html#pybaselines.api.Baseline.mixture_model\n"
                                   " - Adaptative_Min_Max. No parameters are required.\n"
                                   "   https://pybaselines.readthedocs.io/en/latest/api/pybaselines"
                                   "/api/index.html#pybaselines.api.Baseline.optimize_extended_range\n"
                                   " - Munich2021 method. Cutoff value of the diff. Recommended value: 0.006\n"
                                   "   https://pubs.acs.org/doi/10.1021/acs.jctc.0c01279?ref=pdf\n",
                              action="store",
                              required=False, default=None)

    exp_spectrum.add_argument("-s", "--smoothpeak", dest="smoothpeak", choices=['SG'],
                              help="R|Smooth peaks.\n"
                                    " - SG: Savitzky-Golay filter (Analytical Chemistry 36. p. 1627-1639)",
                              action="store",
                              required=False, default=None)

    exp_spectrum.add_argument("--smoothSGparam", dest="smoothSGparam", nargs=2, type=int,
                              help="R|Parameters:\n"
                                   " - Window_length \n"
                                   " - Polyorder\n" 
                                   "   (polyorder must be less than window_length).\n"
                                   "  See parameters for the Savitzky-Golay method:\n"
                                   "     https://docs.scipy.org/doc/scipy/reference/"
                                   "generated/scipy.signal.savgol_filter.html",
                              action="store",
                              required=False, default=None)

    # # Arguments for resp_prep command
    # resp_prep.add_argument("--")
    resp_prep.add_argument("--log", dest="log", type=str,
                           help="Name of the file to write logs from this command",
                           action="store", required=False, default="gecos_resp_prep.log")

    resp_prep.add_argument("-p", "--package", dest="qmpackage", choices=['g16'],
                           help="QM package to analyse (Default: Gaussian16).\n ",
                           action="store", required=False, default="g16")

    resp_prep.add_argument("-d", "--logfolder", dest="logfolder",
                           help="Folder containing the QM outputs.\n ",
                           action="store", required=True, default=None)

    resp_prep.add_argument("--resp_software", dest="respsoftware", choices=['forpyred', 'multiwfn'],
                           help="R|Software to prepare inputs to get RESP charges:\n"
                                " - forpyred: RED python server\n"
                                " - multiwfn: Use Multiwfn software.\n",
                           action="store", required=True)
    resp_prep_exlusive = resp_prep.add_mutually_exclusive_group()
    resp_prep_exlusive.add_argument("--deltaenerthresh", dest="deltaenerthresh",
                                    help="Delta energy threshold to select conformers in kcal/mol",
                                    action="store", required=False, default=None)
    resp_prep_exlusive.add_argument("--nconfthresh", dest="nconfthresh",
                                    help="Delta energy threshold to select conformers in kcal/mol",
                                    action="store", required=False, default=None)



    # Arguments for clustering command
    clustering.add_argument("--log", dest="log", type=str,
                            help="Name of the file to write logs from this command",
                            action="store", required=False, default="gecos_clustering.log")

    clustering.add_argument("-p", "--package", dest="qmpackage", choices=['g16'],
                            help="QM package to analyse (Default: Gaussian16).\n ",
                            action="store", required=False, default="g16")

    clustering.add_argument("-d", "--logfolder", dest="logfolder",
                            help="Folder containing the QM outputs.\n ",
                            action="store", required=True, default=None)

    clustering.add_argument("--indxfile", dest="indxfile",
                            help="A file using the ndx format from "
                            "GROMACS to extract internal coordinates from the log"
                            "files.",
                            action="store", required=False, default=None)

    # Arguments for similarity command
    similarity.add_argument("--log", dest="log", type=str,
                            help="Name of the file to write logs from this command",
                            action="store", required=False, default="gecos_similarity.log")

    similarity.add_argument("-r", "--ref", dest="reference", type=str,
                            help="File name of the spestrum to be used as reference.",
                            action="store", required=True, default=None)

    similarity.add_argument("-t", "--target", dest="target", type=str,
                            help="File name(s) of the spestrum to be compared with the reference.",
                            action="store", required=True, default=None)

    args = parser.parse_args()

    if args.command == "spectrum" or args.command == "clustering":
        if not os.path.isdir(args.logfolder):
            msg = "\nERROR: Directory {} does not exist\n".format(args.logfolder)
            parser.error(msg)
            exit()

    # Default values for cutoffs in the hydrogen bond
    if args.command == "energy":
        if args.hbonds is not None:
            if len(args.hbonds) == 1:
                args.hbonds = [float(args.hbonds[0]), 140.0]
            elif len(args.hbonds) != 2:
                args.hbonds = [3.2, 140.0]
            else:
                args.hbonds = [float(args.hbonds[0]), float(args.hbonds[1])]

        if args.indxfile is not None:
            if not os.path.isfile(args.indxfile):
                msg = "ERROR: If you required extract internal coordinates a ndx file is needed.\n"
                msg += "ERROR: The file {} does not exist in the current directory.".format(args.indxfile)
                parser.error(msg)
                exit()

    if args.command == "exp_spectrum":
        if not os.path.isfile(args.fileexpspectrum):
            msg = "ERROR: The file {} does not exist in the current directory.".format(args.fileexpspectrum)
            parser.error(msg)
            exit()

        if args.smoothpeak is not None:
            if args.smoothpeak.upper() == "SG":
                if not args.smoothSGparam:
                    msg = "ERROR: SG is invoked. The smoothSG option must be present"
                    parser.error(msg)
                    exit()

        if args.baseline is not None:

            if args.baseline.upper() == "ASLS":
                if not args.baselineparams:
                    msg = "ERROR: Baseline ASLS method is invoked. The baselineparams option must be present"
                    parser.error(msg)
                    exit()
                elif len(args.baselineparams) != 2:
                    msg = "ERROR: Baseline ASLS method is invoked. The baselineparams option must have two values"
                    parser.error(msg)
                    exit()

            if args.baseline.upper() == "MIXTURE_MODEL":
                if not args.baselineparams:
                    msg = "ERROR: Baseline MIXTURE_MODEL method is invoked. The baselineparams option must be present"
                    parser.error(msg)
                    exit()
                elif len(args.baselineparams) != 2:
                    msg = "ERROR: Baseline MIXTURE_MODEL method is invoked. " \
                          "The baselineparams option must have two values"
                    parser.error(msg)
                    exit()

            if args.baseline.upper() == "2021MUNICH":
                if not args.baselineparams:
                    msg = "ERROR: Baseline 2021Munich method is invoked. The baselineparams option must be present"
                    parser.error(msg)
                    exit()

    if args.command == "similarity":
        if not os.path.isfile(args.reference):
            msg = "ERROR: The file {} does not exist in the current directory.".format(args.reference)
            parser.error(msg)
            exit()

    return args, parser


# =============================================================================
def main_app():

    # Parse arguments
    args, parser = parse_arguments()

    # Setup log file
    log = utils.init_logger(
        "GaussianGecos",
        fileoutput=args.log,
        append=False, inscreen=True)

    # Print header
    utils.print_header_analysis(logger=log)

    if args.command in ["energy", "spectrum", "resp_prep", "clustering"]:
        if args.qmpackage == "g16":
            workdir = args.logfolder
            g16 = GaussianGecos(workdir, ext="log", logger=log)

            if args.command == "energy":
                g16.extract_energy()
                g16.extract_rmsd()
                if args.clusterconf:
                    g16.cluster_conformers(energy_thr=float(args.clusterparam[0]),
                                           rot_constant_thr=float(args.clusterparam[1]),
                                           rmsd_thr=float(args.clusterparam[2]),
                                           window_energy=float(args.clusterparam[3]),
                                           rmsd_only_heavy=bool(args.clusterparam[4]))
                g16.close_contacts(args)
                if args.rg:
                    g16.radius_of_gyration()
                g16.moment_of_inertia()
                if args.indxfile is not None:
                    g16.extract_internalcoords(args)
                g16.write_to_log(workdir, generate_data_gnuplot=True)

            elif args.command == "spectrum":
                g16.extract_vibrational_ir()
                g16.extract_rmsd()
                if args.clusterconf:
                    _, conftype = g16.cluster_conformers(energy_thr=float(args.clusterparam[0]),
                                                         rot_constant_thr=float(args.clusterparam[1]),
                                                         rmsd_thr=float(args.clusterparam[2]),
                                                         window_energy=float(args.clusterparam[3]),
                                                         rmsd_only_heavy=bool(args.clusterparam[4]))
                else:
                    conftype = None
                vf = g16.getvibfreqs()
                vi = g16.getvibirs()
                temp = g16.gettemperature()
                deltag = g16.getdeltag()
                nspectra = len(vf)
                spectra_ir = GecosIRSpectra(nspectra, vf, vi, temp, deltag,
                                            scale=args.scale,
                                            start=args.start,
                                            end=args.end,
                                            npoints=args.npoints,
                                            width=args.width,
                                            function=args.function,
                                            similarity=args.similarity,
                                            normalize=args.normalize,
                                            logger=log)
                spectra_ir.calculate_spectrum(conftype=conftype, avg=args.averaged)
                g16.write_vib_to_log(workdir, generate_data_gnuplot=True)
                if args.similarity:
                    spectra_ir._write_similarity_to_log()

            elif args.command == "resp_prep":

                try:
                    g16.extract_vibrational_ir()
                    iswrite_deltag = True
                except AttributeError:
                    g16.extract_energy()
                    iswrite_deltag = False
                g16.extract_rmsd()


                if iswrite_deltag:
                    g16.write_vib_to_log(workdir, generate_data_gnuplot=True)
                else:
                    g16.write_resp_to_log(workdir, generate_data_gnuplot=True)

                is_prepareinputs = True
                resp_software = args.respsoftware
                if args.deltaenerthresh:
                    GecosRespAnalysis(workdir, g16, is_prepareinputs, resp_software,
                                      deltaenergythreshold=float(args.deltaenerthresh), logger=log)
                elif args.nconfthresh:
                    GecosRespAnalysis(workdir, g16, is_prepareinputs, resp_software,
                                      nconfsthreshold=int(args.nconfthresh), logger=log)

            elif args.command == "clustering":
                g16.extract_energy()
                g16.extract_rmsd()
                if args.indxfile is not None:
                    if not os.path.isfile(args.indxfile):
                        msg = "ERROR: If you required extract internal coordinates a ndx file is needed.\n"
                        msg += "ERROR: The file {} does not exist in the current directory.".format(args.indxfile)
                        parser.error(msg)
                        exit()
                    g16.extract_internalcoords(args)

                data = np.array([g16._df['RMSDwithoutHs'], g16._df['DeltaEnergy(kcal/mol)']]).transpose()
                cnnclustering = CnnClusteringWrapper(data, logger=log)
                cnnclustering.withdata(radius_cutoff=0.5, similarity_cutoff=0, isdraw="png")

                # distances = pairwise_distances(data)
                # print(distances)
                #
                #
                #
                #
                #
                # #clustering = cluster.Clustering(data)
                # clustering = cluster.Clustering(distances, recipe="distances")
                # print(clustering.fit(radius_cutoff=2.0, similarity_cutoff=2))



                g16.write_clustering_to_log(workdir)

    elif args.command in ["exp_spectrum"]:
        _ = GecosExpIRSpectra(args.fileexpspectrum,
                              resample_list=args.resample_list,
                              normalize=args.normalize,
                              baseline=args.baseline,
                              baselineparam=args.baselineparams,
                              peaksmoothing=args.smoothpeak,
                              peaksmoothingparam=args.smoothSGparam,
                              logger=log)

    elif args.command in ["similarity"]:
        GecosCompareSpectrum(args.reference,
                             args.target,
                             logger=log)

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\n\t\tJob  Done at {} ============\n".format(now)
    print(m) if log is None else log.info(m)


# =============================================================================
if __name__ == "__main__":
    main_app()
