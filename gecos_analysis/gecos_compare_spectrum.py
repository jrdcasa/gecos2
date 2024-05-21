import glob
import os.path
import scipy
import numpy as np
from collections import defaultdict


class GecosCompareSpectrum:

    """
    This class compares reference experimental or theoretical spectrum to others spectra.
    """

    # =========================================================================
    def __init__(self, ref, target, logger=None):

        """
        Initialize the instance

        :param logger: Logger to write results
        """

        self._logger = logger
        self._filename_ref = ref
        self._filenamelist_target = glob.glob(target)

        self._inte_spectra_ref = None
        self._freq_spectra_ref = None

        self._ref_intensity = []
        self._ref_frequency = []
        self._npoints_ref = None

        self._similarity_matrix_pearson = defaultdict(list)
        self._similarity_matrix_spearman = defaultdict(list)

        # Calculate
        self._similarity()
        self._write_similarity_to_log()

    # =========================================================================
    def _similarity(self):

        """
        Check if the dimensions and the boundaries of the spectra are equals.

        :return:
        """

        # Reference spectrum ==========================

        with open(self._filename_ref, 'r') as fref:
            lines = fref.readlines()
            for iline in lines:
                if iline.count("#") != 0:
                    continue
                try:
                    freq, inten, _ = iline.split()
                except ValueError:
                    freq, inten = iline.split()
                self._ref_frequency.append(float(freq))
                self._ref_intensity.append(float(inten))

        if len(self._ref_frequency) != len(self._ref_intensity):
            msg = "\n\t\tThe number of points for intensity and fequency must be the same.\n"
            msg += "\t\t  File = {}\n".format(self._filename_ref)
            msg += "\t\t  Frq number points = {}\n".format(len(self._ref_frequency))
            msg += "\t\t  Int number points = {}\n".format(len(self._ref_intensity))
            print(m) if self._logger is None else self._logger.error(msg)
            exit()

        self._npoints_ref = len(self._ref_frequency)

        # Target spectra ================================
        nspectra = len(self._filenamelist_target)

        for ispectra_name in self._filenamelist_target:
            target_intensity = []
            target_frequency = []
            key = os.path.splitext(os.path.split(ispectra_name)[-1])[0]
            with open(ispectra_name, 'r') as ftarget:
                lines = ftarget.readlines()
                for iline in lines:
                    if iline.count("#") != 0:
                        continue
                    try:
                        freq, inten, _ = iline.split()
                    except ValueError:
                        freq, inten = iline.split()
                    target_frequency.append(float(freq))
                    target_intensity.append(float(inten))

                self._check_dimensions_spectra(ispectra_name, target_frequency, target_intensity)

                pearson = scipy.stats.pearsonr(target_intensity, self._ref_intensity)
                spearman = scipy.stats.spearmanr(target_intensity, self._ref_intensity)

                self._similarity_matrix_pearson[key] = pearson.statistic
                self._similarity_matrix_spearman[key] = spearman.correlation

    # =========================================================================
    def _check_dimensions_spectra(self, ispectra_name, target_frequency, target_intensity):

        """
        Check if the dimensions and the boundaries of the spectra are equals.

        :return:
        """

        if len(target_frequency) != len(target_intensity):
            msg = "\n\t\tThe number of points for intensity and frequency must be the same.\n"
            msg += "\t\t  File = {}\n".format(ispectra_name)
            msg += "\t\t  Frq number points = {}\n".format(len(target_frequency))
            msg += "\t\t  Int number points = {}\n".format(len(target_intensity))
            print(m) if self._logger is None else self._logger.error(msg)
            exit()

        if len(target_frequency) != len(self._ref_frequency):
            msg = "\n\t\tThe number of points for frequency in the reference " \
                  "and target spectra must be the same.\n"
            msg += "\t\t  File reference = {}\n".format(self._filename_ref)
            msg += "\t\t  File target    = {}\n".format(ispectra_name)
            msg += "\t\t  Frq number points ref.   = {}\n".format(len(self._ref_frequency))
            msg += "\t\t  Frq number points target = {}\n".format(len(target_frequency))
            print(m) if self._logger is None else self._logger.error(msg)
            exit()

        lll = [target_frequency[i] - self._ref_frequency[i] for i in range(0, self._npoints_ref)]
        if any(lll) > 1e-01:
            msg = "\n\t\tFrequencies in the reference " \
                  "and target spectra must be alligned.\n"
            msg += "\t\t  File reference = {}\n".format(self._filename_ref)
            msg += "\t\t  File target    = {}\n".format(ispectra_name)
            msg += "\t\t  Frq number points ref.   = {}\n".format(len(self._ref_frequency))
            msg += "\t\t  Frq number points target = {}\n".format(len(target_frequency))
            print(m) if self._logger is None else self._logger.error(msg)
            exit()

    # =========================================================================
    def _write_similarity_to_log(self):

        """
        Write similarity indices to log

        :return:
        """

        # Write table
        m = "\n\t\t ====== SIMILARITY OF CALCULATED AND EXPERIMENTAL SPECTRA =======\n"
        m += "\t\t{0:1s} {1:^40s} {2:^22s} {3:^22s}\n".\
            format('#', 'ID', 'Pearson correlation', "Spearman correlation")

        lenm = len(m)
        m += "\t\t# " + len(m) * "=" + "\n"
        yidx = 0
        for xkey, xvalue in self._similarity_matrix_pearson.items():
            line = "\t\t{0:^40s} {1:^26.6f} {2:^20.6f}\n".\
                format(xkey,
                       self._similarity_matrix_pearson[xkey],
                       self._similarity_matrix_spearman[xkey])
            yidx += 1

            m += line

        m = m[:-1]  # Remove last \n
        print(m) if self._logger is None else self._logger.info(m)

        m1 = "\t\t# " + lenm * "=" + "\n"
        print(m1) if self._logger is None else self._logger.info(m1)
