import glob
import os.path
from collections import defaultdict
import numpy as np
import pandas as pd
import scipy
import datetime


class GecosIRSpectra:

    """
    This class is use to get IR spectra from QM calculations.
    """

    # =========================================================================
    def __init__(self, nspectra, freqs, iractivity, temperature, deltag=None,
                 scale=None, start=0.0, end=4000.0, npoints=500, width=10.0,
                 function="lorentzian", similarity=True, normalize=None, logger=None):

        """
        Initialize the instance

        :param nspectra
        :param freqs
        :param iractivity
        :param temperature
        :param deltag
        :param scale
        :param start
        :param end
        :param npoints
        :param width
        :param function
        :param logger: Logger to write results
        """

        self._logger = logger
        self._nspectra = nspectra
        self._freqs_dict = freqs
        self._iractiivity_dict = iractivity
        self._temperature = temperature
        self._deltag = deltag
        self._start = start
        self._end = end
        self._npoints = npoints
        self._width = width
        self._function = function
        self._similarity = similarity
        self._isscaled = False
        self._nspectra = len(self._freqs_dict)
        self._scale = scale
        self._spectrum = defaultdict()
        self._xvalues = None
        self._scaled_boltz_spectrum = None
        self._scaled_simple_spectrum = None
        self._scaled_similarity_spectrum = None
        self._similarity_matrix_pearson = None
        self._similarity_matrix_spearman = None
        self._tonormalize = normalize
        self._isnormalized = False
        self._nameID_to_idx = defaultdict()
        self._idx_to_nameID = defaultdict()

        self._pearson_spearman_boltzmann_lowestenergy = None
        self._pearson_spearman_simple_lowestenergy = None
        self._pearson_spearman_boltzmann_simple = None
        self._pearson_spearman_similiraty_boltzmann = None
        self._pearson_spearman_similiraty_simple = None
        self._pearson_spearman_similiraty_lowestenergy = None

        # Scale frequencies
        if scale is not None:
            self._scale_frequencies()

    # =========================================================================
    def _scale_frequencies(self):

        """
        Scale frequency positions.

        :return:
        """

        for ikey, ivalue in self._freqs_dict.items():
            self._freqs_dict[ikey] = ivalue * self._scale
        self._isscaled = True

    # =========================================================================
    def calculate_spectrum(self, conftype=None, avg=None):

        """
        Calculate spectrum

        :return:
        """

        # Zip in a list the freq and intensity for each spectra
        peaks = defaultdict()
        for ikey, freqvalue in self._freqs_dict.items():
            irvalue = self._iractiivity_dict[ikey]
            peaks[ikey] = list(zip(freqvalue, irvalue))

        if self._tonormalize is not None:
            now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
            m = "\n\t\t ======== NORMALIZING DATA ======== ({})\n".format(now)
            m += "\t\t Normalizing data using the method: {}".format(self._tonormalize)
            print(m) if self._logger is None else self._logger.info(m)

        # Create the spectrum and save the data
        ind = 0
        for ikey, ivalue in self._freqs_dict.items():
            self._nameID_to_idx[ikey] = ind
            self._idx_to_nameID[ind] = ikey
            self._spectrum[ikey] = np.zeros(self._npoints)
            if self._function.upper() == "LORENTZIAN":
                self._xvalues = np.arange(self._npoints)*float(self._end-self._start)/(self._npoints-1)+self._start
            elif self._function.upper() == "LORENTZIAN_GAUSSIAN":
                self._xvalues = np.arange(self._npoints)*float(self._end-self._start)/self._npoints+self._start
            for i in range(self._npoints):
                x = self._xvalues[i]
                for (pos, height) in peaks[ikey]:
                    self._spectrum[ikey][i] = self._spectrum[ikey][i] + \
                                           self._formula(self._function, x, pos, height, self._width)

            # Normalize spectra
            if self._tonormalize is not None:
                self._normalize(method=self._tonormalize)

            filename = os.path.splitext(os.path.split(ikey)[-1])[0]+"_spectrum.dat"
            with open(filename, 'w') as fdat:
                line = ""
                max_intensity = np.max(self._spectrum[ikey])
                range_freq = float(self._end-self._start)
                for i in range(self._npoints):
                    if self._function.upper() == "LORENTZIAN":
                        realx = range_freq*(i+1)/self._npoints+self._start
                    elif self._function.upper() == "LORENTZIAN_GAUSSIAN":
                        # realx = self._xvalues[i]
                        realx = range_freq*(i+1)/self._npoints+self._start
                    line += "{} {} {}\n".format(realx, self._spectrum[ikey][i], self._spectrum[ikey][i]/max_intensity)
                fdat.writelines(line)

            ind += 1

        # GNUPLOT templates
        self._gnuplot_template()

        # Average spectrum
        self._avg_spectrum(avg, conftype)

        # Similary of the spectra
        # The spectra with lowest energy is used as template
        if self._similarity:
            self._calculate_similarity()
            self._avg_spectrum_by_similarity(conftype)

    # =========================================================================
    def _calculate_similarity(self):

        """
        Assesment of the similarity of the calculated QM spectra using the
        method discussed in https://dx.doi.org/10.1021/acs.jctc.0c00126

        In self._similarity_matrix the element [i,j] represents the similarity
        coefficient of the spectra i and j.
        """

        nspectra = len(self._freqs_dict)
        self._similarity_matrix_pearson = np.zeros([nspectra, nspectra], dtype=np.float32)
        self._similarity_matrix_spearman = np.zeros([nspectra, nspectra], dtype=np.float32)
        xidx = 0
        for xkey, xvalue in self._spectrum.items():
            yidx = 0
            for ykey, yvalue in self._spectrum.items():
                pearson = scipy.stats.pearsonr(xvalue, yvalue)
                spearman = scipy.stats.spearmanr(xvalue, yvalue)
                self._similarity_matrix_pearson[xidx, yidx] = pearson.statistic
                self._similarity_matrix_spearman[xidx, yidx] = spearman.correlation
                yidx += 1
            xidx += 1

        return None

    # =========================================================================
    def _write_similarity_to_log(self):

        """
        Write similarity indices to log

        :return:
        """

        # Write table
        m = "\n\t\t{0:1s} {1:^40s} {2:^22s} {3:^22s}\n".\
            format('#', 'ID', 'Pearson correlation', "Spearman correlation")

        lenm = len(m)
        m += "\t\t# " + len(m) * "=" + "\n"
        yidx = 0
        for xkey, xvalue in self._freqs_dict.items():
            line = "\t\t{0:^40s} {1:^26.6f} {2:^20.6f}\n".\
                format(xkey,
                       self._similarity_matrix_pearson[0, yidx],
                       self._similarity_matrix_spearman[0, yidx])
            yidx += 1

            m += line

        m = m[:-1]  # Remove last \n
        print(m) if self._logger is None else self._logger.info(m)

        m1 = "\t\t# " + lenm * "=" + "\n"
        print(m1) if self._logger is None else self._logger.info(m1)

        # Averaged similarities
        m = '\n'
        if self._pearson_spearman_boltzmann_lowestenergy is not None:
            m += '\t\t Pearson and Spearman coefficients between Boltzmman and Lowest energy  = {0:^9.6f} {1:^9.6f}\n'\
                .format(self._pearson_spearman_boltzmann_lowestenergy[0],
                        self._pearson_spearman_boltzmann_lowestenergy[1])
        if self._pearson_spearman_simple_lowestenergy is not None:
            m += '\t\t Pearson and Spearman coefficients between Simple and Lowest energy     = {0:^9.6f} {1:^9.6f}\n'\
                .format(self._pearson_spearman_simple_lowestenergy[0],
                        self._pearson_spearman_simple_lowestenergy[1])

        if self._pearson_spearman_boltzmann_simple is not None:
            m += '\t\t Pearson and Spearman coefficients between Boltzmman and Simple average = {0:^9.6f} {1:^9.6f}\n'\
                .format(self._pearson_spearman_boltzmann_simple[0],
                        self._pearson_spearman_boltzmann_simple[1])

        if self._pearson_spearman_similiraty_boltzmann is not None:
            m += '\t\t Pearson and Spearman coefficients between Similarity and Boltzmman     = {0:^9.6f} {1:^9.6f}\n'\
                .format(self._pearson_spearman_similiraty_boltzmann[0],
                        self._pearson_spearman_similiraty_boltzmann[1])

        if self._pearson_spearman_similiraty_simple is not None:
            m += '\t\t Pearson and Spearman coefficients between Similarity and Simple avg.   = {0:^9.6f} {1:^9.6f}\n'\
                .format(self._pearson_spearman_similiraty_simple[0],
                        self._pearson_spearman_similiraty_simple[1])

        if self._pearson_spearman_similiraty_lowestenergy is not None:
            m += '\t\t Pearson and Spearman coefficients between Similarity and Lowest energy = {0:^9.6f} {1:^9.6f}\n'\
                .format(self._pearson_spearman_similiraty_lowestenergy[0],
                        self._pearson_spearman_similiraty_lowestenergy[1])

        lenm = len(m.split("\n")[1])
        m1 = "\t\t# " + lenm * "-"
        print(m1+m+m1+"\n") if self._logger is None else self._logger.info(m1+m+m1+"\n")

    # =========================================================================
    def _formula(self, name, x, peak, height, width):
        """The lorentzian curve.

        f(x) = a/(1+a)

        where a is FWHM**2/4
        """

        if name.upper() == "LORENTZIAN":
            a = width**2./4.
            return float(height)*a/((peak-x)**2 + a)
        elif name.upper() == "LORENTZIAN_GAUSSIAN":
            """ See DOI: 10.13140/RG.2.1.4181.6160"""
            return float(2. / np.pi) * (100 / np.log(10)) * float(height) * (
                        width / ((4 * (peak - x) ** 2) + width ** 2))
        else:
            m = "\t\t ERROR: {} formula is not implemented".format(name)
            print(m) if self._logger is not None else self._logger.error(m)

    # =========================================================================
    def _gnuplot_template(self, nx=3, ny=3):

        max_intensity = 0
        for ikey, ir in self._spectrum.items():
            if np.max(ir) > max_intensity:
                max_intensity = np.max(ir)
        if not self._isnormalized:
            max_intensity += 10
        max_intensity = np.floor(max_intensity)
        if max_intensity < 1:
            max_intensity = 1

        linegnuplot = "reset\n"
        linegnuplot += 'set xlabel "Freq (cm^-1)"\n'
        if not self._isnormalized:
            linegnuplot += 'set ylabel "Intensity"\n'
        else:
            linegnuplot += 'set ylabel "Normalized Intensity"\n'
        linegnuplot += 'set grid\n'
        linegnuplot += 'set xrange [{0:d}:{1:d}]\n'.format(int(self._start), int(self._end))
        linegnuplot += 'set yrange [{0:d}:{1:d}]\n'.format(0, int(max_intensity))
        linegnuplot += "# "+50*'*'+"\n"

        iwxt = 1
        iplot = 0
        for ikey, ir_values in self._spectrum.items():
            ifilename = os.path.splitext(os.path.split(ikey)[-1])[0]+"_spectrum.dat"
            if iplot % (nx*ny) == 0:
                linegnuplot += "unset multiplot\n"
                iheight = int(500 * ny)
                iwidth = int(300 * nx)
                linegnuplot += 'set term wxt {0:d} enhanced dashed size {1:d},{2:d} ' \
                               'font "Arial,8"\n'.format(iwxt, iheight, iwidth)
                linegnuplot += 'set multiplot layout {0:d},{1:d}\n'.format(nx, ny)
                iwxt += 1
            ititle = ifilename.replace("_", "-")
            linegnuplot += 'set title "{0:s}"\n'.format(ititle)
            linegnuplot += 'p "./{0:s}" u 1:2 w l notitle\n'.format(ifilename)
            iplot += 1

        with open("template_gnu.plot", "w") as fgnuplot:
            fgnuplot.writelines(linegnuplot)

    # =========================================================================
    def _avg_spectrum(self, avg, conftype):

        """
        Averaging spectra
        """

        # ===============================================
        def boltzmann_weight(conformertype):

            """
            Calculate the boltzmann averaged spectrum
            :return:
            """

            nspectra_avg = 0
            kb_t = 1.987204259/1000  # kcal/mol
            sum_weight = np.zeros(self._npoints)
            sum_bfactor = 0.0

            for iikey, ir_values in self._spectrum.items():
                if conftype is not None:
                    if conformertype[iikey].count("Rotamer") != 0 or conformertype[iikey].count("Identical") != 0:
                        continue
                boltzmann_factor = np.exp(-self._deltag[iikey] / (kb_t * self._temperature[iikey]))
                sum_weight += boltzmann_factor * self._spectrum[iikey]
                sum_bfactor += boltzmann_factor
                nspectra_avg += 1

            self._scaled_boltz_spectrum = sum_weight / sum_bfactor

            filename = "boltzmann_spectrum.dat"
            title = "# Boltzmann Averaged spectrum #spectra: {}\n".format(nspectra_avg)

            range_freq = self._end - self._start
            with open(filename, 'w') as fdat:
                line = title
                for i in range(self._npoints):
                    realx = range_freq * (i + 1) / self._npoints + self._start
                    # line += "{} {}\n".format(self._xvalues[i], self._scaled_boltz_spectrum[i])
                    line += "{} {}\n".format(realx, self._scaled_boltz_spectrum[i])
                fdat.writelines(line)

            return nspectra_avg

        # ===============================================
        def simple_weight(conformertype):

            """
            Calculate the arithmethic averaged spectrum
            :return:
            """

            nspectra_avg = 0
            sum_weight = np.zeros(self._npoints)

            for iikey, ir_values in self._spectrum.items():
                if conftype is not None:
                    if conformertype[iikey].count("Rotamer") != 0 or conformertype[iikey].count("Identical") != 0:
                        continue
                boltzmann_factor = 1.0
                sum_weight += boltzmann_factor * self._spectrum[iikey]
                nspectra_avg += 1

            self._scaled_simple_spectrum = sum_weight / nspectra_avg

            filename = "averaged_spectrum.dat"
            title = "# Averaged spectrum (Non boltzmann) # spectra: {}\n".format(nspectra_avg)

            range_freq = self._end - self._start
            with open(filename, 'w') as fdat:
                line = title
                for i in range(self._npoints):
                    realx = range_freq * (i + 1) / self._npoints + self._start
                    # line += "{} {}\n".format(self._xvalues[i], self._scaled_simple_spectrum[i])
                    line += "{} {}\n".format(realx, self._scaled_simple_spectrum[i])
                fdat.writelines(line)

            return nspectra_avg

        # ===============================================
        def gnu_template_avg(nspectra_avg, nx=1, ny=1):

            tmpdict = defaultdict(list)
            if avg == "boltzmann" and self._deltag is not None:
                tmpdict["boltzmann"] = self._scaled_boltz_spectrum
            elif avg == "simple":
                tmpdict["averaged"] = self._scaled_simple_spectrum
            else:
                tmpdict["boltzmann"] = self._scaled_boltz_spectrum
                tmpdict["averaged"] = self._scaled_simple_spectrum

            max_intensity = 0
            for _, ir in tmpdict.items():
                if np.max(ir) > max_intensity:
                    max_intensity = np.max(ir)
            if not self._isnormalized:
                max_intensity += 10
            max_intensity = np.floor(max_intensity)
            if max_intensity < 1:
                max_intensity = 1

            linegnuplot = "reset\n"
            linegnuplot += 'set xlabel "Freq (cm^-1)"\n'
            if not self._isnormalized:
                linegnuplot += 'set ylabel "Intensity"\n'
            else:
                linegnuplot += 'set ylabel "Normalized Intensity"\n'
            linegnuplot += 'set grid\n'
            linegnuplot += 'set xrange [{0:d}:{1:d}]\n'.format(int(self._start), int(self._end))
            linegnuplot += 'set yrange [{0:d}:{1:d}]\n'.format(0, int(max_intensity))
            linegnuplot += "# "+50*'*'+"\n"

            iwxt = 1
            iplot = 0
            for iikey, ir_values in tmpdict.items():
                ifilename = os.path.splitext(os.path.split(iikey)[-1])[0]+"_spectrum.dat"
                if iplot % (nx*ny) == 0:
                    iheight = int(500 * ny)
                    iwidth = int(300 * nx)
                    linegnuplot += 'set term wxt {0:d} enhanced dashed size {1:d},{2:d} ' \
                                   'font "Arial,8"\n'.format(iwxt, iheight, iwidth)
                    linegnuplot += 'set multiplot layout {0:d},{1:d}\n'.format(nx, ny)
                    iwxt += 1
                ititle = ifilename.replace("_", "-")
                linegnuplot += 'set title "{0:s} # spectra {1:d}"\n'.format(ititle, nspectra_avg)
                linegnuplot += 'p "./{0:s}" u 1:2 w l notitle\n'.format(ifilename)
                linegnuplot += "unset multiplot\n"
                iplot += 1

            with open("template_avg_gnu.plot", "w") as fgnuplot:
                fgnuplot.writelines(linegnuplot)

        # ========================================================
        key_lowest_energy = None
        for ikey, _ in self._deltag.items():
            key_lowest_energy = ikey
            break

        if avg == "boltzmann" and self._deltag is not None:
            xvalue = self._scaled_boltz_spectrum
            nspectra = boltzmann_weight(conftype)
            zvalue = self._spectrum[key_lowest_energy]
            pearson = scipy.stats.pearsonr(xvalue, zvalue)
            spearman = scipy.stats.spearmanr(xvalue, zvalue)
            self._pearson_spearman_boltzmann_lowestenergy = [pearson.statistic, spearman.correlation]
        elif avg == "simple":
            nspectra = simple_weight(conftype)
            yvalue = self._scaled_simple_spectrum
            zvalue = self._spectrum[key_lowest_energy]
            pearson = scipy.stats.pearsonr(yvalue, zvalue)
            spearman = scipy.stats.spearmanr(yvalue, zvalue)
            self._pearson_spearman_simple_lowestenergy = [pearson.statistic, spearman.correlationn]
        else:
            boltzmann_weight(conftype)
            nspectra = simple_weight(conftype)
            xvalue = self._scaled_boltz_spectrum
            yvalue = self._scaled_simple_spectrum
            pearson = scipy.stats.pearsonr(xvalue, yvalue)
            spearman = scipy.stats.spearmanr(xvalue, yvalue)
            self._pearson_spearman_boltzmann_simple = [pearson.statistic, spearman.correlation]
            zvalue = self._spectrum[key_lowest_energy]
            pearson = scipy.stats.pearsonr(xvalue, zvalue)
            spearman = scipy.stats.spearmanr(xvalue, zvalue)
            self._pearson_spearman_boltzmann_lowestenergy = [pearson.statistic, spearman.correlation]
            pearson = scipy.stats.pearsonr(yvalue, zvalue)
            spearman = scipy.stats.spearmanr(yvalue, zvalue)
            self._pearson_spearman_simple_lowestenergy = [pearson.statistic, spearman.correlation]

        gnu_template_avg(nspectra)

    # ===============================================
    def _avg_spectrum_by_similarity(self, conformertype, threshold=0.95):

        nspectra_avg = 0
        sum_weight = np.zeros(self._npoints)

        # Take indices of the spectra
        nrows, ncols = self._similarity_matrix_pearson.shape
        to_add = [0]
        excluded = np.zeros(nrows, dtype=np.bool)
        for irow in range(0, nrows):
            for icol in range(irow+1, ncols):
                if excluded[icol]:
                    continue
                if self._similarity_matrix_pearson[irow, icol] >= threshold:
                    excluded[icol] = True
                    if icol in to_add:
                        to_add.remove(icol)
                    continue
                if icol not in to_add:
                    to_add.append(icol)

        # Calculate the average similariry spectra
        lines = ""
        for i in to_add:
            iikey = self._idx_to_nameID[i]
            if conformertype is not None:
                if conformertype[iikey].count("Rotamer") != 0 or conformertype[iikey].count("Identical") != 0:
                    continue
            boltzmann_factor = 1.0
            sum_weight += boltzmann_factor * self._spectrum[iikey]
            nspectra_avg += 1
            lines += iikey +"\n"
        self._scaled_similarity_spectrum = sum_weight / nspectra_avg

        with open("similarity_avg_spectrum_used.dat", 'w') as f:
            f.writelines(lines)

        # Write the spectra
        tmpdict = defaultdict(list)
        tmpdict["similarity"] = self._scaled_similarity_spectrum

        filename = "similarity_avg_spectrum.dat"
        title = "# Averaged spectrum by similarity (Non boltzmann) # spectra: {}\n".format(nspectra_avg)

        with open(filename, 'w') as fdat:
            line = title
            for i in range(self._npoints):
                line += "{} {}\n".format(self._xvalues[i], self._scaled_similarity_spectrum[i])
            fdat.writelines(line)

        iwxt = 3
        nx = 1
        ny = 1
        iheight = int(500 * ny)
        iwidth = int(300 * nx)
        linegnuplot = 'set term wxt {0:d} enhanced dashed size {1:d},{2:d} ' \
                      'font "Arial,8"\n'.format(iwxt, iheight, iwidth)
        linegnuplot += 'set multiplot layout {0:d},{1:d}\n'.format(nx, ny)
        ititle = filename.replace("_", "-")
        linegnuplot += 'set title "{0:s} # spectra {1:d}"\n'.format(ititle, nspectra_avg)
        linegnuplot += 'p "./{0:s}" u 1:2 w l notitle\n'.format(filename)

        with open("template_avg_gnu.plot", "a") as fgnuplot:
            fgnuplot.writelines(linegnuplot)

        key_lowest_energy = None
        for ikey, _ in self._deltag.items():
            key_lowest_energy = ikey
            break

        xvalue = self._scaled_boltz_spectrum
        yvalue = self._scaled_simple_spectrum
        zvalue = self._spectrum[key_lowest_energy]
        avalue = self._scaled_similarity_spectrum

        pearson = scipy.stats.pearsonr(xvalue, avalue)
        spearman = scipy.stats.spearmanr(xvalue, avalue)
        self._pearson_spearman_similiraty_boltzmann = [pearson.statistic, spearman.correlation]

        pearson = scipy.stats.pearsonr(yvalue, avalue)
        spearman = scipy.stats.spearmanr(yvalue, avalue)
        self._pearson_spearman_similiraty_simple = [pearson.statistic, spearman.correlation]

        pearson = scipy.stats.pearsonr(zvalue, avalue)
        spearman = scipy.stats.spearmanr(zvalue, avalue)
        self._pearson_spearman_similiraty_lowestenergy = [pearson.statistic, spearman.correlation]

    # =========================================================================
    def _normalize(self, method="minmax"):

        method = method.upper()
        self._isnormalized = True

        for ikey, item in self._spectrum.items():
            nx_data = len(self._xvalues)
            ny_data = len(self._spectrum[ikey])
            assert (nx_data == ny_data)
            data = np.array(list(zip(self._xvalues, self._spectrum[ikey])))

            # if method == "MINMAX":
            #     norm_method = "Min-Max normalization (mM)"
            # elif method == "VECTORNORM":
            #     norm_method = "Vector normalization (VN)"
            # elif method == "1-NORM":
            #     norm_method = "1-norm (1-n)"
            # elif method == "SNV":
            #     norm_method = "Standard Normal Variate (SNV)"
            # else:
            #     norm_method = None

            if method == "MINMAX":
                data = self.__normalize_minmax_method(data)
            elif method == "VECTORNORM":
                data = self.__normalize_vectornorm_method(data)
            elif method == "1-NORM":
                data = self.__normalize_onenorm_method(data)
            elif method == "SNV":
                data = self.__normalize_snv_method(data)
            else:
                data = None

            self._spectrum[ikey] = data[:, 1]

    # ========================================================================
    @staticmethod
    def __normalize_minmax_method(data):

        """
        Min-Max normalization

        y_norm(i) = (y(i) - y_min) / (y_max - y_min)

        :return: Normalized data
        """

        shape_x, shape_y = data.shape

        data_norm = np.zeros([shape_x, shape_y], dtype=np.float32)

        y_max = np.max(data, axis=0)[1]
        y_min = np.min(data, axis=0)[1]
        delta_y = y_max - y_min

        for idx, i in enumerate(data):
            data_norm[idx, 0] = data[idx, 0]
            data_norm[idx, 1] = (i[1] - y_min)/delta_y

        return data_norm

    # ========================================================================
    @staticmethod
    def __normalize_vectornorm_method(data):

        """
        Vector normalization

        y_norm(i) = y(i) / sqrt(y(1)^2+y(2)^2+...+y(n)^2), n is the number of data

        :return: Normalized data
        """

        shape_x, shape_y = data.shape

        data_norm = np.zeros([shape_x, shape_y], dtype=np.float32)

        value = np.linalg.norm(data[:, 1], ord=2)

        for idx, i in enumerate(data):
            data_norm[idx, 0] = data[idx, 0]
            data_norm[idx, 1] = i[1] / value

        return data_norm

    # ========================================================================
    @staticmethod
    def __normalize_onenorm_method(data):

        """
        1-norm (1-n)

        y_mean = Sum(1..n) y(i) / n
        y_norm(i) = (y(i) - y_mean) / Sum(1..n)|y(i) - y_mean|, n is the number of data

        :return: Normalized data
        """

        shape_x, shape_y = data.shape

        data_norm = np.zeros([shape_x, shape_y], dtype=np.float32)

        mean = np.mean(data[:, 1])

        denominator_sum = 0.0
        for idx, i in enumerate(data):
            denominator_sum += np.abs(i[1]-mean)

        for idx, i in enumerate(data):
            data_norm[idx, 0] = data[idx, 0]
            data_norm[idx, 1] = (i[1] - mean) / denominator_sum

        return data_norm

    # ========================================================================
    @staticmethod
    def __normalize_snv_method(data):

        """
        1-norm (1-n)

        SD = sqrt( [1/N-1]Sum(1,N) (Ai-A)^2)
        y_norm(i) = (y(i) - y_mean) / SD, n is the number of data

        :return: Normalized data
        """

        shape_x, shape_y = data.shape

        data_norm = np.zeros([shape_x, shape_y], dtype=np.float32)

        mean = np.mean(data[:, 1])
        std = np.std(data[:, 1])

        for idx, i in enumerate(data):
            data_norm[idx, 0] = data[idx, 0]
            data_norm[idx, 1] = (i[1] - mean) / std

        return data_norm
