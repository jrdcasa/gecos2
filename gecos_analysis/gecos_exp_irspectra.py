import numpy as np
import datetime
import os
from scipy.interpolate import interp1d
from scipy import sparse
from scipy import signal
from pybaselines.whittaker import asls, iasls
from pybaselines.spline import mixture_model
from pybaselines.optimizers import optimize_extended_range


# noinspection PyTypeChecker
class GecosExpIRSpectra:

    """
    This class is use to get IR spectra from QM calculations.
    """

    # =========================================================================
    def __init__(self, filename, resample_list=None, normalize=None,
                 baseline=None, baselineparam=None,
                 peaksmoothing=None, peaksmoothingparam=None, logger=None):

        """
        Initialize the instance and perform the calculations on the spectrum.

        :param filename: Name of the file containing the experimental spectra. It is assumed that the file
                         has two columns: the first column is, the second column.
        :param delta_x: Jump between x values in the resample method.
        :param normalize: Method to normalize.
        :param baseline: Method to create the baseline.
        :param peaksmothing: Method to smooth peaks.
        :param logger: Logger to write results

        """

        # Logger
        self._logger = logger
        # Name of file containing the original spectra
        self._filename = filename
        # Numpy object containing the original data
        self._original_data = None
        # Pattern to be used in the name of new files
        self._pattern = os.path.splitext(os.path.split(filename)[-1])[0]
        # Resample the data using this delta_x
        self._normalize_method = normalize
        self._resample_list = resample_list
        # Baseline methods
        self._baseline_method = baseline
        self._baseline_data = None
        self._baselineparams = baselineparam
        # Peak Smoothing
        self._peaksmoothing_method = peaksmoothing
        self._peaksmoothSGparam = peaksmoothingparam

        # Read the file  --> Setup self._original_data
        self._current_data = self._get_data_from_file()

        # Resample spectrum
        if self._resample_list is not None:
            self._current_data = self._resample_data()
            a = np.int(self._resample_list[0])
            b = np.int(self._resample_list[1])
            c = np.int(self._resample_list[2])
            fname_tmp = "00_tmp_"+self._pattern+"_resampled_{}_{}_{}".format(a, b, c)+".dat"
            np.savetxt(fname_tmp, self._current_data, delimiter=" ", fmt='%5.1f %.18f')

        # Normalize spectrum
        if self._normalize_method is not None:
            self._current_data = self._normalize_data(self._current_data)
            fname_tmp = "01_tmp_"+self._pattern+"_normalize_{}".format(self._normalize_method)+".dat"
            np.savetxt(fname_tmp, self._current_data, delimiter=" ", fmt='%5.1f %.18f')

        # Baseline correction
        if self._baseline_method is not None:
            self._baseline_data = self._create_baseline(self._current_data, baseline,
                                                        params=self._baselineparams)
            # Get the final spectra
            ndata = self._current_data.shape[0]
            for i in range(0, ndata):
                point = self._current_data[i, 1] - self._baseline_data[i, 1]
                if point < 0.0:
                    self._current_data[i, 1] = 0.0
                else:
                    self._current_data[i, 1] = point
            fname_tmp = "04-tmp_baselinecorrected_{}.dat".format(self._baseline_method)
            np.savetxt(fname_tmp, self._current_data, delimiter=" ", fmt='%5.1f %.18f')

        # Peak Smoothing
        if self._peaksmoothing_method is not None:
            y = self._current_data[:, 1]
            window_width = int(self._peaksmoothSGparam[0])
            polyorder = int(self._peaksmoothSGparam[1])
            self._current_data[:, 1] = signal.savgol_filter(y, window_width, polyorder)
            fname_tmp = "05-tmp_peaksmoothing_{}.dat".format(self._peaksmoothing_method)
            np.savetxt(fname_tmp, self._current_data, delimiter=" ", fmt='%5.1f %.18f')

        fname_final = self._pattern + "_corrected.dat"
        np.savetxt(fname_final, self._current_data, delimiter=" ", fmt='%5.1f %.18f')

        # Create gnuplot template
        self._gnuplot_template(nx=2, ny=3)

    # =========================================================================
    def _get_data_from_file(self):

        """
        Extract data from the file and setup in the self._original_data array
        """

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\n\t\t ======== GETTING DATA ======== ({})\n".format(now)
        m += "\t\t Getting data from {}".format(self._filename)
        print(m) if self._logger is None else self._logger.info(m)

        ndata = 0

        with open(self._filename, 'r') as fin:
            lines = fin.readlines()
            # Find number of data
            for iline in lines:
                # Comment lines or lines with spaces
                if len(iline[:-1].strip()) == 0:
                    continue
                if iline.count("#") == 0:
                    ndata += 1
            self._original_data = np.zeros([ndata, 2], dtype=np.float32)
            idx = 0
            for iline in lines:
                try:
                    if iline.count("#") == 0:
                        freq, intensity = iline.split()
                        self._original_data[idx, 0] = float(freq)
                        self._original_data[idx, 1] = float(intensity)
                        idx += 1
                except ValueError:  # Blank lines
                    pass

        return self._original_data

    # ========================================================================
    def _resample_data(self):
        """

        Get an array of dimension (ndata,2) and resample it using dx

        :return: A resampled array
        """

        if self._original_data.shape[1] != 2:
            m = "\n\t\t Resampling cannot be done for a array of [{},{}] dimensions.\n".\
                format(self._original_data.shape[0], self._original_data.shape[1])
            m += "\t\t Only arrays of dimension [:, 2] can be resampled."
            print(m) if self._logger is None else self._logger.error(m)
            exit()

        x_min_orig = self._original_data.min(axis=0)[0]
        x_max_orig = self._original_data.max(axis=0)[0]
        ndata_orig = self._original_data.shape[0]
        deltax_orig = (x_max_orig - x_min_orig) / (ndata_orig - 1)

        x_min_resample = self._resample_list[0]
        x_max_resample = self._resample_list[1]
        npoints_resample = self._resample_list[2]

        x_data_orig = self._original_data[:, 0]
        y_data_orig = self._original_data[:, 1]
        funcinterp = interp1d(x_data_orig, y_data_orig, fill_value="extrapolate")
        x_data_tmp = np.linspace(x_min_resample, x_max_resample, int(npoints_resample))

        # ============REAL X VALUES ==================
        nrows = x_data_tmp.shape[0]
        range_freq = x_max_resample - x_min_resample
        x_data_new = np.zeros(nrows, dtype=np.float)
        for i in range(0, nrows):
            x_data_new[i] = range_freq * (i + 1) / npoints_resample + x_min_resample

        y_data_new = funcinterp(x_data_new)
        deltax_new = (np.max(x_data_new) - np.min(x_data_new)) / (x_data_new.shape[0] - 1)

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\n\t\t ======== RESAMPLING THE DATA ======== ({})\n".format(now)
        m += "\t\t Original  data --> Number of points: {0:7d} DeltaX: {1:>10.2f}\n"\
            .format(ndata_orig, deltax_orig)
        m += "\t\t Resampled data --> Number of points: {0:7d} DeltaX: {1:>10.2f}"\
            .format(x_data_new.shape[0], deltax_new)
        print(m) if self._logger is None else self._logger.info(m)

        return np.stack((x_data_new, y_data_new), axis=1)

    # ========================================================================
    def _normalize_data(self, data):

        if self._normalize_method.upper() == "MINMAX":
            norm_method = "Min-Max normalization (mM)"
        elif self._normalize_method.upper() == "VECTORNORM":
            norm_method = "Vector normalization (VN)"
        elif self._normalize_method.upper() == "1-NORM":
            norm_method = "1-norm (1-n)"
        elif self._normalize_method.upper() == "SNV":
            norm_method = "Standard Normal Variate (SNV)"
        else:
            norm_method = None

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\n\t\t ======== NORMALIZING DATA ======== ({})\n".format(now)
        m += "\t\t Normalizing data using the method: {}".format(norm_method)
        print(m) if self._logger is None else self._logger.info(m)

        if self._normalize_method.upper() == "MINMAX":
            data = self.__normalize_minmax_method(data)
        elif self._normalize_method.upper() == "VECTORNORM":
            data = self.__normalize_vectornorm_method(data)
        elif self._normalize_method.upper() == "1-NORM":
            data = self.__normalize_onenorm_method(data)
        elif self._normalize_method.upper() == "SNV":
            data = self.__normalize_snv_method(data)
        else:
            data = None

        return data

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

    # ========================================================================
    def _create_baseline(self, data, baseline_method, params=None):

        ref1 = "von der Esch et al.  J. Chem. Theory Comput. 2021, 17, 985−995"
        ref2 = "https://pybaselines.readthedocs.io/en/latest/api/pybaselines" \
               "/api/index.html#pybaselines.api.Baseline.asls"
        ref3 = "https://pybaselines.readthedocs.io/en/latest/api/pybaselines"\
               "/api/index.html#pybaselines.api.Baseline.mixture_model"
        ref4 = "https://pybaselines.readthedocs.io/en/latest/api/pybaselines" \
               "/api/index.html#pybaselines.api.Baseline.optimize_extended_range"

        if baseline_method.upper() == "2021MUNICH":
            method = "2021Munich. \n\t\t({})".format(ref1)
        elif baseline_method.upper() == "ASLS":
            method = "ASLS. \n\t\t({})".format(ref2)
        elif baseline_method.upper() == "MIXTURE_MODEL":
            method = "MIXTURE_MODEL. \n\t\t({})".format(ref3)
        elif baseline_method.upper() == "OPTIMIZE_EXTENDED_RANGE":
            method = "OPTIMIZE_EXTENDED_RANGE. \n\t\t({})".format(ref4)
        else:
            method = None

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\n\t\t ======== BASELINE CORRECTION ======== ({})\n".format(now)
        m += "\t\t Baseline correction using the method: {}".format(method)
        print(m) if self._logger is None else self._logger.info(m)

        if baseline_method.upper() == "2021MUNICH":
            threshold = float(params[0])
            data = self.__creating_baseline_munich(data, threshold)
        elif baseline_method.upper() == "ASLS":
            lamb, p = params[0:2]
            ydata = asls(data[:, 1], lam=float(lamb), p=float(p))
            data = np.stack((data[:, 0], ydata[0]), axis=1)
            np.savetxt("03-tmp_baseline.dat", data, delimiter=" ")
        elif baseline_method.upper() == "MIXTURE_MODEL":
            lamb, p = params[0:2]
            ydata = mixture_model(data[:, 1], lam=float(lamb), p=float(p))
            data = np.stack((data[:, 0], ydata[0]), axis=1)
            np.savetxt("03-tmp_baseline.dat", data, delimiter=" ")
        elif baseline_method.upper() == "OPTIMIZE_EXTENDED_RANGE":
            ydata = optimize_extended_range(data[:, 1])
            data = np.stack((data[:, 0], ydata[0]), axis=1)
            np.savetxt("03-tmp_baseline.dat", data, delimiter=" ")
        else:
            data = None

        return data

    # ========================================================================
    @staticmethod
    def __creating_baseline_munich(data, threshold, debug=False):

        # Step 1: Differencite the spectrum in data
        dydx = np.gradient(data[:, 1])
        dydx = np.stack((data[:, 0], dydx), axis=1)
        np.savetxt("02-tmp_diff_spectra.dat", dydx, delimiter=" ")

        # Step 2: Get the region of interest using the threshold value
        # If value in diff funnction are greater to threshold then True
        mask_diff_old = np.where(np.abs(dydx[:, 1]) > threshold, True, False)
        mask_diff_new = np.zeros(mask_diff_old.shape, dtype=np.bool)

        # Step 3: Add three neighboring data points in both directions
        for i in range(0, len(mask_diff_old)):
            if mask_diff_old[i]:
                for ineigh in range(-3, 4):
                    idx = i + ineigh
                    if idx < 0 or idx > len(mask_diff_old)-1:
                        continue
                    mask_diff_new[i + ineigh] = 1

        # Step 4: The region of no interest is formed by those marked as False in the mask_diff_new array
        xlist_nointerest_points = list()
        ylist_nointerest_points = list()
        for i in range(0, len(mask_diff_new)):
            if not mask_diff_new[i]:
                xlist_nointerest_points.append(data[i, 0])
                ylist_nointerest_points.append(data[i, 1])

        # Step 5: Interpolate the points to non-interest to create the baseline.
        # The baseline must have the same x-values that the data set
        x_interp = np.linspace(xlist_nointerest_points[0], xlist_nointerest_points[-1], len(data[:, 0]))
        x_interp = data[:, 0]
        y_interp = np.interp(x_interp, xlist_nointerest_points, ylist_nointerest_points)

        if debug:
            points = np.stack((xlist_nointerest_points, ylist_nointerest_points), axis=1)
            np.savetxt("03-tmp_debug_noninterest_points.dat", points, delimiter=" ")

        baseline = np.stack((x_interp, y_interp), axis=1)
        np.savetxt("03-tmp_baseline.dat", baseline, delimiter=" ")

        return baseline

    # =========================================================================
    def _gnuplot_template(self, nx=3, ny=3):

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\n\t\t ======== GNU TEMPLATE ======== ({})\n".format(now)
        m += "\t\t A gnuplot template has been generated:\n"
        m += "\t\t\t To visualize run the following command:\n"
        m += "\t\t\t     gnuplot -p 99-template_expIR.gnuplot"
        print(m) if self._logger is None else self._logger.info(m)

        start = np.min(self._current_data[:, 0])
        end = np.max(self._current_data[:, 0])

        linegnuplot = "reset\n"
        linegnuplot += 'set xlabel "Freq (cm^-1)"\n'
        linegnuplot += 'set ylabel "Intensity"\n'
        linegnuplot += 'set grid\n'
        linegnuplot += 'set xrange [{0:d}:{1:d}]\n'.format(int(start), int(end))
        linegnuplot += 'set yrange [{0:d}:{1:d}]\n'.format(0, 1)
        linegnuplot += 'set format y "%.1f"\n'
        linegnuplot += "# "+50*'*'+"\n"
        linegnuplot += "\n"

        iwxt = 1
        iheight = int(500 * ny)
        iwidth = int(300 * nx)
        linegnuplot += 'set term wxt {0:d} enhanced dashed size {1:d},{2:d} ' \
                       'font "Arial,8"\n'.format(iwxt, iheight, iwidth)
        linegnuplot += 'set multiplot layout {0:d},{1:d}\n'.format(nx, ny)

        linegnuplot += "set title \"Original\"\n"
        linegnuplot += 'p "./{0:s}" u 1:2 w l notitle lc rgb \"black\"\n\n'.format(self._filename)
        if self._normalize_method is not None:
            linegnuplot += "set title \"Normalize {}\"\n".format(self._normalize_method)
            fname_tmp1 = "01_tmp_" + self._pattern + "_normalize_{}".format(self._normalize_method) + ".dat"
            linegnuplot += 'p "./{0:s}" u 1:2 w l notitle lc rgb \"black\" \n\n'.format(fname_tmp1)
        if self._baseline_method is not None:
            linegnuplot += "set title \"Baseline {}\"\n".format(self._baseline_method)
            if self._baseline_method.upper() == "2021MUNICH":
                fname_tmp = "02-tmp_diff_spectra.dat"
                linegnuplot += "set yrange [-0.2:0.2]\n"
                linegnuplot += 'p "./{0:s}" u 1:2 w l notitle lc rgb \"black\", ' \
                               '0.006 w l notitle, -0.006 w l notitle\n\n'.format(fname_tmp)
            if self._baseline_method.upper() == "2021MUNICH" or \
                self._baseline_method.upper() == "ASLS" or \
                    self._baseline_method.upper() == "MIXTURE_MODEL" or \
                    self._baseline_method.upper() == "OPTIMIZE_EXTENDED_RANGE":

                linegnuplot += 'set yrange [{0:f}:{1:f}]\n'.format(0, 1)
                linegnuplot += 'p "./{0:s}" u 1:2 w l notitle lc rgb \"black\", ' \
                               '"./{1:s}" u 1:2 w l ls 6 notitle\n\n'.\
                    format(fname_tmp1, "03-tmp_baseline.dat")

        if self._baseline_method is not None and self._peaksmoothing_method is not None:
            linegnuplot += "set title \"Peak Smoothing {}\"\n".format(self._peaksmoothing_method)
            fname_tmp = "04-tmp_baselinecorrected_{}.dat".format(self._baseline_method)
            linegnuplot += 'p "./{0:s}" u 1:2 w l notitle lc rgb \"black\", ' \
                           '"./{1:s}" u 1:2 w l lc rgb \"blue\" title "Smoothed"\n\n'.\
                           format(fname_tmp, "05-tmp_peaksmoothing_{}.dat".format(self._peaksmoothing_method))
        if self._baseline_method is None and self._peaksmoothing_method is not None:
            linegnuplot += "set title \"Peak Smoothing {}\"\n".format(self._peaksmoothing_method)
            linegnuplot += 'p "./{0:s}" u 1:2 w l lc rgb \"blue\" title "Smoothed"\n\n'. \
                format("05-tmp_peaksmoothing_{}.dat".format(self._peaksmoothing_method))

        linegnuplot += "set title \"Final Spectrum\"\n"
        fname_final = self._pattern + "_corrected.dat"

        linegnuplot += 'p "./{0:s}" u 1:2 w l lc rgb \"blue\" notitle \n\n'. \
                       format(fname_final)

        linegnuplot += "unset multiplot\n"

        linegnuplot += "# "+50*'*'+"\n"
        linegnuplot += "\n"
        iwxt += 1
        linegnuplot += "\n"
        linegnuplot += 'set term wxt {0:d} enhanced dashed size {1:d},{2:d} ' \
                       'font "Arial,8"\n'.format(iwxt, 600, 400)
        linegnuplot += 'set multiplot layout {0:d},{1:d}\n'.format(1, 1)

        linegnuplot += '\nunset yrange\n'
        linegnuplot += 'p\\\n'
        linegnuplot += '    "./{0:s}" u 1:2 w l lc rgb \"black\" title "Original",\\\n'. \
                       format(self._filename)
        linegnuplot += '    "./{0:s}" u 1:2 w l lc rgb \"blue\" title "Corrected"\n\n'. \
                       format(fname_final)

        linegnuplot += "unset multiplot\n"

        linegnuplot += "# "+50*'*'+"\n"
        linegnuplot += "\n"
        iwxt += 1
        linegnuplot += "\n"
        linegnuplot += 'set term wxt {0:d} enhanced dashed size {1:d},{2:d} ' \
                       'font "Arial,8"\n'.format(iwxt, 600, 400)
        linegnuplot += 'set multiplot layout {0:d},{1:d}\n'.format(1, 1)

        linegnuplot += '\nunset yrange\n'
        linegnuplot += 'p\\\n'
        linegnuplot += '    "./{0:s}" u 1:2 w l lc rgb \"blue\" title "Corrected"\n\n'. \
                       format(fname_final)

        with open("99-template_expIR.gnuplot", "w") as fgnuplot:
            fgnuplot.writelines(linegnuplot)
