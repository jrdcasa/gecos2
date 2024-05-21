import os.path
from commonnn import cluster
import sklearn
from sklearn.metrics import pairwise_distances
import matplotlib as mpl
import matplotlib.pyplot as plt


class CnnClusteringWrapper:

    # ==========================================================
    def __init__(self, data=None, logger=None):

        """
        Init
        :param data: An array (N-points, N-dimensions(features))with the data to make clusters
        """

        if data is not None:
            self.setup_data(data)
        else:
            self._data = None

        if logger is not None:
            self._logger = logger
        else:
            self._logger = None

        self._ax_props = None
        self._dot_props = None
        self._plot_defaults()

    # ==========================================================
    def setup_data(self, data):

        """
        Setup the date for clustering
        :param data: An array (N-points, N-dimensions(features))with the data to make clusters

        """

        self._data = data

    # ==========================================================
    def withdata(self, radius_cutoff=2.0, similarity_cutoff=0, isdraw=None):

        """
        The neighbourhood of a specific point will be collected brute-force by computing
        the (Euclidean) distances to all other points and comparing them to the radius cutoff.
        This will be the slowest possible approach but it has fairly conservative memory usage.

        :return:
        """

        list_isdraw_options = ["screen", "png"]
        if isdraw is not None:
            if isdraw not in list_isdraw_options:
                msg = "\n\t\tOptions allowed for draw are {}!!!!\n".format(list_isdraw_options)
                msg += "\t\tYou have passed {} as option!!!!\n".format(isdraw)
                msg += "\t\tClustering is not plotted!!!!\n".format(isdraw)
                print(msg) if self._logger is None else self._logger.warn(msg)
                isdraw = None

        if self._data is None:
            msg = "There is not data for clustering, No clustering is done!!!!"
            print(msg) if self._logger is None else self._logger.warn(msg)
            return None

        clustering = cluster.Clustering(self._data)
        fit = clustering.fit(radius_cutoff=radius_cutoff, similarity_cutoff=similarity_cutoff)
        print(fit) if self._logger is None else self._logger.info(fit)

        if isdraw is not None:
            fig, ax = plt.subplots(1, 2)
            ax[0].set_title("original")
            clustering.evaluate(
                ax=ax[0], original=True,
                ax_props=self._ax_props, plot_props=self._dot_props
                )

            ax[1].set_title("clustered")
            clustering.evaluate(
                ax=ax[1],
                ax_props=self._ax_props, plot_props=self._dot_props
                )
            fig.tight_layout()
            if isdraw == "screen":
                fig.show()
            elif isdraw == "png":
                fig.savefig("figure.png")

        return clustering.labels

    # ==========================================================
    def _plot_defaults(self):

        path = os.path.split(__file__)[0]
        # ======================== PLOT =======================#
        # Matplotlib configuration
        mpl.rc_file(os.path.join(path, "matplotlibrc"), use_default_template=False)

        # Axis property defaults for the plots
        self._ax_props = {
            "aspect": "auto"
        }

        # Property defaults for plotted lines
        self._dot_props = {
            "marker": "o",
            "markeredgecolor": "k"
        }


