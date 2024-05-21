from collections import defaultdict
import numpy as np
from MDAnalysis.analysis import distances
from MDAnalysis.lib.mdamath import angle as mdangle
from MDAnalysis.lib.mdamath import norm as mddist
from utils.atomic_data import element_vdw_truhlar_radius


class CloseContacts:

    """
    This class needs a data structure as the following:

    info_frag: A dictionary containing the info of fragments
                          required to calculate the VdW interactions.
                        info_frag [0] = [[atomic index, symbol, X, Y, Z, file],
                                         [atomic index, symbol, X, Y, Z, file]
                                         ...
                                         [atomic index, symbol, X, Y, Z, file]]
                                         One sublist for each atom in the fragment
                        ...
                        info_frag[numfrags-1] = ....

    At the moment, the class only works with 1 or 2 fragments. A fragment is defined as
    two diconnected molecules.

    Types of contacts:
        vdwcontacts (based on distance):
            (i,j) is retained if d_ij < vdw(i)+vdw(j)+delta (https://github.com/getcontacts/getcontacts)
        scorecontacts Pycontact criteria:
            (https://www.sciencedirect.com/science/article/pii/S0006349517350518?via%3Dihub)
            (https://www.sciencedirect.com/science/article/pii/S0005273614002946)
            (i,j) is stored if  sij = 1 / [ 1 +exp(5.0*(dij-4.0))] if dij<cutoff; 0 otherwise
    """

    # =========================================================================
    def __init__(self, logger=None):

        """
        Initialize the instance

        :param logger: Logger to write results
        """

        self._logger = logger

        # Dict containing the index of VdW interactions for the QM
        # file matching with the dictionary key
        self._delta = None
        self._vdwcontacts_intermol = defaultdict(list)
        self._vdwcontacts_intramol = defaultdict(list)
        self._scorecontacts_intermol = defaultdict(list)
        self._scorecontacts_intramol = defaultdict(list)
        self._hbond_list = defaultdict(list)
        self._max_cutoff_vdw = None  # Angstroms

    # =========================================================================
    def _calc_distances_fragments(self, info_frag, typeinteract):

        """
        Calculate the distances between atoms belonging to different fragments
        :param info_frag:
        :param typeinteract: "intra" or "inter"
        :return:
        """

        frag_coords = defaultdict(list)
        frag_idx = defaultdict(list)
        frag_symbol = defaultdict(list)

        dist_arr = defaultdict()

        if typeinteract != "intra" and typeinteract != "inter":
            m = '\t\tERROR. {} is not allowed. Allowed values = ["intra", "inter"]'.format(typeinteract)
            print(m) if self._logger is None else self._logger.error(m)
            exit()

        # Compile information
        numfrag = len(info_frag)
        for idx_frag in range(numfrag):
            for idx, iatom in enumerate(info_frag[idx_frag]):
                frag_coords[idx_frag].append(iatom[2])
                frag_idx[idx_frag].append(iatom[0])
                frag_symbol[idx_frag].append(iatom[1])

        # Intermolecular
        if typeinteract == "inter" and numfrag == 2:
            arr1 = np.array(frag_coords[0])
            arr2 = np.array(frag_coords[1])
            dist_arr[0] = distances.distance_array(arr1, arr2)
        else:
            dist_arr[0] = []

        # Intramolecular
        if typeinteract == "intra":
            for idxfrag, ifrag in enumerate(frag_coords):
                arr1 = np.array(frag_coords[idxfrag])
                dist_arr[idxfrag] = distances.distance_array(arr1, arr1)

        return dist_arr

    # =========================================================================
    def _calculate_maxcutoff_vdw(self, info_frag):

        """
        Calculate the vdw cutoff based on maximum vdw radius of the atoms
        in the fragments as:
            vdw_max_1 + vdw_max_2 + self._delta

        :param info_frag:
        :return:
        """

        # Default value
        # if self._delta is None:
        #     self._delta = 0.25

        # Calculate the maximum VdW for the atoms in the selections
        numfrag = len(info_frag)
        max_vdw = numfrag * [0.0]
        for idx_frag in range(numfrag):
            for idx, iatom in enumerate(info_frag[idx_frag]):
                symbol = iatom[1]
                if element_vdw_truhlar_radius[symbol] > max_vdw[idx_frag]:
                    max_vdw[idx_frag] = element_vdw_truhlar_radius[symbol]

        if len(max_vdw) == 1:
            self._max_cutoff_vdw = max_vdw[0] + max_vdw[0] + self._delta
        else:
            self._max_cutoff_vdw = max_vdw[0] + max_vdw[1] + self._delta

    # =========================================================================
    @staticmethod
    def score(dij):
        """
        Calculate score as (https://www.sciencedirect.com/science/article/pii/S0006349517350518?via%3Dihub)

        :param dij:
        :return:
        """

        return 1.0 / (1 + np.exp(5.0 * (dij - 4.0)))

    # =========================================================================
    def vanderwaals(self, id_label, info_frag, list_nb_pairs, delta=0.25, noignore_hs=False):

        """

        A VdW interaction is defined as:
                r_vdw(at1) + r_wdw(at2) + delta

        :param info_frag: A dictionary containing the info of fragments
                          required to calculate the VdW interactions.
                        info_frag [0] = [[atomic index, symbol, X, Y, Z, file],
                                         [atomic index, symbol, X, Y, Z, file]
                                         ...
                                         [atomic index, symbol, X, Y, Z, file]]
                                         One sublist for each atom in the fragment
                        ...
                        info_frag[numfrags-1] = ....
        :param delta:
        :param list_nb_pairs
        :param id_label
        :param noignore_hs
        :return:
        """

        # Delta value
        self._delta = delta

        # Check for only 1 or two fragments
        numfrags = len(info_frag)
        if numfrags != 2 and numfrags != 1:
            m = "\t\tERROR. VdW interactions can only be calculated for\n " \
                "\t\tone or two fragments in the current implementation\n " \
                "\t\t({} fragments found)".format(numfrags)
            print(m) if self._logger is None else self._logger.error(m)
            exit()

        # Get cutoff_vdw if is None
        if self._max_cutoff_vdw is None:
            self._calculate_maxcutoff_vdw(info_frag)

        # Calculate the distances and get the contact pairs for vdw (intermolecular)
        if numfrags == 2:
            distance_arr_inter = self._calc_distances_fragments(info_frag, typeinteract="inter")
            x_axis, y_axis = distance_arr_inter[0].shape
            for ilocal in range(0, x_axis):
                for jlocal in range(0, y_axis):
                    symbol1 = info_frag[0][ilocal][1]
                    symbol2 = info_frag[1][jlocal][1]
                    idx_i = info_frag[0][ilocal][0]
                    idx_j = info_frag[1][jlocal][0]
                    if not noignore_hs and (symbol1 == "H" or symbol2 == "H"):
                        continue
                    dij = distance_arr_inter[0][ilocal, jlocal]
                    if dij < self._max_cutoff_vdw:
                        vdw1 = element_vdw_truhlar_radius[symbol1]
                        vdw2 = element_vdw_truhlar_radius[symbol2]
                        cutoff_vdw = vdw1 + vdw2 + self._delta
                        if dij < cutoff_vdw:
                            # Debug
                            # print(info_frag[0][ilocal][0], info_frag[1][jlocal][0],
                            #       distance_arr_inter[ilocal, jlocal], cutoff_vdw)
                            self._vdwcontacts_intermol[id_label].append([idx_i, idx_j])
                    if dij < 5.0:
                        self._scorecontacts_intermol[id_label].append([idx_i, idx_j, self.score(dij)])

            # Calculate the distances and get the contact pairs for vdw (intramolecular)
            distance_arr_intra = self._calc_distances_fragments(info_frag, typeinteract="intra")
            for idx_ifrag, dist_frag in enumerate(distance_arr_intra):
                x_axis, y_axis = distance_arr_intra[idx_ifrag].shape
                for ilocal in range(0, x_axis):
                    for jlocal in range(ilocal, y_axis):
                        symbol1 = info_frag[idx_ifrag][ilocal][1]
                        symbol2 = info_frag[idx_ifrag][jlocal][1]
                        idx_i = info_frag[idx_ifrag][ilocal][0]
                        idx_j = info_frag[idx_ifrag][jlocal][0]
                        if not noignore_hs and (symbol1 == "H" or symbol2 == "H"):
                            continue
                        if idx_i == idx_j:
                            continue
                        # Exclude bonded pair of atoms
                        if not {idx_i, idx_j} in list_nb_pairs[idx_ifrag]:
                            continue
                        dij = distance_arr_intra[idx_ifrag][ilocal, jlocal]
                        if dij < self._max_cutoff_vdw:
                            vdw1 = element_vdw_truhlar_radius[symbol1]
                            vdw2 = element_vdw_truhlar_radius[symbol2]
                            cutoff_vdw = vdw1 + vdw2 + self._delta
                            if dij < cutoff_vdw:
                                # Debug
                                # print(info_frag[idx_ifrag][ilocal][0], info_frag[idx_ifrag][jlocal][0],
                                #       distance_arr_intra[ilocal, jlocal], cutoff_vdw)
                                self._vdwcontacts_intramol[id_label].append([idx_i, idx_j])
                        if dij < 5.0:
                            self._scorecontacts_intramol[id_label].append([idx_i, idx_j, self.score(dij)])

        elif numfrags == 1:
            distance_arr_intra = self._calc_distances_fragments(info_frag, typeinteract="intra")
            for idx_ifrag, dist_frag in enumerate(distance_arr_intra):
                x_axis, y_axis = distance_arr_intra[idx_ifrag].shape
                for ilocal in range(0, x_axis):
                    for jlocal in range(ilocal, y_axis):
                        symbol1 = info_frag[idx_ifrag][ilocal][1]
                        symbol2 = info_frag[idx_ifrag][jlocal][1]
                        idx_i = info_frag[idx_ifrag][ilocal][0]
                        idx_j = info_frag[idx_ifrag][jlocal][0]
                        if not noignore_hs and (symbol1 == "H" or symbol2 == "H"):
                            continue
                        if idx_i == idx_j:
                            continue
                        # Exclude bonded pair of atoms
                        if not {idx_i, idx_j} in list_nb_pairs[idx_ifrag]:
                            continue
                        dij = distance_arr_intra[idx_ifrag][ilocal, jlocal]
                        if dij < self._max_cutoff_vdw:
                            vdw1 = element_vdw_truhlar_radius[symbol1]
                            vdw2 = element_vdw_truhlar_radius[symbol2]
                            cutoff_vdw = vdw1 + vdw2 + self._delta
                            if dij < cutoff_vdw:
                                # Debug
                                # print(info_frag[idx_ifrag][ilocal][0], info_frag[idx_ifrag][jlocal][0],
                                #       distance_arr_intra[ilocal, jlocal], cutoff_vdw)
                                self._vdwcontacts_intramol[id_label].append([idx_i, idx_j])
                        if dij < 5.0:
                            self._scorecontacts_intramol[id_label].append([idx_i, idx_j, self.score(dij)])
        # DEBUG
        # print(self._vdwcontacts_intermol)
        # print(self._scorecontacts_intermol)
        # print(self._vdwcontacts_intramol)
        # print(self._scorecontacts_intramol)

    # =========================================================================
    def hbonds(self, id_label, info_frag, coords,
               hd_atom_idx_list, hd_d_idx_dict,
               acc_atom_idx_list, hbond_dist_cutoff, hbond_angle_cutoff):

        """

        Calculate the hydrogen bonds

        Dn--H .... :Ac

        Jeffrey[2] categorizes H bonds with donor-acceptor
        distances of 2.2-2.5 Å as “strong, mostly covalent”,
        2.5-3.2 Å as “moderate, mostly electrostatic”,
        3.2-4.0 Å as “weak, electrostatic” (page 12).
        Energies are given as 40-14, 15-4, and <4 kcal/mol respectively.
        (Jeffrey, George A.; An introduction to hydrogen bonding, Oxford University Press, 1997.)

        Interactions with D–H⋯A angle in the range 120–140° are seen to have substantially
        reduced stabilisation energies and angles below 120° are generally unlikely to correspond
        to significant interactions (CrystEngComm, 2009,11, 1563-1571).

        :param info_frag: A dictionary containing the info of fragments
                          required to calculate the VdW interactions.
                        info_frag [0] = [[atomic index, symbol, X, Y, Z, file],
                                         [atomic index, symbol, X, Y, Z, file]
                                         ...
                                         [atomic index, symbol, X, Y, Z, file]]
                                         One sublist for each atom in the fragment
                        ...
                        info_frag[numfrags-1] = ....
        :param hd_atom_idx_list
        :param hd_d_idx_dict
        :param coords
        :param acc_atom_idx_list
        :param hbond_dist_cutoff
        :param hbond_angle_cutoff
        :param id_label
        :return:
        """

        # Check for only 1 or two fragments
        numfrags = len(info_frag)
        if numfrags != 2 and numfrags != 1:
            m = "\t\tERROR. VdW interactions can only be calculated for\n " \
                "\t\tone or two fragments in the current implementation\n " \
                "\t\t({} fragments found)".format(numfrags)
            print(m) if self._logger is None else self._logger.error(m)
            exit()

        # Check distance and angle between HD and A atoms and make list of H
        # bonds (atoms involved, angle and distance)
        for hd_atom_idx in hd_atom_idx_list:
            d_atom_idx = hd_d_idx_dict[hd_atom_idx]
            d_atom_coords = coords[d_atom_idx]
            hd_atom_coords = coords[hd_atom_idx]
            for acc_atom_idx in acc_atom_idx_list:
                if acc_atom_idx == d_atom_idx:        # don't consider the H atom's D atom as A atom
                    continue
                else:
                    acc_atom_coords = coords[acc_atom_idx]

                v1 = np.array(d_atom_coords-hd_atom_coords)
                v2 = np.array(hd_atom_coords-acc_atom_coords)
                hbond_dist = mddist(v2)
                hbond_angle = 180.-mdangle(v1, v2)*180/np.pi
                if hbond_angle > hbond_angle_cutoff and hbond_dist < hbond_dist_cutoff:
                    self._hbond_list[id_label].append([d_atom_idx, hd_atom_idx, acc_atom_idx, hbond_dist, hbond_angle])

        # DEBUG
        # print(self._hbond_list)


