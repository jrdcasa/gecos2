import os
from collections import defaultdict


def pdbtosdf(filepdbname):

    new_atom_line = []
    new_chg_line = []
    element_list = []
    atom_formal_charges = []

    # Read all ATOM lines in the PDB file.
    with open(filepdbname, 'r') as f:
        idx_atom = 1
        while True:
            line = f.readline()
            if line.find("ATOM") != -1 or line.find("HETATM") != -1 :
                name = line[12:16].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                element = line[76:78].strip()
                element_list.append(element)
                charge = line[78:80].strip()
                if charge:
                    signs = ["+", "-"]
                    if charge[0] in signs:
                        chg = "{0:3s} {1:1s}{2:1s}".format(str(idx_atom), str(charge[0]), str(charge[1]))
                        s = str(charge[0])
                        c = int(charge[1])
                    else:
                        chg = "{0:3s} {1:1s}{2:1s}".format(str(idx_atom), str(charge[1]), str(charge[0]))
                        s = str(charge[1])
                        c = int(charge[0])
                    new_chg_line.append("{}".format(chg))
                    if s == "+":
                        atom_formal_charges.append(c)
                    else:
                        atom_formal_charges.append(-c)
                else:
                    atom_formal_charges.append(0)
                new_atom_line.append("{0:10.4f}{1:10.4f}{2:10.4f} {3:<3s} {4:3s}{5:3s}{6:3s}{7:3s}{8:3s}"
                                     "{9:3s}{10:3s}{11:3s}{12:3s}{13:3s}{14:3s}{15:3s}\n".
                                     format(x, y, z, element,
                                            "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0"))
                idx_atom += 1
            if not line:
                break
    natoms = len(new_atom_line)

    # Read all CONECT lines in the PDB file.
    bond_dict = defaultdict(list)
    with open(filepdbname, 'r') as f:
        while True:
            line = f.readline()
            if line.find("CONECT") != -1:
                spl = line.split()
                iat = str(int(spl[1])-1)
                for ipos in range(2,len(spl)):
                    jat = str(int(spl[ipos])-1)
                    bond_dict[iat].append(jat)
            if not line:
                break

    bond_list = []
    for key, values in bond_dict.items():
        for ival in values:
            element = sorted([int(key), int(ival)])
            if element not in bond_list:
                bond_list.append(element)

    nbonds = len(bond_list)
    bond_orders, bond_orders_resonance = bond_perception(bond_list, element_list, atom_formal_charges)

    # Choose a resonance structure consequent with the formal charges
    idxs_charged = [atom_formal_charges.index(idx) for idx in atom_formal_charges if idx != 0]
    nresonances = len(bond_orders_resonance)
    for iresonance in range(nresonances):
        isfound = False
        idx_ibond = 0
        for ibond in bond_list:
            iat = ibond[0]
            jat = ibond[1]
            if iat in idxs_charged and jat in idxs_charged and bond_orders_resonance[iresonance][idx_ibond] != 1.0:
                isfound = True
            idx_ibond += 1
        if not isfound:
            break

    # write SDF
    path = os.path.split(filepdbname)[0]
    fileext = os.path.split(filepdbname)[1]
    filenamesdf = os.path.splitext(fileext)[0]+"_frompdb.sdf"
    with open(filenamesdf, 'w') as fout:
        fout.writelines("{}".format(filenamesdf))
        fout.writelines("\n\n\n".format(filenamesdf))
        # First column: natoms, Second column: nbonds, Fihth column: 0 no chiral, 1 chiral, Others: 0
        fout.writelines("{0:3d}{1:3d}{2:3d}{3:3d}{4:3d}               999 V2000\n".format(natoms, nbonds, 0, 0, 0))
        # Atoms
        for iline in new_atom_line:
            fout.writelines(iline)
        # Bonds
        idx_bond = 0
        for ibond in bond_list:
            bo = int(bond_orders_resonance[iresonance][idx_bond])
            line = "{0:3d}{1:3d}{2:3d}{3:3d}{4:3d}{5:3d}{6:3d}\n".format(ibond[0]+1, ibond[1]+1, bo, 0, 0, 0, 0)
            fout.writelines(line)
            idx_bond += 1
        # Charges
        for icharge in new_chg_line:
            line = "M  CHG{0:3s}  {1:7s}\n".format("  1", icharge)
            fout.writelines(line)
        fout.writelines("M  END\n")
        fout.writelines("$$$$")


# =========================================================================
def bond_perception(bond_list, elements, formal_charges=None, total_charge=0):

    """
    This function assigns bond orders to the bonds according to the algorithm
    reported in:

    "Automated simultaneous assignment of bond orders and formal charges"
    Ivan D. Welsh and Jane R. Allison
    J. Cheminform (2019) 11:18

    https://doi.org/10.1186/s13321-019-0340-0

    The function uses the external software ``indigo-bondorders`` (located in thirdparty/indigo-bondorder).
    This code is compiled and installed in thirdparty/indigox


    .. warning::
        The structure to assign bonds needs to have all hydrogen bonds correctly placed.
        United atom models do not work with this function.

    """

    import indigox as ix

    # Periodic Table data from indigox
    PT = ix.PeriodicTable()

    # Build a molecule in the indigox framework
    mol = ix.Molecule()
    mol.SetTotalCharge(total_charge)

    # Add all atoms in a dictionary and get the bonds in the
    # framework of indigox program
    all_atoms = dict()
    for i, j in bond_list:
        if i not in all_atoms:
            # Element of i
            e = elements[i]
            all_atoms[i] = mol.NewAtom(PT[e])
            index = all_atoms[i].SetIndex(i)
            name = e + str(index)
            all_atoms[i].SetName(name)
            all_atoms[i].SetFormalCharge(formal_charges[i])

        if j not in all_atoms:
            # Element of j
            e = elements[j]
            all_atoms[j] = mol.NewAtom(PT[e])
            index = all_atoms[j].SetIndex(j)
            name = e+str(index)
            all_atoms[j].SetName(name)
            all_atoms[j].SetFormalCharge(formal_charges[j])

        mol.NewBond(all_atoms[i], all_atoms[j])

    # Setup to use the FPT algorithm with single electrons without preplacing
    # to calculate bond orders and formal charges
    opts = ix.Options.AssignElectrons
    opts.ALGORITHM = opts.Algorithm.FPT
    opts.FPT.ADD_EDGES_TO_TD = False
    opts.FPT.MINIMUM_PROPAGATION_DEPTH = 1
    opts.USE_ELECTRON_PAIRS = False

    # Calculate bond orders and formal charges.
    # Count have the ttotal number of resonance structures
    nresonances = mol.AssignElectrons()
    # print("{} resonace structure(s) calculated "
    #       "with a score of {}.".format(nresonances, mol.GetMinimumElectronAssignmentScore()))

    # Sum all order bonds for the resonace structures.
    nbonds = len(bond_list)
    order_bonds = defaultdict(float)
    order_bonds_resonance_dict = defaultdict(list)
    for iresonance in range(nresonances):

        mol.ApplyElectronAssignment(iresonance)
        index_bond = 0
        for ibond in mol.GetBonds():
            i = ibond.GetSourceAtom().GetIndex()
            j = ibond.GetTargetAtom().GetIndex()
            bo = ibond.GetOrder()
            if bo == bo.SINGLE_BOND:
                order_bonds[index_bond] += 1.
                order_bonds_resonance_dict[iresonance].append(1.)
            elif bo == bo.DOUBLE_BOND:
                order_bonds[index_bond] += 2.
                order_bonds_resonance_dict[iresonance].append(2.)
            elif bo == bo.TRIPLE_BOND:
                order_bonds[index_bond] += 3.
                order_bonds_resonance_dict[iresonance].append(3.)
            else:
                m = "Bond order cannot be assigned between {} and {} atoms".format(i, j)
                m + "Bond order: {}".format(bo)
                print(m)
            index_bond += 1

    # Correct for aromaticity
    for ibond in range(nbonds):
        m = order_bonds[ibond] % nresonances
        if m != 0:
            order_bonds[ibond] = 1.5
        else:
            order_bonds[ibond] /= nresonances

    return order_bonds, order_bonds_resonance_dict

