import psutil
import numpy as np


# #############################################################################
# noinspection PyUnboundLocalVariable
def distance_array(ref, conf, openmp=False):

    """
    Calculate the distances among the points in two arrays using C and Cython

    PBCs are not corrected in the function

    Parameters:
        * ``ref``: (type: ndarray): Reference array containing a number of points with shape [npoints1, dim],
        where dim is the dimension in 2-D or 3-D points
        * ``conf``: (type: ndarray): Reference array containing a number of points with shape [npoints2, dim],
        where dim is the dimension in 2-D or 3-D points
        * ``openmp``: (type: boolean) Parallel or not

    Returns:
        * ``dist``: (type: 2D-ndarray). Distances in Angstroms between the ith atom (row) and the jth atoms
        * ``rijx``: (type: 2D-ndarray). Component-x of the vector i,j
        * ``rijy``: (type: 2D-ndarray). Component-y of the vector i,j
        * ``rijz``: (type: 2D-ndarray). Component-z of the vector i,j

    ``Examples``:


    """

    from ext_libc.c_distances_openmp import calc_minimum_image_neighbours

    rows = ref.shape[0]
    cols = conf.shape[0]

    # Check memory
    vmem_available = psutil.virtual_memory()[1]  # bytes
    # 4 arrays of rows*cols --> float64
    mem_req = 4*rows*cols*8  # bytes
    if 4 * rows * cols * 8 > vmem_available:
        print("There is not sufficient memory to calculate distances...")
        print("Internal_coordinates.py --> distance_array fucntion.")
        print("Memory available: {0:.1f} Mb".format(vmem_available/1024/1024))
        print("Memory required : {0:.1f} Mb".format(mem_req/1024/1024))
        print("Aborting program...")
        exit()

    calc_minimum_image_neighbours()

    return None

# #############################################################################
def dihedral_py(c_at1, c_at2, c_at3, c_at4, units="radians"):

    x1, y1, z1 = c_at1[:]
    x2, y2, z2 = c_at2[:]
    x3, y3, z3 = c_at3[:]
    x4, y4, z4 = c_at4[:]

    #  //1st bond
    delx1 = x1 - x2
    dely1 = y1 - y2
    delz1 = z1 - z2

    # //2nd bond
    delx2 = x3 - x2
    dely2 = y3 - y2
    delz2 = z3 - z2

    delx2m = -delx2
    dely2m = -dely2
    delz2m = -delz2

    # //3rd bond
    delx3 = x4 - x3
    dely3 = y4 - y3
    delz3 = z4 - z3

    # //c,s calculation
    ax = dely1*delz2m - delz1*dely2m
    ay = delz1*delx2m - delx1*delz2m
    az = delx1*dely2m - dely1*delx2m
    bx = dely3*delz2m - delz3*dely2m
    by = delz3*delx2m - delx3*delz2m
    bz = delx3*dely2m - dely3*delx2m

    rasq = ax*ax + ay*ay + az*az
    rbsq = bx*bx + by*by + bz*bz
    rgsq = delx2m*delx2m + dely2m*dely2m + delz2m*delz2m
    rg = np.sqrt(rgsq)

    rginv = ra2inv = rb2inv = 0.0
    if rg > 0:
        rginv = 1.0/rg
    if rasq > 0:
        ra2inv = 1.0/rasq
    if rbsq > 0:
        rb2inv = 1.0/rbsq
    rabinv = np.sqrt(ra2inv*rb2inv)

    c = (ax*bx + ay*by + az*bz)*rabinv
    s = rg*rabinv*(ax*delx3+ay*dely3+az*delz3)

    if c > 1.0:
        c = 1.0
    if c < -1.0:
        c = -1.0

    if units == "degree":
        factor = 57.295779  # // 180.0 / PI;
        ad = factor * np.arctan2(s, c)
    else:
        ad = np.arctan2(s, c)

    return ad
