# import both numpy and the Cython declarations for numpy
import cython
import numpy as np
cimport numpy as np
from cython.parallel import prange

# declare the interface to the C code
cdef extern from "calc_distance.c":
    #cdef bint USED_OPENMP
    void c_distance_array (double* ref, double* conf, double* dist,
                           double* rijx, double* rijy, double* rijz,
                           int n, int m, int dim)
    void c_distance_diagonal_array (double* ref, double* conf, double* dist,
                                    double* rijx, double* rijy, double* rijz,
                                    int n, int dim)

#OPENMP_ENABLED = True if USED_OPENMP else False

# ====================================================================================
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def calc_minimum_image_neighbors_com(float[:,:] com_coords,
                                     float[:] box, float radius
                                    ):

    cdef int i, j
    cdef float dx, dy, dz, ddsq, radiussq
    cdef int imagex, imagey, imagez
    cdef int nmols = com_coords.shape[0]
    cdef np.ndarray neigh_distances_sq_python
    cdef np.ndarray neigh_images_python
    cdef float[:, :] neigh_distances_sq = np.zeros([nmols, nmols], dtype=np.float32)
    cdef int[:, :] neigh_imagesx = np.zeros([nmols, nmols], dtype=np.int32)
    cdef int[:, :] neigh_imagesy = np.zeros([nmols, nmols], dtype=np.int32)
    cdef int[:, :] neigh_imagesz = np.zeros([nmols, nmols], dtype=np.int32)

    radiussq = radius*radius
    # Calculate all distances between coms.
    # using minimum image convention
    for i in prange(0, nmols, nogil=True):
        for j in range(i+1, nmols):
            imagex = 0
            dx = com_coords[i, 0] - com_coords[j, 0]
            if dx > box[0]*0.5: dx = dx - box[0]; imagex = -1
            if dx <= -box[0]*0.5: dx = dx + box[0]; imagex = +1

            imagey = 0
            dy = com_coords[i, 1] - com_coords[j, 1]
            if dy > box[1]*0.5: dy = dy - box[1]; imagey = -1
            if dy <= -box[1]*0.5: dy = dy + box[1]; imagey = +1

            imagez = 0
            dz = com_coords[i, 2] - com_coords[j, 2]
            if dz > box[2]*0.5: dz = dz - box[2]; imagez = -1
            if dz <= -box[2]*0.5: dz = dz + box[2]; imagez = +1

            ddsq = dx*dx+dy*dy+dz*dz
            neigh_distances_sq[i, j] = ddsq
            neigh_imagesx[i, j] = imagex
            neigh_imagesy[i, j] = imagey
            neigh_imagesz[i, j] = imagez

    pair_distsq_dict = dict()
    pair_image_dict = dict()
    for i in range(0, nmols):
        for j in range(i+1, nmols):
            if neigh_distances_sq[i, j] != 0.0 and neigh_distances_sq[i, j] < radiussq:
                pair_distsq_dict[(i,j)] = neigh_distances_sq[i, j]
                pair_image_dict[(i,j)] = [neigh_imagesx[i, j], neigh_imagesy[i, j], neigh_imagesz[i, j]]

    return pair_distsq_dict, pair_image_dict

# =============================================================================
def calc_distance_array_openmp(np.ndarray[np.float64_t, ndim=2, mode="c"] ref,
                               np.ndarray[np.float64_t, ndim=2, mode="c"] conf,
                               np.ndarray[np.float64_t, ndim=2, mode="c"] dist,
                               np.ndarray[np.float64_t, ndim=2, mode="c"] rijx,
                               np.ndarray[np.float64_t, ndim=2, mode="c"] rijy,
                               np.ndarray[np.float64_t, ndim=2, mode="c"] rijz):

    cdef int rows = ref.shape[0]
    cdef int cols = conf.shape[0]
    dim1 = ref.shape[1]
    dim2 = conf.shape[1]
    if dim1 != dim2:
        print ("ERROR. The dimension of points to calculate distances must be equal!!!!")
        return dist

    c_distance_array(&ref[0,0], &conf[0,0], &dist[0,0], &rijx[0,0],
                     &rijy[0,0], &rijz[0,0], rows, cols, dim1)
    return None

# =============================================================================
def calc_distance_diagonal_openmp(np.ndarray[np.float64_t, ndim=2, mode="c"] ref,
                                  np.ndarray[np.float64_t, ndim=2, mode="c"] conf,
                                  np.ndarray[np.float64_t, ndim=1, mode="c"] dist,
                                  np.ndarray[np.float64_t, ndim=1, mode="c"] rijx,
                                  np.ndarray[np.float64_t, ndim=1, mode="c"] rijy,
                                  np.ndarray[np.float64_t, ndim=1, mode="c"] rijz):

    cdef int rows = ref.shape[0]
    cdef int cols = conf.shape[0]
    dim1 = ref.shape[1]
    dim2 = conf.shape[1]
    if dim1 != dim2:
        print ("ERROR. The dimension of points to calculate distances must be equal!!!!")
        return dist

    c_distance_diagonal_array(&ref[0,0], &conf[0,0], &dist[0], &rijx[0],
                     &rijy[0], &rijz[0], rows, dim1)

    return None
