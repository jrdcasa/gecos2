#ifndef __DISTANCES_H
#define __DISTANCES_H

#include <stdio.h>
#include <math.h>

#ifdef PARALLEL
  #include <omp.h>
#endif


static void c_distance_array(double* ref,
                             double* conf,
                             double* dist,
                             double* rijx,
                             double* rijy,
                             double* rijz,
                             int row,
                             int col,
                             int dim)
{

    // It calculates the distance between all points in the ref array and conf array
    //
    //
    //  Example:
    //   ref=np.ndarray(shape=(2,3), buffer=np.array([[1.0, 1.1, 1.2],
    //                                                [2.1, 2.2, 2.3]]))
    //   conf=np.ndarray(shape=(4,3), buffer=np.array([[11.0, 11.1, 11.2],
    //                                                 [22.0, 22.1, 22.2],
    //                                                 [33.0, 33.1, 33.2],
    //                                                 [44.0, 44.1, 44.2]]))
    //  Result:
    //   dist = [[17.32050808 36.37306696 55.42562584 74.47818473]
    //          [15.41525219 34.46781107 53.52036995 72.57292884]]
    //
    // ref array must have the shape [row,3]. Row is the number of points in ref
    // conf array must have the shape [col,3]. Col is the number of points in conf
    // dim is the dimension of the point.

    double x, y, z;
    x = y = z = 0.0;
    int r, c;

    // row : Number of points in the array ref
    // col : Number of points in the array conf

    #ifdef PARALLEL
      int procs = omp_get_num_procs();
      omp_set_num_threads(procs);
    #endif

    #ifdef PARALLEL
       #pragma omp parallel for private(r, c, x, y, z) shared(dist, rijx, rijy, rijz) collapse(2)
    #endif
    for (r=0; r<row; r++) {
        for (c=0; c<col; c++) {
            //#ifdef PARALLEL
            //    printf("ithread: %d \n", omp_get_thread_num());
            //#endif
            x = conf[c*dim+0]-ref[r*dim+0];
            y = conf[c*dim+1]-ref[r*dim+1];
            z = conf[c*dim+2]-ref[r*dim+2];
            rijx[r*col+c] = x;
            rijy[r*col+c] = y;
            rijz[r*col+c] = z;
            dist[r*col+c] = sqrt((x*x)+(y*y)+(z*z));
        }
    }

}

static void c_distance_diagonal_array(double* ref,
                                      double* conf,
                                      double* dist,
                                      double* rijx,
                                      double* rijy,
                                      double* rijz,
                                      int row,
                                      int dim)
{

    // It calculates the distance in the diagonal between points in the ref array and conf array
    //
    //
    //  Example:
    //   ref=np.ndarray(shape=(2,3), buffer=np.array([[1.0, 1.1, 1.2],
    //                                                [2.1, 2.2, 2.3]]))
    //   conf=np.ndarray(shape=(4,3), buffer=np.array([[11.0, 11.1, 11.2],
    //                                                 [22.0, 22.1, 22.2],
    //                                                 [33.0, 33.1, 33.2],
    //                                                 [44.0, 44.1, 44.2]]))
    //  Result:
    //   dist = [[17.32050808     0.0      0.0 0.0]
    //          [    0.0     34.46781107   0.0 0.0]]
    //
    // ref array must have the shape [row,3]. Row is the number of points in ref
    // conf array must have the shape [col,3]. Col is the number of points in conf
    // dim is the dimension of the point.

    double x, y, z;
    x = y = z = 0.0;
    int r;
    // row : Number of points in the array ref
    // col : Number of points in the array conf

    #ifdef PARALLEL
        int procs = omp_get_num_procs();
        omp_set_num_threads(procs);
    #endif

    #ifdef PARALLEL
        #pragma omp parallel for private(r, x, y, z) shared(dist, rijx, rijy, rijz)
    #endif
    for (r=0; r<row; r++) {
//            #ifdef PARALLEL
//            printf("ithread: %d \n", omp_get_thread_num());
//            #endif
        x = conf[r*dim+0]-ref[r*dim+0];
        y = conf[r*dim+1]-ref[r*dim+1];
        z = conf[r*dim+2]-ref[r*dim+2];
        rijx[r] = x;
        rijy[r] = y;
        rijz[r] = z;
        dist[r] = sqrt((x*x)+(y*y)+(z*z));
    }

}

#endif