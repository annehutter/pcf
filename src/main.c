#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> 
#include <assert.h>
#include <unistd.h>
#include <time.h>
#ifdef __MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#ifdef __SPRNG
#include <sprng.h>
#endif
#else
#include <fftw3.h>
#endif

#include "utils.h"
#include "phys_const.h"
#include "cosmology.h"
#include "domain.h"
#include "correlation_function.h"
#include "cross_correlation_function.h"

/* --------------------------------------------------------- */
/* EXMPLE MAIN                                               */
/* --------------------------------------------------------- */

int main (int argc, /*const*/ char * argv[]) { 
    int size = 1;
    int myRank = 0;

    double t1, t2;
    
#ifdef __MPI    
    MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &size); 
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank); 
    
    t1 = MPI_Wtime();
    
    fftw_mpi_init();
#else
    t1 = time(NULL);
#endif
    
    int whichCorrFunc = 10;      // 0 = 3D correlation function, 1 = angular correlation function
    
    int gridsize = 256;         // only requied for domain decomposition according to FFTW (not relevant for correlation function)
    
    /* COSMOLOGY */
    double h = 0.7;
    double omega_m = 0.27;
    double omega_l = 0.73;
    double redshift = 7.;
    
    /* BOX PARAMETERS */
    double boxsize = 109.;
    double xboxsize = 0.6;       // in units of boxsize!
    double yboxsize = 0.6;
    double zboxsize = 0.6;
    double xlow = 0.1;
    double ylow = 0.1;
    double zlow = 0.1;
    int dir = 2;                /* 0 : x-axis; 1 : y-axis; 2 : z-axis */

    /* INITIALIZATION OF STRUCTS */
    domain_t *domain = initDomain(gridsize, myRank, size);
    cosmparams_t *cosm = read_params_to_cosmparams(h, omega_m, omega_l, redshift);
    boxparams_t *box = read_params_to_boxparams(boxsize, xboxsize, yboxsize, zboxsize, xlow, ylow, zlow, dir);
        
    /* GENERATING SOURCES (DATA) */
    int numSources_total = 1000;
    int numSources = numSources_total/size;
    
    double *xsource = allocate_array_double(numSources, "xsource");
    double *ysource = allocate_array_double(numSources, "ysource");
    double *zsource = allocate_array_double(numSources, "zsource");
    
    /* FOR DOMAIN DECOMPOSITION IN GENERATING DATA (need to span between 0 and 1 across all domains) */
    double lowLimit = (double)myRank / (double) size;
    double upLimit = (double)(myRank + 1) / (double) size;
    
#ifdef __SPRNG
    int *stream = init_sprng(1, myRank, size, SEED, SPRNG_DEFAULT);
    for(int source=0; source<numSources; source++)
    {
        xsource[source] = box->xlow + sprng(stream) * (box->xup - box->xlow);
        ysource[source] = box->ylow + sprng(stream) * (box->zup - box->ylow);
        zsource[source] = zlow + (upLimit - lowLimit) * (box->zup - box->zlow) * (myRank + sprng(stream));
    }
    free_sprng(stream);
#else
    srand(time(NULL)*(domain->thisRank+1)*getpid());    
    for(int source=0; source<numSources; source++)
    {
        xsource[source] = box->xlow + randdouble() * (box->xup - box->xlow);
        ysource[source] = box->ylow + randdouble() * (box->zup - box->ylow);
        zsource[source] = zlow + (upLimit - lowLimit) * (box->zup - box->zlow) * (myRank + randdouble());
    }
#endif
    
#ifdef __DEBUG_CORRFUNC
    for(int i=0; i<numSources; i++) printf("rank %d: x,y,z = %e %e %e\n", myRank, xsource[i], ysource[i], zsource[i]);
#endif
    
    /* SELECTING SOURCES IN A SUBBOX */
    boxparams_t *subbox = read_params_to_boxparams(boxsize, xboxsize, yboxsize, zboxsize, xlow, ylow, zlow, dir);
    select_subbox_cartesian(subbox, &numSources, &xsource, &ysource, &zsource);

    /* CALCULATION OF THE CORRELATION FUNCTION */
    if(whichCorrFunc == 0)
        calc_3Dcorrfunc_cartesian(numSources, xsource, ysource, zsource, subbox, domain, "test_3D_corr.dat");
    else if (whichCorrFunc == 1)
        calc_ACF_cartesian(numSources, xsource, ysource, zsource, subbox, cosm, domain, "test_ACF_from_cartesian.dat");
    else if (whichCorrFunc == 10)
        calc_cross_3Dcorrfunc_cartesian(numSources, xsource, ysource, zsource, numSources, xsource, ysource, zsource, subbox, domain, "test_3D_crosscorr.dat");
    else if (whichCorrFunc == 11)
        calc_cross_ACF_cartesian(numSources, xsource, ysource, zsource, numSources, xsource, ysource, zsource, subbox, cosm, domain, "test_crossACF_from_cartesian.dat");
    else
        printf("Other options (e.g. ACF from cosTheta, phi values) are not included in this example main.\n");

    free(xsource);
    free(ysource);
    free(zsource);
    
    deallocate_boxparams(subbox);
    deallocate_boxparams(box);
    deallocate_cosmparams(cosm);
    deallocate_domain(domain);

#ifdef __MPI
    fftw_mpi_cleanup();
        
    t2 = MPI_Wtime();
    printf("Execution took %f s\n", t2-t1);
    MPI_Finalize();
#else
    fftw_cleanup();
    
    t2 = time(NULL);
    printf("Execution took %f s\n", t2-t1);
#endif
    
    return 0;
}


//     int distributeAccordingToDomains = 0;
//     double *xsource = NULL, *ysource = NULL, *zsource = NULL;
//     double *xtmp = NULL, *ytmp = NULL, *ztmp = NULL;
// #ifdef __MPI
//     MPI_Status status;
//     int mpiStat;
// 
//     srand(100);
// 
//     if(myRank == 0)
//     {        
//         xsource = allocate_array_double(numSources_total, "xsource");
//         ysource = allocate_array_double(numSources_total, "ysource");
//         zsource = allocate_array_double(numSources_total, "zsource");
// 
//         for(int source=0; source<numSources_total; source++)
//         {
//             xsource[source] = randdouble();
//             ysource[source] = randdouble();
//             zsource[source] = randdouble();
//         }
//     }
//     
//     
//     if(distributeAccordingToDomains == 0)
//     {
//         if(myRank == 0)
//         {
//             /* sort array */
//             xtmp = allocate_array_double(numSources_total, "xsource");
//             ytmp = allocate_array_double(numSources_total, "ysource");
//             ztmp = allocate_array_double(numSources_total, "zsource");
//             int *rank = allocate_array_int(numSources_total, "rank");
//             int *numSources_rank = allocate_array_int(size, "numSources_rank");
//             int *offset_rank = allocate_array_int(size, "numSources_rank");
//             int *counter = allocate_array_int(size, "counter");
// 
//             for(int source=0; source<numSources_total; source++)
//             {
//                 rank[source] = (int)(zsource[source] * (double)size);
//                 numSources_rank[rank[source]] += 1;
//             }
//             
//             int sum = 0;
//             for(int i=0; i<size; i++)
//             {
//                 offset_rank[i] += sum;
//                 sum += numSources_rank[i];
//             }
//             
//             for(int source=0; source<numSources_total; source++)
//             {
//                 xtmp[offset_rank[rank[source]] + counter[rank[source]]] = xsource[source];
//                 ytmp[offset_rank[rank[source]] + counter[rank[source]]] = ysource[source];
//                 ztmp[offset_rank[rank[source]] + counter[rank[source]]] = zsource[source];
//                 counter[rank[source]]++;
//             }
//             
//             free(xsource);
//             free(ysource);
//             free(zsource);
//             
//             for(int destRank=1; destRank<size; destRank++)
//             {
//                 MPI_Send(&numSources_rank[destRank], 1, MPI_INT, destRank, destRank, MPI_COMM_WORLD);
//                 MPI_Send(&(xtmp[offset_rank[destRank]]), numSources_rank[destRank], MPI_DOUBLE, destRank, destRank*10, MPI_COMM_WORLD);
//                 assert(mpiStat == MPI_SUCCESS);
//                 MPI_Send(&(ytmp[offset_rank[destRank]]), numSources_rank[destRank], MPI_DOUBLE, destRank, destRank*20, MPI_COMM_WORLD);
//                 assert(mpiStat == MPI_SUCCESS);
//                 MPI_Send(&(ztmp[offset_rank[destRank]]), numSources_rank[destRank], MPI_DOUBLE, destRank, destRank*30, MPI_COMM_WORLD);
//                 assert(mpiStat == MPI_SUCCESS);
//             }
//             
//             numSources = numSources_rank[0];
//             xsource = allocate_array_double(numSources, "xsource");
//             ysource = allocate_array_double(numSources, "ysource");
//             zsource = allocate_array_double(numSources, "zsource");
//             
//             for(int source=0; source<numSources; source++)
//             {
//                 xsource[source] = xtmp[source];
//                 ysource[source] = ytmp[source];
//                 zsource[source] = ztmp[source];
//             }
//             
//             free(xtmp);
//             free(ytmp);
//             free(ztmp);
//             
//             free(rank);
//             free(numSources_rank);
//             free(offset_rank);
//             free(counter);
//         }
//         else
//         {
//             MPI_Recv(&numSources, 1, MPI_INT, 0, myRank, MPI_COMM_WORLD, &status);
//             
//             xsource = allocate_array_double(numSources, "xsource");
//             ysource = allocate_array_double(numSources, "ysource");
//             zsource = allocate_array_double(numSources, "zsource");
// 
//             MPI_Recv(xsource, numSources, MPI_DOUBLE, 0, myRank*10, MPI_COMM_WORLD, &status);
//             MPI_Recv(ysource, numSources, MPI_DOUBLE, 0, myRank*20, MPI_COMM_WORLD, &status);
//             MPI_Recv(zsource, numSources, MPI_DOUBLE, 0, myRank*30, MPI_COMM_WORLD, &status);
//         }
//     }
//     else
//     {
//         if(myRank == 0)
//         {
//             for(int destRank=1; destRank<size; destRank++)
//             {
//                 MPI_Send(&(xsource[numSources*destRank]), numSources, MPI_DOUBLE, destRank, destRank*10, MPI_COMM_WORLD);
//                 assert(mpiStat == MPI_SUCCESS);
//                 MPI_Send(&(ysource[numSources*destRank]), numSources, MPI_DOUBLE, destRank, destRank*20, MPI_COMM_WORLD);
//                 assert(mpiStat == MPI_SUCCESS);
//                 MPI_Send(&(zsource[numSources*destRank]), numSources, MPI_DOUBLE, destRank, destRank*30, MPI_COMM_WORLD);
//                 assert(mpiStat == MPI_SUCCESS);
//             }
//     //         xsource = realloc(xsource, numSources);
//     //         ysource = realloc(ysource, numSources);
//     //         zsource = realloc(zsource, numSources);
// 
//     //         xsource = xtmp;
//     //         ysource = ytmp;
//     //         zsource = ztmp;
//         }
//         else
//         {
//             xsource = allocate_array_double(numSources, "xsource");
//             ysource = allocate_array_double(numSources, "xsource");
//             zsource = allocate_array_double(numSources, "xsource");
// 
//             MPI_Recv(xsource, numSources, MPI_DOUBLE, 0, myRank*10, MPI_COMM_WORLD, &status);
//             MPI_Recv(ysource, numSources, MPI_DOUBLE, 0, myRank*20, MPI_COMM_WORLD, &status);
//             MPI_Recv(zsource, numSources, MPI_DOUBLE, 0, myRank*30, MPI_COMM_WORLD, &status);
//         }
//     }
// #endif
