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

#define PI 3.14159265358979323846
#define SEED 985456376
#define SEED2 0


cosmparams_t *initCosmparams()
{
    cosmparams_t *newCosmparams = malloc(sizeof(cosmparams_t));
    if(newCosmparams == NULL)
    {
        fprintf(stderr, "ERROR: initCosmparams: Not enough memory to allocate newCosmparams.\n");
        exit(EXIT_FAILURE);
    }
    
    newCosmparams->h = 0.;
    newCosmparams->omega_m = 0.;
    newCosmparams->omega_l = 0.;
    newCosmparams->redshift = 0.;
    
    return newCosmparams;
}

cosmparams_t *read_params_to_cosmparams(double h, double omega_m , double omega_l, double redshift)
{
    cosmparams_t * newCosmparams = initCosmparams();
    
    newCosmparams->h = h;
    newCosmparams->omega_m = omega_m;
    newCosmparams->omega_l = omega_l;
    newCosmparams->redshift = redshift;
    
    return newCosmparams;
}

void deallocate_cosmparams(cosmparams_t *thisCosmparams)
{
    if(thisCosmparams != NULL) free(thisCosmparams);
}

boxparams_t *initBoxparams()
{
    boxparams_t *newBoxparams = malloc(sizeof(boxparams_t));
    if(newBoxparams == NULL)
    {   
        fprintf(stderr, "ERROR: initBoxparams: Not enough memory to allocate newBoxparams.\n");
        exit(EXIT_FAILURE);
    }
    
    newBoxparams->boxsize = 1.;
    
    newBoxparams->xboxsize = 1.;
    newBoxparams->yboxsize = 1.;
    newBoxparams->zboxsize = 1.;
    
    newBoxparams->xlow = 0.;
    newBoxparams->ylow = 0.;
    newBoxparams->zlow = 0.;
    
    newBoxparams->xup = 1.;
    newBoxparams->yup = 1.;
    newBoxparams->zup = 1.;
    
    newBoxparams->xc = 0.5;
    newBoxparams->yc = 0.5;
    newBoxparams->zc = 0.5;

    newBoxparams->dir = 0;
    
    return newBoxparams;
}

boxparams_t *read_params_to_boxparams(double boxsize, double xboxsize, double yboxsize, double zboxsize, double xlow, double ylow, double zlow, int dir)
{
    boxparams_t *newBoxparams = initBoxparams();
    
    newBoxparams->boxsize = boxsize;    // h^-1 Mpc

    /* all following values are in units of boxsize */
    newBoxparams->xboxsize = xboxsize;  
    newBoxparams->yboxsize = yboxsize;
    newBoxparams->zboxsize = zboxsize;
    
    newBoxparams->xlow = xlow;
    newBoxparams->ylow = ylow;
    newBoxparams->zlow = zlow;
    
    newBoxparams->xup = xlow + xboxsize;
    newBoxparams->yup = ylow + yboxsize;
    newBoxparams->zup = zlow + zboxsize;
    
    newBoxparams->xc = xlow + 0.5 * xboxsize;
    newBoxparams->yc = ylow + 0.5 * yboxsize;
    newBoxparams->zc = zlow + 0.5 * zboxsize;

    newBoxparams->dir = dir;
    
    return newBoxparams;
}

void deallocate_boxparams(boxparams_t *thisBoxparams)
{
    if(thisBoxparams != NULL) free(thisBoxparams);
}

/* --------------------------------------------------------- */
/* SELECT SUBBOX IN SIMULATION BOX                           */
/* --------------------------------------------------------- */

void select_subbox_cartesian(boxparams_t *thisBoxparams, int *numSources, double **xsource, double **ysource, double **zsource)
{
    int tmpNumSources = *numSources;
    double *tmpXsource = allocate_array_double(*numSources, "tmpXsource");
    double *tmpYsource = allocate_array_double(*numSources, "tmpYsource");
    double *tmpZsource = allocate_array_double(*numSources, "tmpZsource");
    
    double xlow = thisBoxparams->xlow;
    double ylow = thisBoxparams->ylow;
    double zlow = thisBoxparams->zlow;
    
    double xup = thisBoxparams->xup;
    double yup = thisBoxparams->yup;
    double zup = thisBoxparams->zup;
    
    int counter = 0;
    for(int source=0; source<tmpNumSources; source++)
    {
        if(xlow <= (*xsource)[source] && (*xsource)[source] < xup && 
           ylow <= (*ysource)[source] && (*ysource)[source] < yup && 
           zlow <= (*zsource)[source] && (*zsource)[source] < zup)
        {
            tmpXsource[counter] = (*xsource)[source];
            tmpYsource[counter] = (*ysource)[source];
            tmpZsource[counter] = (*zsource)[source];
            counter++;
        }
    }
    
#ifdef __DEBUG_CORRFUNC
    printf("Selected %d from %d sources\n", counter, *numSources);
#endif
    
    free(*xsource);
    free(*ysource);
    free(*zsource);
    
    *numSources = counter;
    *xsource = realloc(tmpXsource, (*numSources) * sizeof(double));
    *ysource = realloc(tmpYsource, (*numSources) * sizeof(double));
    *zsource = realloc(tmpZsource, (*numSources) * sizeof(double));
}

/* --------------------------------------------------------- */
/* GENERATING RANDOMS                                        */
/* --------------------------------------------------------- */

/* RANDOMS FOR 3D GRID: x, y, z */

void generate_randoms_cartesian(domain_t *thisDomain, boxparams_t *thisBoxparams, int numSources, int *numRand, double **xrand, double **yrand, double **zrand)
{
    int thisRank = thisDomain->thisRank;
    int recvNumSources = 0;
    
#ifdef __MPI
    MPI_Allreduce(&numSources, &recvNumSources, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    *numRand = 10 * recvNumSources / thisDomain->size;
#else
    *numRand = 10 * numSources;
#endif
    
    double *xrand_tmp = allocate_array_double(*numRand, "xrand");
    double *yrand_tmp = allocate_array_double(*numRand, "yrand");
    double *zrand_tmp = allocate_array_double(*numRand, "zrand");
    
    double lowLimit = (double)(thisDomain->local_0_start) / (double)(thisDomain->nbins);
    double upLimit = (double)(thisDomain->local_0_start + thisDomain->local_n0) / (double)(thisDomain->nbins);
    
    double xlow = thisBoxparams->xlow;
    double ylow = thisBoxparams->ylow;
    double zlow = thisBoxparams->zlow;
    double xup = thisBoxparams->xup;
    double yup = thisBoxparams->yup;
    double zup = thisBoxparams->zup;
    
#ifdef __SPRNG
    int *stream = init_sprng(4, thisDomain->thisRank, thisDomain->size, SEED2, SPRNG_DEFAULT);
    for(int source=0; source<*numRand; source++)
    {
        xrand_tmp[source] = xlow + sprng(stream) * (xup - xlow);
        yrand_tmp[source] = ylow + sprng(stream) * (yup - ylow);
        zrand_tmp[source] = zlow + (upLimit - lowLimit) * (zup - zlow) * (thisRank + sprng(stream));
    }
    free_sprng(stream);
#else
    srand(time(NULL)*(thisDomain->thisRank + 1)*getpid());
//     srand(186389 * (thisDomain->thisRank + 1));
    
    for(int source=0; source<*numRand; source++)
    {
        xrand_tmp[source] = xlow + randdouble() * (xup - xlow);
        yrand_tmp[source] = ylow + randdouble() * (yup - ylow);
        zrand_tmp[source] = zlow + (upLimit - lowLimit) * (zup - zlow) * (thisRank + randdouble());
    }
#endif
    
#ifdef __DEBUG_CORRFUNC
    printf("rank %d: numRand = %d\n", thisRank, *numRand);
    for(int i=0; i<*numRand; i++) printf("rank %d: random: x,y,z = %e %e %e\n", thisRank, xrand_tmp[i], yrand_tmp[i], zrand_tmp[i]);
#endif

    *xrand = xrand_tmp;
    *yrand = yrand_tmp;
    *zrand = zrand_tmp;
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void deallocate_randoms_cartesian(double *xrand, double *yrand, double *zrand)
{
    if(xrand != NULL) free(xrand);
    if(yrand != NULL) free(yrand);
    if(zrand != NULL) free(zrand);
}

/* RANDOMS FOR ACF: cos(theta), phi */

void generate_randoms_ra_dec(domain_t *thisDomain, cosmparams_t *thisCosmparams, boxparams_t *thisBoxparams, int numSources, int *numRand, double **rRand, double **cosThetaRand, double **phiRand)
{    
    double *xrand = NULL;
    double *yrand = NULL;
    double *zrand = NULL;
    
    generate_randoms_cartesian(thisDomain, thisBoxparams, numSources, numRand, &xrand, &yrand, &zrand);

    *rRand = NULL;
    *cosThetaRand = NULL;
    *phiRand = NULL;
    
    generate_ra_dec(thisCosmparams, thisBoxparams, *numRand, xrand, yrand, zrand, rRand, cosThetaRand, phiRand);

    free(xrand);
    free(yrand);
    free(zrand);
}

void deallocate_randoms_ra_dec(double *rRand, double *cosThetaRand, double *phiRand)
{
    if(rRand != NULL) free(rRand);
    if(cosThetaRand != NULL) free(cosThetaRand);
    if(phiRand != NULL) free(phiRand);
}

/* --------------------------------------------------------- */
/* GENERATING COSTHETA, PHI FROM X,Y,Z                       */
/* --------------------------------------------------------- */

void generate_ra_dec(cosmparams_t *thisCosmparams, boxparams_t *thisBoxparams, int numSources, double *xsource, double *ysource, double *zsource, double **r, double **cosTheta, double **phi)
{
    double comDistInMpc = calcComDist(thisCosmparams->h, thisCosmparams->omega_m, thisCosmparams->omega_l, thisCosmparams->redshift, 0.)/Mpc_cm;
    double comDistInBoxsize = comDistInMpc / thisBoxparams->boxsize * thisCosmparams->h;
    
    double xc = thisBoxparams->xc;
    double yc = thisBoxparams->yc;
    double zc = thisBoxparams->zc;

    double *r_tmp = allocate_array_double(numSources, "r");
    double *cosTheta_tmp = allocate_array_double(numSources, "cosTheta");
    double *phi_tmp = allocate_array_double(numSources, "phi");
    
    if(thisBoxparams->dir == 0)
        xc = xc - comDistInBoxsize;
    else if(thisBoxparams->dir == 1)
        yc = yc - comDistInBoxsize;
    else
        zc = zc -comDistInBoxsize;
    
    int counter = 0;
    for(int source=0; source<numSources; source++)
    {
        r_tmp[counter] = sqrt((xsource[source]-xc) * (xsource[source]-xc) + (ysource[source]-yc) * (ysource[source]-yc) + (zsource[source]-zc) * (zsource[source]-zc));
        cosTheta_tmp[counter] = (zsource[source]-zc) / r_tmp[counter];
        phi_tmp[counter] = atan2(ysource[source]-yc, xsource[source]-xc);
        
        if(phi_tmp[counter] < 0.)
            phi_tmp[counter] = phi_tmp[counter] + 2.*PI;
        
        if(r_tmp[counter] == 0.)
        {
            cosTheta_tmp[counter] = 0.;
            phi_tmp[counter] = 0.;
        }
                    
        counter++;
    }
    
    assert(numSources == counter);
    *r = r_tmp;
    *cosTheta = cosTheta_tmp;
    *phi = phi_tmp;
}

void deallocate_ra_dec(double *r, double *cosTheta, double *phi)
{
    if(r != NULL) free(r);
    if(cosTheta != NULL) free(cosTheta);
    if(phi != NULL) free(phi);
}

/* --------------------------------------------------------- */
/* MPI COMMUNICATION ROUTINES                                */
/* --------------------------------------------------------- */

#ifdef __MPI
void send_recv_array(int thisRank, int destRank, int num, double *array, int sendingRank, int *recvNum, double **recvArray)
{
    MPI_Request send_request;
    MPI_Status status;
    int mpiStat;
    
    int tmp_recvNum = *recvNum;
    
    /* sending & receiving number of entries that are going to be sent */
    mpiStat = MPI_Isend(&num, 1, MPI_INT, destRank, thisRank*10, MPI_COMM_WORLD, &send_request);
    assert(mpiStat == MPI_SUCCESS);
    mpiStat = MPI_Recv(&tmp_recvNum, 1, MPI_INT, sendingRank, sendingRank*10, MPI_COMM_WORLD, &status);
    assert(mpiStat == MPI_SUCCESS);
        
    /* allocating memory for receiving array */
    *recvArray = allocate_array_double(tmp_recvNum, "array");
        
    /* sending & receiving array */
    mpiStat = MPI_Isend(array, num, MPI_DOUBLE, destRank, thisRank*20, MPI_COMM_WORLD, &send_request);
    assert(mpiStat == MPI_SUCCESS);
    mpiStat = MPI_Recv(*recvArray, tmp_recvNum, MPI_DOUBLE, sendingRank, sendingRank*20, MPI_COMM_WORLD, &status);
    assert(mpiStat == MPI_SUCCESS);
        
    *recvNum = tmp_recvNum;
}
#endif

#ifdef __MPI
void send_recv_array_separated(int thisRank, int destRank, int num, double *array, int sendingRank, int *recvNum, double **recvArray)
{
    MPI_Status status;
    int mpiStat;
    
    int tmp_recvNum = *recvNum;
    
    if(thisRank < destRank)
    {
        /* sending number of entries that are going to be sent */
        mpiStat = MPI_Send(&num, 1, MPI_INT, destRank, thisRank*10, MPI_COMM_WORLD);
        assert(mpiStat == MPI_SUCCESS);
        
        /* sending array */
        mpiStat = MPI_Send(array, num, MPI_DOUBLE, destRank, thisRank*20, MPI_COMM_WORLD);
        assert(mpiStat == MPI_SUCCESS);
    }
    
    if(thisRank > sendingRank)
    {
        /* receive number of entries that are going to be received */
        mpiStat = MPI_Recv(&tmp_recvNum, 1, MPI_INT, sendingRank, sendingRank*10, MPI_COMM_WORLD, &status);
        assert(mpiStat == MPI_SUCCESS);
    
        /* allocating memory for receiving array */
        *recvArray = allocate_array_double(tmp_recvNum, "array");
        
        /* receving array */
        mpiStat = MPI_Recv(*recvArray, tmp_recvNum, MPI_DOUBLE, sendingRank, sendingRank*20, MPI_COMM_WORLD, &status);
        assert(mpiStat == MPI_SUCCESS);
        
        *recvNum = tmp_recvNum;
    }
}
#endif

/* --------------------------------------------------------- */
/* COMPUTING DISTANCE HISTOGRAM                              */
/* --------------------------------------------------------- */

/* FOR 3D GRIDS */
double measure_distance_3D(double x1, double y1, double z1, double x2, double y2, double z2)
{
    return sqrt((x1-x2) * (x1-x2) + (y1-y2) * (y1-y2) + (z1-z2) * (z1-z2));
}

int get_distance_r_bin(double distance, int distanceInLog, int numBins)
{
    // good values for numBins: 10 for log and 30 for linear scaling

    int distanceBin = 0;
    if(distanceInLog != 1)
    {
        distanceBin = (int)(0.5 * distance * numBins);
    }
    else
    {
        distanceBin = (int)(0.5 * log10(distance) / 3. * (double)numBins);
    }
    
    return distanceBin;
}

void calc_distance_distribution_DD_r(int num, double *xsource, double *ysource, double *zsource, int numBins, long int **hist, int distanceInLog)
{
    long int *tmpHist = *hist;
    double distance = 0.;
    int distanceBin = 0;
    
    for(int source1=0; source1<num; source1++)
    {
        for(int source2=source1+1; source2<num; source2++)
        {
            distance = measure_distance_3D(xsource[source1], ysource[source1], zsource[source1], xsource[source2], ysource[source2], zsource[source2]);
            distanceBin = get_distance_r_bin(distance, distanceInLog, numBins);
            
            if(distanceBin >=0 && distanceBin<numBins)
                tmpHist[distanceBin] = tmpHist[distanceBin] + 2;
            else
                printf("distanceBin = %d\n", distanceBin);
        }
    }
}

void calc_distance_distribution_DR_r(int num1, double *xsource1, double *ysource1, double *zsource1, int num2, double *xsource2, double *ysource2, double *zsource2, int numBins, long int **hist, int distanceInLog)
{
    long int *tmpHist = *hist;
    double distance = 0.;
    int distanceBin = 0;
    
    for(int source1=0; source1<num1; source1++)
    {
        for(int source2=0; source2<num2; source2++)
        {
            distance = measure_distance_3D(xsource1[source1], ysource1[source1], zsource1[source1], xsource2[source2], ysource2[source2], zsource2[source2]);
            distanceBin = get_distance_r_bin(distance, distanceInLog, numBins);

            if(distanceBin >=0 && distanceBin<numBins)
                tmpHist[distanceBin]++;
        }
    }
}

void calc_distance_distribution_DD_r_over_ranks(int num1, double *xsource1, double *ysource1, double *zsource1, int num2, double *xsource2, double *ysource2, double *zsource2, int numBins, long int **hist, int distanceInLog)
{
    long int *tmpHist = *hist;
    double distance = 0.;
    int distanceBin = 0;
    
    for(int source1=0; source1<num1; source1++)
    {
        for(int source2=0; source2<num2; source2++)
        {
            distance = measure_distance_3D(xsource1[source1], ysource1[source1], zsource1[source1], xsource2[source2], ysource2[source2], zsource2[source2]);
            distanceBin = get_distance_r_bin(distance, distanceInLog, numBins);

            if(distanceBin >=0 && distanceBin<numBins)
                tmpHist[distanceBin] = tmpHist[distanceBin] + 2;
        }
    }
}

/* FOR ACF */

double acos_fast(double x)
{
    return sqrt(2*(1.-x) + (1.-x)*(1.-x)/3. + (1.-x)*(1.-x)*(1.-x)*4./45.);
}

double measure_distance(double cosTheta1, double phi1, double cosTheta2, double phi2)
{
//     return acos(sin(acos(cosTheta1)) * sin(acos(cosTheta2)) * cos(phi1-phi2) + cosTheta1 * cosTheta2);
    return acos_fast(sqrt(1.-cosTheta1*cosTheta1) * sqrt(1.-cosTheta2*cosTheta2) * cos(phi1-phi2) + cosTheta1 * cosTheta2);
}

int get_distance_theta_bin(double distanceRad, int distanceInLog, int numBins)
{
    // good values for numBins: 10 for log and 30 for linear scaling

    int distanceBin = 0;
    if(distanceInLog != 1)
    {
        distanceBin = distanceRad / PI * 180. * numBins;
    }
    else
    {
        distanceBin = (int)((log10(distanceRad / PI * 180.) + 3.) * (double)numBins);
    }
    
    return distanceBin;
}

void calc_distance_distribution_DD_cosTheta(int num, double *cosTheta, double *phi, int numBins, long int **hist, int numTheta, int distanceInLog)
{
    long int *tmpHist = *hist;
    double distanceRad = 0.;
    int distanceBin = 0;
    
    for(int source1=0; source1<num; source1++)
    {
        for(int source2=source1+1; source2<num; source2++)
        {
            distanceRad = measure_distance(cosTheta[source1], phi[source1], cosTheta[source2], phi[source2]);
            distanceBin = get_distance_theta_bin(distanceRad, distanceInLog, numTheta);
            
            if(distanceBin >=0 && distanceBin<numBins)
                tmpHist[distanceBin] = tmpHist[distanceBin] + 2;
        }
    }
}

void calc_distance_distribution_DR_cosTheta(int num1, double *cosTheta1, double *phi1, int num2, double *cosTheta2, double *phi2, int numBins, long int **hist, int numTheta, int distanceInLog)
{
    long int *tmpHist = *hist;
    double distanceRad = 0.;
    int distanceBin = 0;
    
    for(int source1=0; source1<num1; source1++)
    {
        for(int source2=0; source2<num2; source2++)
        {
            distanceRad = measure_distance(cosTheta1[source1], phi1[source1], cosTheta2[source2], phi2[source2]);
            distanceBin = get_distance_theta_bin(distanceRad, distanceInLog, numTheta);
            
            if(distanceBin >=0 && distanceBin<numBins)
                tmpHist[distanceBin]++;
        }
    }
}

void calc_distance_distribution_DD_cosTheta_over_ranks(int num1, double *cosTheta1, double *phi1, int num2, double *cosTheta2, double *phi2, int numBins, long int **hist, int numTheta, int distanceInLog)
{
    long int *tmpHist = *hist;
    double distanceRad = 0.;
    int distanceBin = 0;
    
    for(int source1=0; source1<num1; source1++)
    {
        for(int source2=0; source2<num2; source2++)
        {
            distanceRad = measure_distance(cosTheta1[source1], phi1[source1], cosTheta2[source2], phi2[source2]);
            distanceBin = get_distance_theta_bin(distanceRad, distanceInLog, numTheta);
            
            if(distanceBin >=0 && distanceBin<numBins)
                tmpHist[distanceBin] = tmpHist[distanceBin] + 2;
        }
    }
}

/* --------------------------------------------------------- */
/* COMMUNICATING (& COMPUTING) DISTANCE HISTOGRAM              */
/* --------------------------------------------------------- */

void count_pairs_r(domain_t *thisDomain, int numSources, double *xsource, double *ysource, double *zsource, int numRand, double *xrand, double *yrand, double *zrand, int numBins, long int **histDD, long int **histDR, long int **histRR, int distanceInLog)
{
    /* INITIALIZATION */
    int thisRank = thisDomain->thisRank;
    int size = thisDomain->size;
    int destRank = thisRank;
    int sendingRank = thisRank;
    
    int recvNumSources = 0;
    double *recvXsource = NULL, *recvYsource = NULL, *recvZsource = NULL;
    
    int recvNumRand = 0;
    double *recvXrand = NULL, *recvYrand = NULL, *recvZrand = NULL;
    
    /* LOOP OVER ALL COMBINATIONS OF RANKS */
    if(thisRank == 0) printf("Computing distances between pairs\n");
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    for(int j=0; j<0.5*size+1; j++)
    {
        /* rank to send galaxies to: destRank is by j different to thisRank */
        destRank = (thisRank + size + j) % size;
        /* rank from which galaxies are received */
        sendingRank = (thisRank + size - j) % size;
                
#ifdef __MPI
        if(destRank != thisRank && sendingRank != thisRank)
        {
            if(size%2 == 0 && j >= size/2)
            {
                send_recv_array_separated(thisRank, destRank, numSources, xsource, sendingRank, &recvNumSources, &recvXsource);
                send_recv_array_separated(thisRank, destRank, numSources, ysource, sendingRank, &recvNumSources, &recvYsource);
                send_recv_array_separated(thisRank, destRank, numSources, zsource, sendingRank, &recvNumSources, &recvZsource);

                send_recv_array_separated(thisRank, destRank, numRand, xrand, sendingRank, &recvNumRand, &recvXrand);
                send_recv_array_separated(thisRank, destRank, numRand, yrand, sendingRank, &recvNumRand, &recvYrand);
                send_recv_array_separated(thisRank, destRank, numRand, zrand, sendingRank, &recvNumRand, &recvZrand);
            }
            else{                
                send_recv_array(thisRank, destRank, numSources, xsource, sendingRank, &recvNumSources, &recvXsource);
                send_recv_array(thisRank, destRank, numSources, ysource, sendingRank, &recvNumSources, &recvYsource);
                send_recv_array(thisRank, destRank, numSources, zsource, sendingRank, &recvNumSources, &recvZsource);

                send_recv_array(thisRank, destRank, numRand, xrand, sendingRank, &recvNumRand, &recvXrand);
                send_recv_array(thisRank, destRank, numRand, yrand, sendingRank, &recvNumRand, &recvYrand);
                send_recv_array(thisRank, destRank, numRand, zrand, sendingRank, &recvNumRand, &recvZrand);
            }
    
            /* compute distances (use DR routine!) between DD', DR', RD', RR' (check for double counts!) */
            if(recvXsource != NULL)
            {
                calc_distance_distribution_DD_r_over_ranks(numSources, xsource, ysource, zsource, recvNumSources, recvXsource, recvYsource, recvZsource, numBins, histDD, distanceInLog);

                calc_distance_distribution_DR_r(numSources, xsource, ysource, zsource, recvNumRand, recvXrand, recvYrand, recvZrand, numBins, histDR, distanceInLog);
                calc_distance_distribution_DR_r(numRand, xrand, yrand, zrand, recvNumSources, recvXsource, recvYsource, recvZsource, numBins, histDR, distanceInLog);
                
                calc_distance_distribution_DD_r_over_ranks(numRand, xrand, yrand, zrand, recvNumRand, recvXrand, recvYrand, recvZrand, numBins, histRR, distanceInLog);
            }

            if(recvXsource != NULL)
            {
                free(recvXsource);
                recvXsource = NULL;
            }
            if(recvYsource != NULL)
            {
                free(recvYsource);
                recvYsource = NULL;
            }
            if(recvZsource != NULL)
            {
                free(recvZsource);
                recvZsource = NULL;
            }
            if(recvXrand != NULL)
            {
                free(recvXrand);
                recvXrand = NULL;
            }
            if(recvYrand != NULL)
            {
                free(recvYrand);
                recvYrand = NULL;
            }
            if(recvZrand != NULL)
            {
                free(recvZrand);
                recvZrand = NULL;
            }
        }
#endif    

        if(j == 0)
        {
            printf("  rank %d: calculating DD\n", thisDomain->thisRank);
            calc_distance_distribution_DD_r(numSources, xsource, ysource, zsource, numBins, histDD, distanceInLog);
            printf("  rank %d: calculating DR\n", thisDomain->thisRank);
            calc_distance_distribution_DR_r(numSources, xsource, ysource, zsource, numRand, xrand, yrand, zrand, numBins, histDR, distanceInLog);
            printf("  rank %d: calculating RR\n", thisDomain->thisRank);
            calc_distance_distribution_DD_r(numRand, xrand, yrand, zrand, numBins, histRR, distanceInLog);
            printf("  rank %d: calculating DD, DR and RR across ranks\n", thisRank);
        }
#ifdef __MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }
}

void count_pairs_cosTheta(domain_t *thisDomain, int numSources, double *cosTheta, double *phi, int numRand, double *cosThetaRand, double *phiRand, int numBins, long int **histDD, long int **histDR, long int **histRR, int numTheta, int distanceInLog)
{
    /* INITIALIZATION */
    int thisRank = thisDomain->thisRank;
    int size = thisDomain->size;
    int destRank = thisRank;
    int sendingRank = thisRank;
    
    int recvNumSources = 0;
    double *recvCosTheta = NULL, *recvPhi = NULL;
    
    int recvNumRand = 0;
    double *recvCosThetaRand = NULL, *recvPhiRand = NULL;
    
    /* LOOP OVER ALL COMBINATIONS OF RANKS */
    if(thisRank == 0) printf("Computing distances between pairs\n");
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    for(int j=0; j<0.5*size+1; j++)
    {
        /* rank to send galaxies to: destRank is by j different to thisRank */
        destRank = (thisRank + size + j) % size;
        /* rank from which galaxies are received */
        sendingRank = (thisRank + size - j) % size;
                
#ifdef __MPI
        if(destRank != thisRank && sendingRank != thisRank)
        {
            if(size%2 == 0 && j >= size/2)
            {
                send_recv_array_separated(thisRank, destRank, numSources, cosTheta, sendingRank, &recvNumSources, &recvCosTheta);
                send_recv_array_separated(thisRank, destRank, numSources, phi, sendingRank, &recvNumSources, &recvPhi);

                send_recv_array_separated(thisRank, destRank, numRand, cosThetaRand, sendingRank, &recvNumRand, &recvCosThetaRand);
                send_recv_array_separated(thisRank, destRank, numRand, phiRand, sendingRank, &recvNumRand, &recvPhiRand);
            }
            else{                
                send_recv_array(thisRank, destRank, numSources, cosTheta, sendingRank, &recvNumSources, &recvCosTheta);
                send_recv_array(thisRank, destRank, numSources, phi, sendingRank, &recvNumSources, &recvPhi);

                send_recv_array(thisRank, destRank, numRand, cosThetaRand, sendingRank, &recvNumRand, &recvCosThetaRand);
                send_recv_array(thisRank, destRank, numRand, phiRand, sendingRank, &recvNumRand, &recvPhiRand);
            }
    
            /* compute distances (use DR routine!) between DD', DR', RD', RR' (check for double counts!) */
            if(recvCosTheta != NULL)
            {
                calc_distance_distribution_DD_cosTheta_over_ranks(numSources, cosTheta, phi, recvNumSources, recvCosTheta, recvPhi, numBins, histDD, numTheta, distanceInLog);

                calc_distance_distribution_DR_cosTheta(numSources, cosTheta, phi, recvNumRand, recvCosThetaRand, recvPhiRand, numBins, histDR, numTheta, distanceInLog);
                calc_distance_distribution_DR_cosTheta(numRand, cosThetaRand, phiRand, recvNumSources, recvCosTheta, recvPhi, numBins, histDR, numTheta, distanceInLog);
                
                calc_distance_distribution_DD_cosTheta_over_ranks(numRand, cosThetaRand, phiRand, recvNumRand, recvCosThetaRand, recvPhiRand, numBins, histRR, numTheta, distanceInLog);
            }

            if(recvCosTheta != NULL)
            {
                free(recvCosTheta);
                recvCosTheta = NULL;
            }
            if(recvPhi != NULL)
            {
                free(recvPhi);
                recvPhi = NULL;
            }
            if(recvCosThetaRand != NULL)
            {
                free(recvCosThetaRand);
                recvCosThetaRand = NULL;
            }
            if(recvPhiRand != NULL)
            {
                free(recvPhiRand);
                recvPhiRand = NULL;
            }
        }
#endif    

        if(j == 0)
        {
            printf("  rank %d: calculating DD\n", thisDomain->thisRank);
            calc_distance_distribution_DD_cosTheta(numSources, cosTheta, phi, numBins, histDD, numTheta, distanceInLog);
            printf("  rank %d: calculating DR\n", thisDomain->thisRank);
            calc_distance_distribution_DR_cosTheta(numSources, cosTheta, phi, numRand, cosThetaRand, phiRand, numBins, histDR, numTheta, distanceInLog);
            printf("  rank %d: calculating RR\n", thisDomain->thisRank);
            calc_distance_distribution_DD_cosTheta(numRand, cosThetaRand, phiRand, numBins, histRR, numTheta, distanceInLog);
            printf("  rank %d: calculating DD, DR and RR across ranks\n", thisRank);
        }

    }
}

#ifdef __MPI
void collapse_histogram_across_ranks(int numBins, long int **hist)
{
    long int *hist_global = allocate_array_long_int(numBins, "hist");
    
    MPI_Allreduce(*hist, hist_global, numBins, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    
    for(int i=0; i<numBins; i++)
        (*hist)[i] = hist_global[i];
    
    free(hist_global);
}
#endif 

#ifdef __MPI
void collapse_pairs_across_ranks(int numBins, long int **histDD, long int **histDR, long int **histRR)
{
    collapse_histogram_across_ranks(numBins, histDD);
    collapse_histogram_across_ranks(numBins, histDR);
    collapse_histogram_across_ranks(numBins, histRR);
}
#endif

#ifdef __MPI
void collapse_numPairs_across_ranks(int *numSources, int *numRand)
{
    int recvNumSources = 0;
    int recvNumRand = 0;
    
    MPI_Allreduce(numSources, &recvNumSources, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(numRand, &recvNumRand, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    *numSources = recvNumSources;
    *numRand = recvNumRand;
}
#endif

/* --------------------------------------------------------- */
/* CALCULATING CORRELATION FUNCTION                          */
/* --------------------------------------------------------- */

/* r bins need to range from 0 to 1 */
double *calc_r_bins(int numBins, int distanceInLog)
{
    double *distance = allocate_array_double(numBins, "distance");
    
    if(distanceInLog != 1)
    {
        for(int i=0; i<numBins; i++)
            distance[i] = 2. * (double)i / (double) numBins;
    }
    else
    {
        for(int i=0; i<numBins; i++)
            distance[i] = 2.*pow(1.e3, (double)i / (double)numBins - 1.);
    }
    
    return distance;
}

/* theta bins need to range from 0 to 180 degrees */
double *calc_theta_bins(int numBins, int distanceInLog, int numTheta)   // gives result in degrees
{
    double *distanceDeg = allocate_array_double(numBins, "distanceRad");
    
    if(distanceInLog != 1)
    {
        for(int i=0; i<numBins; i++)
            distanceDeg[i] = (double)i / (double)numBins;
    }
    else
    {
        for(int i=0; i<numBins; i++)
            distanceDeg[i] = pow(10., (double)i / (double)numTheta - 3.);
    }
    
    return distanceDeg;
}

void generate_histograms_r(int numBins, long int **histDD, long int **histDR, long int **histRR, double **distanceArray, int distanceInLog)
{
    *histDD = allocate_array_long_int(numBins, "histDD");
    *histRR = allocate_array_long_int(numBins, "histRR");
    *histDR = allocate_array_long_int(numBins, "histDR");
    
    *distanceArray = calc_r_bins(numBins, distanceInLog);
}

void generate_histograms_theta(int numBins, long int **histDD, long int **histDR, long int **histRR, double **distanceArray, int distanceInLog, int numTheta)
{
    *histDD = allocate_array_long_int(numBins, "histDD");
    *histRR = allocate_array_long_int(numBins, "histRR");
    *histDR = allocate_array_long_int(numBins, "histDR");
    
    *distanceArray = calc_theta_bins(numBins, distanceInLog, numTheta);
}

void deallocate_histograms(double *distanceArray, long int *histDD, long int *histDR, long int *histRR)
{
    if(distanceArray != NULL) free(distanceArray);
    if(histDD != NULL) free(histDD);
    if(histDR != NULL) free(histDR);
    if(histRR != NULL) free(histRR);
}

double *calc_correlation_landy_szalay(int numBins, long int **histDD, long int **histDR, long int **histRR, int numData, int numRand)
{
    double *correlation = allocate_array_double(numBins, "correlation");
    
    double fact1 = (double)(numRand - 1) / (double)(numData - 1);
    double fact2 = (double)numRand / (double)(numData - 1);
    
    long int *tmpHistDD = *histDD;
    long int *tmpHistDR = *histDR;
    long int *tmpHistRR = *histRR;
    
    for(int i=0; i<numBins; i++)
    {
        correlation[i] = fact1 * fact1 * (double)tmpHistDD[i] / (double)tmpHistRR[i] - 2. * fact2 * (double)tmpHistDR[i] / (double)tmpHistRR[i] + 1.; 
    }
    
    return correlation;
}

void write_corrfunc_to_file(domain_t *thisDomain, int numBins, double *distanceArray, double *correlation, char *type, char *filename)
{
    FILE *f;
    
    if(thisDomain->thisRank == 0)
    {
        f = fopen(filename, "w");
        
        if(strcmp(type, "3D") == 0)
        {
            fprintf(f, "# distance [h^-1 Mpc]\t3D correlation function");
            printf("distance [h^-1 Mpc]\t3D correlation function\n");
        }
        else if(strcmp(type, "ACF") == 0)
        {
            fprintf(f, "# theta [deg]\tangular correlation function");
            printf("theta [deg]\tangular correlation function\n");
        }
        
        for(int i=0; i<numBins; i++)
        {
            printf("%e\t%e\n", distanceArray[i], correlation[i]);
            fprintf(f, "%e\t%e\n", distanceArray[i], correlation[i]);
        }
        
        fclose(f);
    }
}

/* --------------------------------------------------------- */

void calc_3Dcorrfunc_cartesian(int numSources, double *xsource, double *ysource, double *zsource, boxparams_t *thisBoxparams, domain_t *thisDomain, char *filename)
{
    int numRand = 0;
    double *xrand = NULL, *yrand = NULL, *zrand = NULL;
    
    int numBins = 30;
    int distanceInLog = 0;
    long int *histDD = NULL, *histDR = NULL, *histRR = NULL;
    double *distance = NULL;
    double *correlation = NULL;
    
    /* generate histograms for pairs */
    generate_histograms_r(numBins, &histDD, &histDR, &histRR, &distance, distanceInLog);
    
    /* generate randoms */
    generate_randoms_cartesian(thisDomain, thisBoxparams, numSources, &numRand, &xrand, &yrand, &zrand);
    
    /* count pairs across ranks */
    count_pairs_r(thisDomain, numSources, xsource, ysource, zsource, numRand, xrand, yrand, zrand, numBins, &histDD, &histDR, &histRR, distanceInLog);

#ifdef __MPI
    collapse_pairs_across_ranks(numBins, &histDD, &histDR, &histRR);
    collapse_numPairs_across_ranks(&numSources, &numRand);
#endif
    correlation = calc_correlation_landy_szalay(numBins, &histDD, &histDR, &histRR, numSources, numRand);
    write_corrfunc_to_file(thisDomain, numBins, distance, correlation, "3D", filename);
    
    free(correlation);
    deallocate_randoms_cartesian(xrand, yrand, zrand);
    deallocate_histograms(distance, histDD, histDR, histRR);
}


void calc_ACF_cartesian(int numSources, double *xsource, double *ysource, double *zsource, boxparams_t *thisBoxparams, cosmparams_t *thisCosmparams, domain_t *thisDomain, char *filename)
{
    double *r = NULL, *cosTheta = NULL, *phi = NULL;
    int numRand = 0;
    double *rRand = NULL, *cosThetaRand = NULL, *phiRand = NULL;
    
    int numBins = 30;
    int numTheta = 10;
    int distanceInLog = 0;
    long int *histDD = NULL, *histDR = NULL, *histRR = NULL;
    double *distanceDeg = NULL;
    double *correlation = NULL;
    
    /* generate histograms for pairs */
    generate_histograms_theta(numBins, &histDD, &histDR, &histRR, &distanceDeg, distanceInLog, numTheta);
    
    /* convert x,y,z into ra,dec*/
    generate_ra_dec(thisCosmparams, thisBoxparams, numSources, xsource, ysource, zsource, &r, &cosTheta, &phi);

    /* generate randoms */
    generate_randoms_ra_dec(thisDomain, thisCosmparams, thisBoxparams, numSources, &numRand, &rRand, &cosThetaRand, &phiRand);

    /* count pairs across ranks */
    count_pairs_cosTheta(thisDomain, numSources, cosTheta, phi, numRand, cosThetaRand, phiRand, numBins, &histDD, &histDR, &histRR, numTheta, distanceInLog);

#ifdef __MPI
    collapse_pairs_across_ranks(numBins, &histDD, &histDR, &histRR);
    collapse_numPairs_across_ranks(&numSources, &numRand);
#endif
    correlation = calc_correlation_landy_szalay(numBins, &histDD, &histDR, &histRR, numSources, numRand);
    write_corrfunc_to_file(thisDomain, numBins, distanceDeg, correlation, "ACF", filename);
    
    free(correlation);
    deallocate_ra_dec(r, cosTheta, phi);
    deallocate_randoms_ra_dec(rRand, cosThetaRand, phiRand);
    deallocate_histograms(distanceDeg, histDD, histDR, histRR);
}


void calc_ACF_ra_dec(int numSources, double *cosTheta, double *phi, boxparams_t *thisBoxparams, cosmparams_t *thisCosmparams, domain_t *thisDomain, char *filename)
{
    int numRand = 0;
    double *rRand = NULL, *cosThetaRand = NULL, *phiRand = NULL;
        
    int numBins = 30;
    int numTheta = 10;
    int distanceInLog = 0;
    long int *histDD = NULL, *histDR = NULL, *histRR = NULL;
    double *distanceDeg = NULL;
    double *correlation = NULL;
    
    /* generate histograms for pairs */
    generate_histograms_theta(numBins, &histDD, &histDR, &histRR, &distanceDeg, distanceInLog, numTheta);

    /* generate randoms */
    generate_randoms_ra_dec(thisDomain, thisCosmparams, thisBoxparams, numSources, &numRand, &rRand, &cosThetaRand, &phiRand);
    
    /* count pairs across ranks */
    count_pairs_cosTheta(thisDomain, numSources, cosTheta, phi, numRand, cosThetaRand, phiRand, numBins, &histDD, &histDR, &histRR, numTheta, distanceInLog);
    
#ifdef __MPI
    collapse_pairs_across_ranks(numBins, &histDD, &histDR, &histRR);
    collapse_numPairs_across_ranks(&numSources, &numRand);
#endif
    correlation = calc_correlation_landy_szalay(numBins, &histDD, &histDR, &histRR, numSources, numRand);
    write_corrfunc_to_file(thisDomain, numBins, distanceDeg, correlation, "ACF", filename);

    free(correlation);
    deallocate_randoms_ra_dec(rRand, cosThetaRand, phiRand);
    deallocate_histograms(distanceDeg, histDD, histDR, histRR);
}
