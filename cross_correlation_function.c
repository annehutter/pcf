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
/* COMMUNICATING (& COMPUTING) DISTANCE HISTOGRAM            */
/* --------------------------------------------------------- */

void count_pairs_DR_r(domain_t *thisDomain, int numSources, double *xsource, double *ysource, double *zsource, int numRand, double *xrand, double *yrand, double *zrand, int numBins, long int **histDR, int distanceInLog)
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
                calc_distance_distribution_DR_r(numSources, xsource, ysource, zsource, recvNumRand, recvXrand, recvYrand, recvZrand, numBins, histDR, distanceInLog);
                calc_distance_distribution_DR_r(numRand, xrand, yrand, zrand, recvNumSources, recvXsource, recvYsource, recvZsource, numBins, histDR, distanceInLog);
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
            printf("  rank %d: calculating DR\n", thisDomain->thisRank);
            calc_distance_distribution_DR_r(numSources, xsource, ysource, zsource, numRand, xrand, yrand, zrand, numBins, histDR, distanceInLog);
        }
#ifdef __MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }
}

void count_pairs_DR_cosTheta(domain_t *thisDomain, int numSources, double *cosTheta, double *phi, int numRand, double *cosThetaRand, double *phiRand, int numBins, long int **histDR, int numTheta, int distanceInLog)
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
                calc_distance_distribution_DR_cosTheta(numSources, cosTheta, phi, recvNumRand, recvCosThetaRand, recvPhiRand, numBins, histDR, numTheta, distanceInLog);
                calc_distance_distribution_DR_cosTheta(numRand, cosThetaRand, phiRand, recvNumSources, recvCosTheta, recvPhi, numBins, histDR, numTheta, distanceInLog);
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
            printf("  rank %d: calculating DR\n", thisDomain->thisRank);
            calc_distance_distribution_DR_cosTheta(numSources, cosTheta, phi, numRand, cosThetaRand, phiRand, numBins, histDR, numTheta, distanceInLog);
        }
    }
}

#ifdef __MPI
void ccf_collapse_pairs_across_ranks(int numBins, long int **histD1D2, long int **histD1R2, long int **histR1D2, long int **histR1R2)
{
    collapse_histogram_across_ranks(numBins, histD1D2);
    collapse_histogram_across_ranks(numBins, histD1R2);
    collapse_histogram_across_ranks(numBins, histR1D2);
    collapse_histogram_across_ranks(numBins, histR1R2);
}
#endif

/* --------------------------------------------------------- */
/* CALCULATING CORRELATION FUNCTION                          */
/* --------------------------------------------------------- */

void ccf_generate_histograms_r(int numBins, long int **histD1D2, long int **histD1R2, long int **histR1D2, long int **histR1R2, double **distanceArray, int distanceInLog)
{
    *histD1D2 = allocate_array_long_int(numBins, "histD1D2");
    *histD1R2 = allocate_array_long_int(numBins, "histD1R2");
    *histR1D2 = allocate_array_long_int(numBins, "histR1D2");
    *histR1R2 = allocate_array_long_int(numBins, "histR1R2");

    *distanceArray = calc_r_bins(numBins, distanceInLog);
}

void ccf_generate_histograms_theta(int numBins, long int **histD1D2, long int **histD1R2, long int **histR1D2, long int **histR1R2, double **distanceArray, int distanceInLog, int numTheta)
{
    *histD1D2 = allocate_array_long_int(numBins, "histD1D2");
    *histD1R2 = allocate_array_long_int(numBins, "histD1R2");
    *histR1D2 = allocate_array_long_int(numBins, "histR1D2");
    *histR1R2 = allocate_array_long_int(numBins, "histR1R2");
    
    *distanceArray = calc_theta_bins(numBins, distanceInLog, numTheta);
}

void ccf_deallocate_histograms(double *distanceArray, long int *histD1D2, long int *histD1R2, long int *histR1D2, long int *histR1R2)
{
    if(distanceArray != NULL) free(distanceArray);
    if(histD1D2 != NULL) free(histD1D2);
    if(histD1R2 != NULL) free(histD1R2);
    if(histR1D2 != NULL) free(histR1D2);
    if(histR1R2 != NULL) free(histR1R2);
}

double *calc_cross_correlation_landy_szalay(int numBins, long int **histD1D2, long int **histD1R2, long int **histR1D2, long int **histR1R2, int numData1, int numData2, int numRand1, int numRand2)
{
    double *correlation = allocate_array_double(numBins, "correlation");
    
    double fact1 = (double)numRand1 / (double)numData1;
    double fact2 = (double)numRand2 / (double)numData2;
    
    long int *tmpHistD1D2 = *histD1D2;
    long int *tmpHistD1R2 = *histD1R2;
    long int *tmpHistR1D2 = *histR1D2;
    long int *tmpHistR1R2 = *histR1R2;
    
    for(int i=0; i<numBins; i++)
    {
        correlation[i] = fact1 * fact2 * (double)tmpHistD1D2[i] / (double)tmpHistR1R2[i] - fact1 * (double)tmpHistD1R2[i] / (double)tmpHistR1R2[i] - fact2 * (double)tmpHistR1D2[i] / (double)tmpHistR1R2[i] + 1.; 
    }
    
    return correlation;
}

/* --------------------------------------------------------- */

void calc_cross_3Dcorrfunc_cartesian(int numSources1, double *x1, double *y1, double *z1, int numSources2, double *x2, double *y2, double *z2, boxparams_t *thisBoxparams, domain_t *thisDomain, char *filename)
{
    int numRand1 = 0, numRand2 = 0;
    double *xr1 = NULL, *yr1 = NULL, *zr1 = NULL;
    double *xr2 = NULL, *yr2 = NULL, *zr2 = NULL;
    
    int numBins = 30;
    int distanceInLog = 0;
    long int *histD1D2 = NULL, *histD1R2 = NULL, *histR1D2 = NULL, *histR1R2 = NULL;
    double *distance = NULL;
    double *correlation = NULL;
    
    /* generate histograms for pairs */
    ccf_generate_histograms_r(numBins, &histD1D2, &histD1R2, &histR1D2, &histR1R2, &distance, distanceInLog);
    
    /* generate randoms */
    generate_randoms_cartesian(thisDomain, thisBoxparams, numSources1, &numRand1, &xr1, &yr1, &zr1);
    generate_randoms_cartesian(thisDomain, thisBoxparams, numSources2, &numRand2, &xr2, &yr2, &zr2);

    /* count pairs across ranks */
    count_pairs_DR_r(thisDomain, numSources1, x1, y1, z1, numSources2, x2, y2, z2, numBins, &histD1D2, distanceInLog);
    count_pairs_DR_r(thisDomain, numSources1, x1, y1, z1, numRand2, xr2, yr2, zr2, numBins, &histD1R2, distanceInLog);
    count_pairs_DR_r(thisDomain, numRand1, xr1, yr1, zr1, numSources2, x2, y2, z2, numBins, &histR1D2, distanceInLog);
    count_pairs_DR_r(thisDomain, numRand1, xr1, yr1, zr1, numRand2, xr2, yr2, zr2, numBins, &histR1R2, distanceInLog);
    
#ifdef __MPI
    ccf_collapse_pairs_across_ranks(numBins, &histD1D2, &histD1R2, &histR1D2, &histR1R2);
    collapse_numPairs_across_ranks(&numSources1, &numRand1);
    collapse_numPairs_across_ranks(&numSources2, &numRand2);
#endif
    correlation = calc_cross_correlation_landy_szalay(numBins, &histD1D2, &histD1R2, &histR1D2, &histR1R2, numSources1, numSources2, numRand1, numRand2);
    write_corrfunc_to_file(thisDomain, numBins, distance, correlation, "3D", filename);
    
    free(correlation);
    deallocate_randoms_cartesian(xr1, yr1, zr1);
    deallocate_randoms_cartesian(xr2, yr2, zr2);
    ccf_deallocate_histograms(distance, histD1D2, histD1R2, histR1D2, histR1R2);
}


void calc_cross_ACF_cartesian(int numSources1, double *x1, double *y1, double *z1, int numSources2, double *x2, double *y2, double *z2, boxparams_t *thisBoxparams, cosmparams_t *thisCosmparams, domain_t *thisDomain, char *filename)
{
    double *r1 = NULL, *cosTheta1 = NULL, *phi1 = NULL;
    double *r2 = NULL, *cosTheta2 = NULL, *phi2 = NULL;
    int numRand1 = NULL, numRand2 = NULL;
    double *rRand1 = NULL, *cosThetaRand1 = NULL, *phiRand1 = NULL;
    double *rRand2 = NULL, *cosThetaRand2 = NULL, *phiRand2 = NULL;

    int numBins = 30;
    int numTheta = 10;
    int distanceInLog = 0;
    long int *histD1D2 = NULL, *histD1R2 = NULL, *histR1D2 = NULL, *histR1R2 = NULL;
    double *distanceDeg = NULL;
    double *correlation = NULL;
    
    /* generate histograms for pairs */
    ccf_generate_histograms_theta(numBins, &histD1D2, &histD1R2, &histR1D2, &histR1R2, &distanceDeg, distanceInLog, numTheta);
    
    /* convert x,y,z into ra,dec*/
    generate_ra_dec(thisCosmparams, thisBoxparams, numSources1, x1, y1, z1, &r1, &cosTheta1, &phi1);
    generate_ra_dec(thisCosmparams, thisBoxparams, numSources2, x2, y2, z2, &r2, &cosTheta2, &phi2);

    /* generate randoms */
    generate_randoms_ra_dec(thisDomain, thisCosmparams, thisBoxparams, numSources1, &numRand1, &rRand1, &cosThetaRand1, &phiRand1);
    generate_randoms_ra_dec(thisDomain, thisCosmparams, thisBoxparams, numSources2, &numRand2, &rRand2, &cosThetaRand2, &phiRand2);

    /* count pairs across ranks */
    count_pairs_DR_cosTheta(thisDomain, numSources1, cosTheta1, phi1, numSources2, cosTheta2, phi2, numBins, &histD1D2, numTheta, distanceInLog);
    count_pairs_DR_cosTheta(thisDomain, numSources1, cosTheta1, phi1, numRand2, cosThetaRand2, phiRand2, numBins, &histD1R2, numTheta, distanceInLog);
    count_pairs_DR_cosTheta(thisDomain, numRand1, cosThetaRand1, phiRand1, numSources2, cosTheta2, phi2, numBins, &histR1D2, numTheta, distanceInLog);
    count_pairs_DR_cosTheta(thisDomain, numRand1, cosThetaRand1, phiRand1, numRand2, cosThetaRand2, phiRand2, numBins, &histR1R2, numTheta, distanceInLog);

#ifdef __MPI
    ccf_collapse_pairs_across_ranks(numBins, &histD1D2, &histD1R2, &histR1D2, &histR1R2);
    collapse_numPairs_across_ranks(&numSources1, &numRand1);
    collapse_numPairs_across_ranks(&numSources2, &numRand2);
#endif
    correlation = calc_cross_correlation_landy_szalay(numBins, &histD1D2, &histD1R2, &histR1D2, &histR1R2, numSources1, numSources2, numRand1, numRand2);
    write_corrfunc_to_file(thisDomain, numBins, distanceDeg, correlation, "ACF", filename);
    
    free(correlation);
    deallocate_ra_dec(r1, cosTheta1, phi1);
    deallocate_ra_dec(r2, cosTheta2, phi2);
    deallocate_randoms_ra_dec(rRand1, cosThetaRand1, phiRand1);
    deallocate_randoms_ra_dec(rRand2, cosThetaRand2, phiRand2);
    ccf_deallocate_histograms(distanceDeg, histD1D2, histD1R2, histR1D2, histR1R2);

}


