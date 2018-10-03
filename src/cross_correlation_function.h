#ifndef CROSS_CORRELATION_FUNCTION_H
#define CROSS_CORRELATION_FUNCTION_H

/* --------------------------------------------------------- */
/* COMMUNICATING (& COMPUTING) DISTANCE HISTOGRAM            */
/* --------------------------------------------------------- */

void count_pairs_DR_r(domain_t *thisDomain, int numSources, double *xsource, double *ysource, double *zsource, int numRand, double *xrand, double *yrand, double *zrand, int numBins, long int **histDR, int distanceInLog);
void count_pairs_DR_cosTheta(domain_t *thisDomain, int numSources, double *cosTheta, double *phi, int numRand, double *cosThetaRand, double *phiRand, int numBins, long int **histDR, int numTheta, int distanceInLog);

#ifdef __MPI
void ccf_collapse_pairs_across_ranks(int numBins, long int **histD1D2, long int **histD1R2, long int **histR1D2, long int **histR1R2);
#endif

/* --------------------------------------------------------- */
/* CALCULATING CORRELATION FUNCTION                          */
/* --------------------------------------------------------- */

void ccf_generate_histograms_r(int numBins, long int **histD1D2, long int **histD1R2, long int **histR1D2, long int **histR1R2, double **distanceArray, int distanceInLog);
void ccf_generate_histograms_theta(int numBins, long int **histD1D2, long int **histD1R2, long int **histR1D2, long int **histR1R2, double **distanceArray, int distanceInLog, int numTheta);
void ccf_deallocate_histograms(double *distanceArray, long int *histD1D2, long int *histD1R2, long int *histR1D2, long int *histR1R2);
double *calc_cross_correlation_landy_szalay(int numBins, long int **histD1D2, long int **histD1R2, long int **histR1D2, long int **histR1R2, int numData1, int numData2, int numRand1, int numRand2);

/* --------------------------------------------------------- */

void calc_cross_3Dcorrfunc_cartesian(int numSources1, double *x1, double *y1, double *z1, int numSources2, double *x2, double *y2, double *z2, boxparams_t *thisBoxparams, domain_t *thisDomain, char *filename);
void calc_cross_ACF_cartesian(int numSources1, double *x1, double *y1, double *z1, int numSources2, double *x2, double *y2, double *z2, boxparams_t *thisBoxparams, cosmparams_t *thisCosmparams, domain_t *thisDomain, char *filename);


#endif
