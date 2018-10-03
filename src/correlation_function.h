#ifndef CORRELATION_FUNCTION_H
#define CORRELATION_FUNCTION_H

typedef struct{
    double h;
    double omega_m;
    double omega_l;
    double redshift;
}cosmparams_t;

typedef struct{
    double boxsize;
    
    double xboxsize;
    double yboxsize;
    double zboxsize;
    
    double xlow;
    double ylow;
    double zlow;
    
    double xup;
    double yup;
    double zup;
    
    double xc;
    double yc;
    double zc;
    
    int dir;
}boxparams_t;


cosmparams_t *initCosmparams();
cosmparams_t *read_params_to_cosmparams(double h, double omega_m , double omega_l, double redshift);
void deallocate_cosmparams(cosmparams_t *thisCosmparams);
boxparams_t *initBoxparams();
boxparams_t *read_params_to_boxparams(double boxsize, double xboxsize, double yboxsize, double zboxsize, double xlow, double ylow, double zlow, int dir);
void deallocate_boxparams(boxparams_t *thisBoxparams);

/* --------------------------------------------------------- */
/* SELECT SUBBOX IN SIMULATION BOX                           */
/* --------------------------------------------------------- */

void select_subbox_cartesian(boxparams_t *thisBoxparams, int *numSources, double **xsource, double **ysource, double **zsource);

/* --------------------------------------------------------- */
/* GENERATING RANDOMS                                        */
/* --------------------------------------------------------- */

/* RANDOMS FOR 3D GRID: x, y, z */

void generate_randoms_cartesian(domain_t *thisDomain, boxparams_t *thisBoxparams, int numSources, int *numRand, double **xrand, double **yrand, double **zrand);
void deallocate_randoms_cartesian(double *xrand, double *yrand, double *zrand);

/* RANDOMS FOR ACF: cos(theta), phi */

void generate_randoms_ra_dec(domain_t *thisDomain, cosmparams_t *thisCosmparams, boxparams_t *thisBoxparams, int numSources, int *numRand, double **rRand, double **cosThetaRand, double **phiRand);
void deallocate_randoms_ra_dec(double *rRand, double *cosThetaRand, double *phiRand);

/* --------------------------------------------------------- */
/* GENERATING COSTHETA, PHI FROM X,Y,Z                       */
/* --------------------------------------------------------- */

void generate_ra_dec(cosmparams_t *thisCosmparams, boxparams_t *thisBoxparams, int numSources, double *xsource, double *ysource, double *zsource, double **r, double **cosTheta, double **phi);
void deallocate_ra_dec(double *r, double *cosTheta, double *phi);

/* --------------------------------------------------------- */
/* MPI COMMUNICATION ROUTINES                                */
/* --------------------------------------------------------- */

#ifdef __MPI
void send_recv_array(int thisRank, int destRank, int num, double *array, int sendingRank, int *recvNum, double **recvArray);
void send_recv_array_separated(int thisRank, int destRank, int num, double *array, int sendingRank, int *recvNum, double **recvArray);
#endif

/* --------------------------------------------------------- */
/* COMPUTING DISTANCE HISTOGRAM                              */
/* --------------------------------------------------------- */

/* FOR 3D GRIDS */
double measure_distance_3D(double x1, double y1, double z1, double x2, double y2, double z2);
int get_distance_r_bin(double distance, int distanceInLog, int numBins);
void calc_distance_distribution_DD_r(int num, double *xsource, double *ysource, double *zsource, int numBins, long int **hist, int distanceInLog);
void calc_distance_distribution_DR_r(int num1, double *xsource1, double *ysource1, double *zsource1, int num2, double *xsource2, double *ysource2, double *zsource2, int numBins, long int **hist, int distanceInLog);
void calc_distance_distribution_DD_r_over_ranks(int num1, double *xsource1, double *ysource1, double *zsource1, int num2, double *xsource2, double *ysource2, double *zsource2, int numBins, long int **hist, int distanceInLog);

/* FOR ACF */
double acos_fast(double x);
double measure_distance(double cosTheta1, double phi1, double cosTheta2, double phi2);
int get_distance_theta_bin(double distanceRad, int distanceInLog, int numBins);
void calc_distance_distribution_DD_cosTheta(int num, double *cosTheta, double *phi, int numBins, long int **hist, int numTheta, int distanceInLog);
void calc_distance_distribution_DR_cosTheta(int num1, double *cosTheta1, double *phi1, int num2, double *cosTheta2, double *phi2, int numBins, long int **hist, int numTheta, int distanceInLog);
void calc_distance_distribution_DD_cosTheta_over_ranks(int num1, double *cosTheta1, double *phi1, int num2, double *cosTheta2, double *phi2, int numBins, long int **hist, int numTheta, int distanceInLog);

/* --------------------------------------------------------- */
/* COMMUNICATING (& COMPUTING) DISTANCE HISTOGRAM              */
/* --------------------------------------------------------- */

void count_pairs_r(domain_t *thisDomain, int numSources, double *xsource, double *ysource, double *zsource, int numRand, double *xrand, double *yrand, double *zrand, int numBins, long int **histDD, long int **histDR, long int **histRR, int distanceInLog);
void count_pairs_cosTheta(domain_t *thisDomain, int numSources, double *cosTheta, double *phi, int numRand, double *cosThetaRand, double *phiRand, int numBins, long int **histDD, long int **histDR, long int **histRR, int numTheta, int distanceInLog);

#ifdef __MPI
void collapse_histogram_across_ranks(int numBins, long int **hist);
void collapse_pairs_across_ranks(int numBins, long int **histDD, long int **histDR, long int **histRR);
void collapse_numPairs_across_ranks(int *numSources, int *numRand);
#endif

/* --------------------------------------------------------- */
/* CALCULATING CORRELATION FUNCTION                          */
/* --------------------------------------------------------- */

/* r bins need to range from 0 to 1 */
double *calc_r_bins(int numBins, int distanceInLog);
/* theta bins need to range from 0 to 180 degrees */
double *calc_theta_bins(int numBins, int distanceInLog, int numTheta);   // gives result in degrees
void generate_histograms_r(int numBins, long int **histDD, long int **histDR, long int **histRR, double **distanceArray, int distanceInLog);
void generate_histograms_theta(int numBins, long int **histDD, long int **histDR, long int **histRR, double **distanceArray, int distanceInLog, int numTheta);
void deallocate_histograms(double *distanceArray, long int *histDD, long int *histDR, long int *histRR);
double *calc_correlation_landy_szalay(int numBins, long int **histDD, long int **histDR, long int **histRR, int numData, int numRand);
void write_corrfunc_to_file(domain_t *thisDomain, int numBins, double *distanceArray, double *correlation, char *type, char *filename);

/* --------------------------------------------------------- */

void calc_3Dcorrfunc_cartesian(int numSources, double *xsource, double *ysource, double *zsource, boxparams_t *thisBoxparams, domain_t *thisDomain, char *filename);
void calc_ACF_cartesian(int numSources, double *xsource, double *ysource, double *zsource, boxparams_t *thisBoxparams, cosmparams_t *thisCosmparams, domain_t *thisDomain, char *filename);
void calc_ACF_ra_dec(int numSources, double *cosTheta, double *phi, boxparams_t *thisBoxparams, cosmparams_t *thisCosmparams, domain_t *thisDomain, char *filename);


#endif
