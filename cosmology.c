#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "phys_const.h"
#include "cosmology.h"

#define SQR(X) ((X) * (X))
#define CUB(X) ((X) * (X) * (X))
#define M_PI acos(-1.0)

double time_from_redshift_flatuniverse(double h, double omega_l, double omega_m, double zmin, double zmax)
{
    double H0 = h * 1.e7 / Mpc_cm;
    double prefactor = 2. / (3 * H0 * sqrt(omega_l));
    double tmp = sqrt(omega_l / omega_m);
    
    return prefactor*(asinh(tmp*pow(1.+zmin, -1.5)) - asinh(tmp*pow(1.+zmax, -1.5)));
}

double hubble_flatuniverse(double h, double omega_l, double omega_m, double z)
{
    double H0 = h * 1.e7 / Mpc_cm;
        
    return H0 * pow(omega_m * pow((1.+z),3) + omega_l,0.5);;
}


double mean_numdensity_flatuniverse(double h, double omega_b, double z)
{
    double H0 = h * 1.e7 / Mpc_cm;
    double mean_numdensity = 3.*SQR(H0) / (8.*M_PI*G) / mp_g * omega_b * CUB(1.+z);
    
    return mean_numdensity;
}

double mean_numdensity_H_flatuniverse(double h, double omega_b, double Y, double z)
{
    double H0 = h * 1.e7 / Mpc_cm;
    double mean_numdensity_H = 3.*SQR(H0) / (8.*M_PI*G) / mp_g * omega_b * CUB(1.+z) * (1.-Y);

    return mean_numdensity_H;
}

double mean_numdensity_He_flatuniverse(double h, double omega_b, double Y, double z)
{
    double H0 = h * 1.e7 / Mpc_cm;
    double mean_numdensity_He = 3.*SQR(H0)/(8.*M_PI*G) / mp_g * omega_b * CUB(1.+z) * 0.25*Y;

    return mean_numdensity_He;
}


double calcComDist(double h, double omega_m, double omega_l, double zmax, double zmin)  // in cm
{
    double z, delta_z;
    int N = 100;
    double sum = 0.;

    delta_z = (zmax - zmin) / (double)N;
    if(delta_z > 1.e-5)
    {
        N = (int)((zmax - zmin)*1.e5);
        delta_z = (zmax-zmin) / (double)N;
    }
    
    z = zmax;
//     printf("H(z=%e) = %e\n", z, hubble_flatuniverse(h, omega_m, omega_l, z));
    for(int i=0; i<N; i++){
        sum += clight_cm / hubble_flatuniverse(h, omega_m, omega_l, z) * delta_z;
        z = z - delta_z;
    }
    
    return sum;
}

double calcPhysDist(double h, double omega_m, double omega_l, double zmax, double zmin) // in cm
{
    double z, delta_z;
    int N = 100;
    double sum = 0.;

    delta_z = (zmax - zmin) / (double)N;
    if(delta_z > 1.e-5)
    {
        N = (int)((zmax - zmin)*1.e5);
        delta_z = (zmax-zmin) / (double)N;
    }
    
    z = zmax;
    for(int i=0; i<N; i++){
        sum += clight_cm /(hubble_flatuniverse(h, omega_m, omega_l, z) * (1.+z)) * delta_z;
        z = z - delta_z;
    }
    
    return sum;
}

double calcLumDist(double h, double omega_m, double omega_l, double zmax, double zmin)
{
    double z, delta_z;
    int N = 100;
    double sum = 0.;

    delta_z = (zmax - zmin) / (double)N;
    if(delta_z > 1.e-5)
    {
        N = (int)((zmax - zmin)*1.e5);
        delta_z = (zmax-zmin) / (double)N;
    }
    
    z = zmax;
    for(int i=0; i<N; i++){
        sum += clight_cm /hubble_flatuniverse(h, omega_m, omega_l, z) * delta_z;
        z = z - delta_z;
    }
    sum = (1. + zmax) * sum;

    return sum;
}
