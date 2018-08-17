#ifndef COSMOLOGY_H
#define COSMOLOGY_H

double time_from_redshift_flatuniverse(double h, double omega_l, double omega_m, double zmin, double zmax);

double hubble_flatuniverse(double h, double omega_l, double omega_m, double z);

double mean_numdensity_flatuniverse(double h, double omega_b, double z);
double mean_numdensity_H_flatuniverse(double h, double omega_b, double Y, double z);
double mean_numdensity_He_flatuniverse(double h, double omega_b, double Y, double z);

double calcComDist(double h, double omega_m, double omega_l, double zmax, double zmin);
double calcPhysDist(double h, double omega_m, double omega_l, double zmax, double zmin);
double calcLumDist(double h, double omega_m, double omega_l, double zmax, double zmin);

#endif
