Description
===========

A library to compute the 3D and angular correlation function using MPI. Data and random points are distributed over domains, and pairs are calculated within and across domains.
A fast serial implementation can be found on the `Corrfunc webpage <https://github.com/manodeep/Corrfunc>`__

Installation
============

Pre-requisities
---------------

Serial run
``````````

1. fftw3 library: ``fftw3 >= 3.3.3``

Parallel run
````````````

1. MPI library
2. fftw3 & fftw3-mpi library: ``fftw3 >= 3.3.3``

FFTW3
'''''

Go to the `FFTW webpage <http://www.fftw.org/download.html>` to install fftw3. Ensure to compile the library with the ``enable-mpi`` flag for parallel runs
::
    
    $ ./configure --enable-mpi
    $ make
    $ make install
    
Note: To create the dynamic libraries, run configure with the ``--enable-shared`` flag. 

Download & Build
----------------

::

    $ git clone https://github.com/annehutter/pcorrfunc.git
    $ make

This will download the code and first test case from the github directory and compile the source code.

Execution
---------

The first test case can then be run by
::

    $ mpiexec -np #PROCESSORS ./pcorrfunc

    
Functions
=========

The library contains currently three functions:
- **calc_3Dcorrfunc_cartesian** *(int numPoints, double *x, double *y, double *z, boxparams_t *thisBoxparams, domain_t *thisDomain, char *outputName)*: computes the 3D correlation function from a list of x, y, z coordinates.
- **calc_ACF_cartesian** *(int numPoints, double *x, double *y, double *z, boxparams_t *thisBoxparams, cosmparams_t *thisCosmparams, domain_t *thisDomain, char *outputName)*: computes the angular correlation function for the box or a set of points located at redshift z from a list of x, y, z coordinates.
- **calc_ACF_ra_dec** *(int numPoints, double *cosTheta, double *phi, boxparams_t *thisBoxparams, cosmparams_t *thisCosmparams, domain_t *thisDomain, char *outputName)*: computes the angular correlation function for the box or a set of points located at redshift z from a list of cos(theta), phi coordinates.

These functions require as input the number of points, *double* arrays of the *x*, *y*, *z* coordinates, and the following structs:
- *boxparams_t*: contains the box measurements, i.e. also galaxies in a subbox can be chosen.
- *cosmparams_t*: contains the cosmological parameters and redshift where the box or set of point are located. This is only required for the angular correlation function.
- *domain_t*: contains the domain decomposition. Currently the domain decomposition agrees with those in FFTW, i.e. slicing the box along the z-axis.
Examples how these structs are initialised and used are in the example ``main.c``.
