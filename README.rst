Description
===========

A library to compute the 3D and angular correlation function using MPI

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
