PFFT - Massively Parallel FFT based on FFTW3
============================================

.. image:: https://api.travis-ci.org/mpip/pfft.svg
    :alt: Build Status
    :target: https://travis-ci.org/mpip/pfft/

Overview
--------
PFFT is a software library for computing massively parallel, fast Fourier transformations
on distributed memory architectures. PFFT can be understood as a generalization of FFTW-MPI
to multidimensional data decomposition.
The library is written in C and MPI. A Fortran interface is also available.
Support for hybrid parallelization based on OpenMP and MPI is under development.

Install
-------
At first, you need an install of FFTW-3.3 with enabled MPI support. 
We highly recommend to install the latest release of FFTW since the MPI code got several
bugfixes in the last releases. Since some fixes are still pending until release FFTW-3.3.5.
we offer some helpful scripts for installing FFTW-3.3.4 together with two patches
(one bugfix and one performance improvement) at our sofware page

  https://www-user.tu-chemnitz.de/~potts/workgroup/pippig/software.php#scripts

The install of PFFT follows the typical steps::

  ./bootstrap.sh
  ./configure
  make
  make install

Optionally, a bunch of test programs can be built with::

  make check

The bootstrap step can be skipped if you delivered a PFFT tarball,
i.e., the file configure was already generated.
Make sure that configure can find a working install of FFTW-3.3 with MPI support,
i.e., set::

  CPPFLAGS=$PATHTOFFTW/include

and::

  LDFLAGS=$PATHTOFFTW/lib64 or LDFLAGS=$PATHTOFFTW/lib

depending on your system architecture.

Documentation
-------------

PFFT tarballs include a detailed user guide at::

  doc/manual.pdf

If you have cloned the PFFT repository, the user manual is automatically built during `make`.
Of course this requires a working LaTeX enviroment.
Alternatively, you can download a recent version of the user manual at

  http://www.tu-chemnitz.de/~potts/workgroup/pippig/software.php.en

Note that using PFFT is very similar to FFTW. The interface is as close to the fftw_mpi
interface as possible. Therefore, it is a good start to read FFTW manual:

  http://www.fftw.org/fftw3_doc

At least you should understand how FFTW deals with distributed memory FFTs:

  http://www.fftw.org/fftw3_doc/Distributed_002dmemory-FFTW-with-MPI.html#Distributed_002dmemory-FFTW-with-MPI

Next, you can have a look at the test programs in directory 'tests' to learn the details of the PFFT interface.

For installation instructions, you can also refer to the file INSTALL
in this directory.

For an theoretical introduction, please read the paper

  ''PFFT - An Extension of FFTW to Massively Parallel Architectures''

available at

  http://www.tu-chemnitz.de/~potts/workgroup/pippig/publikationen.php.en

This is the most current general paper, and the one that we recommend if you wish 
to cite PFFT.

Python interface
----------------

A python interface is available at a distinct repository

  https://github.com/rainwoodman/pfft-python

Many thanks to Yu Feng for his great work.


Directory structure
-------------------

==================  ============================================================
aclocal.m4          Macros for configure script                                 
------------------  ------------------------------------------------------------
api (dir)           Source code for user interface                              
------------------  ------------------------------------------------------------
AUTHORS             Information about the authors of PFFT                       
------------------  ------------------------------------------------------------
bootstrap.sh        Bootstrap shell script that call Autoconf and friends       
                    in order to generate configure                              
------------------  ------------------------------------------------------------
build-aux (dir)     Used by configure script                                    
------------------  ------------------------------------------------------------
ChangeLog           A short version history                                     
------------------  ------------------------------------------------------------
config.h.in         Used by configure script                                    
------------------  ------------------------------------------------------------
configure           Configure script build from configure.ac by bootstrap.sh    
------------------  ------------------------------------------------------------
configure.in        Autoconf configure script template                          
------------------  ------------------------------------------------------------
CONVENTIONS         Makro naming conventions for developers                     
------------------  ------------------------------------------------------------
COPYING             Information about redistributing PFFT                       
------------------  ------------------------------------------------------------
doc (dir)           User and developer documentation                            
------------------  ------------------------------------------------------------
fconfig.h.in        Used by configure script (Fortran definitions)              
------------------  ------------------------------------------------------------
gcell (dir)         Source code for ghost cell support                          
------------------  ------------------------------------------------------------
include (dir)       Header files                                                
------------------  ------------------------------------------------------------
INSTALL             Installation instructions                                   
------------------  ------------------------------------------------------------
kernel (dir)        Source code for core library routines                       
------------------  ------------------------------------------------------------
m4 (dir)            Contains macros for configure script                        
------------------  ------------------------------------------------------------
Makefile.am         Automake Makefile template                                  
------------------  ------------------------------------------------------------
Makefile.in         Makefile template generated from Makefile.am,               
                    processed by configure script                               
------------------  ------------------------------------------------------------
NEWS                New and noteworthy                                          
------------------  ------------------------------------------------------------
pfft.pc.in          Template for PFFT package information                       
------------------  ------------------------------------------------------------
README              This file                                                   
------------------  ------------------------------------------------------------
scripts (dir)       A collection of useful script files                         
------------------  ------------------------------------------------------------
tests (dir)         Simples examples for using PFFT routines                    
------------------  ------------------------------------------------------------
TODO                Current work to be done                                     
------------------  ------------------------------------------------------------
util (dir)          Source code for auxilliary routines                         
==================  ============================================================

Feedback
--------
Your comments are welcome! This is the first version of the library and may
not be as robust or well documented as it should be. Please keep track of bugs
or missing/confusing instructions and report them to

  Michael Pippig <michael.pippig.tuc@gmail.com>

If you find PFFT useful, we would be delighted to hear about what application
you are using PFFT for!
 
