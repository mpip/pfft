[2]ifpackageloaded#1#2 [2]ifpackageloaded#1#2 [3]ifpackageloaded#1#2#3

#1

Installation and linking
========================

The install of PFFT is based on the Autotools and follows the typical
workflow

::

    ./configure
    make
    make install

Install of the latest official FFTW release
-------------------------------------------

PFFT depends on Release 3.3.3 of the FFTW library . For the sake of
completeness, we show the command line based install procedure in the
following. However, note that we provide install scripts on
`{www.tu-chemnitz.de/~mpip}/software.php <{www.tu-chemnitz.de/~mpip}/software.php>`__that
simplify the install a lot. We highly recommend to use these install
scripts, since they additionally apply several performance patches and
bugfixes that have been submitted to the FFTW developers but are not yet
included in the official FFTW releases.

::

    wget http://www.fftw.org/fftw-§\fftwversionsl§.tar.gz
    tar xzvf fftw-§\fftwversion§.tar.gz
    cd fftw-§\fftwversion§
    ./configure --enable-mpi --prefix=$HOME/local/fftw3_mpi §\label{lst:fftw:conf}§
    make
    make install

The MPI algorithms of FFTW must be build with a MPI C compiler. Add the
statement ``MPICC=\$MPICCOMP`` at the end of line [lst:fftw:conf] if the
``configure`` script fails to determine the right MPI C compiler
``\$MPICCOMP``. Similarly, the MPI Fortran compiler ``\$MPIFCOMP`` is
set by ``MPIFC=\$MPIFCOMP``.

Install of the PFFT library
---------------------------

In the simplest case, the hardware platform and the -3.3.3 library are
recognized by the PFFT configure script automatically, so all we have to
do is

::

    wget http://www.tu-chemnitz.de/~mpip/software/pfft-§\pfftversionsl§.tar.gz
    tar xzvf pfft-§\pfftversion§.tar.gz
    cd pfft-§\pfftversion§
    ./configure
    make
    make check
    make install

Hereby, the optional call ``make check`` builds the test programs. If
the -3.3.3 software library is already installed on your system but not
found by the configure script, you can provide the FFTW installation
directory ``\$FFTWDIR`` to configure by

.. code:: bash

    ./configure --with-fftw3=$FFTWDIR

This call implies that the FFTW header files are located in
``\$FFTWDIR/include`` and the FFTW library files are located in
``\$FFTWDIR/lib``. Otherwise, one should specify the FFTW include path
``\$FFTWINC`` and the FFTW library path ``\$FFTWLIB`` separately by

::

    ./configure --with-fftw3-includedir=$FFTWINC --with-fftw3-libdir=$FFTWLIB

At the end, this is equivalent to

::

    ./configure CPPFLAGS=-I$FFTWINC LDFLAGS=-L$FFTWLIB

which is more common to experienced users of the Autotools. To install
PFFT in a user specified directory ``\$PFFTINSTDIR`` call configure with
the option

::

    ./configure --prefix=$PFFTINSTDIR

However, this option is mandatory whenever you do not have root
permissions on your machine, since the default install paths of
``configure`` are not accessible by standard users. The PFFT library
must be built with a MPI compiler. In Section [sec:fftw\ :sub:`i`\ nst]
we already described how to hand the right compilers to the
``configure`` script. Some more options are

:code:\`[\`keywords=]–enable-float: Produces a single-precision version
of PFFT (float) instead of the default double-precision (double); see
[sec:prec].

:code:\`[\`keywords=]–enable-long-double: Produces a long-double
precision version of PFFT (long double) instead of the default
double-precision (double); see [sec:prec].

``--disable-fortran``: Disables inclusion of Fortran wrapper routines in
the standard PFFT libraries.

``--disable-tests``: Disables build of test programs.

For more details on the options of the ``configure`` script call

::

    ./configure --help

How to include PFFT in your program
-----------------------------------

All programs using PFFT should include its header file

::

    #include <pfft.h>

This header includes the FFTW headers ``fftw.h``, ``fftw-mpi.h``
automatically. Make sure that the compiler can find them by setting the
include flags appropriately. You must also link to the PFFT, FFTW and
FFTW-MPI libraries. On Unix, this means adding
``-lpfft -lfftw3_mpi -lfftw3 -lm`` at the end of the link command. For
example, to build ``pfft_test.c`` use the following compiler invocation

::

    mpicc pfft_test.c -I$PFFTINC -I$FFTWINC -L$PFFTLIB -L$FFTWLIB -lpfft -lfftw3_mpi -lfftw3 -lm

Substitute ``mpicc`` by any other MPI C compiler if you like.
``\$PFFTINC``, ``\$FFTWINC``, ``\$PFFTLIB``, and ``\$FFTWLIB`` denote
the PFFT and FFTW include and library paths, respectively. If you use
the install scripts mentioned in Sect. [sec:pfft-inst], these paths will
be

::

    PFFTINC = $HOME/local/pfft-§\pfftversion§/include
    FFTWINC = $HOME/local/fftw-§\fftwversion§/include
    PFFTINC = $HOME/local/pfft-§\pfftversion§/lib
    FFTWINC = $HOME/local/fftw-§\fftwversion§/lib

