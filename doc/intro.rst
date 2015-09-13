[2]ifpackageloaded#1#2 [2]ifpackageloaded#1#2 [3]ifpackageloaded#1#2#3

#1

Introduction
============

A popular software library for computing FFTs is FFTW . This library
also includes a parallel FFT implementation (FFTW-MPI) based on the
Message Passing Interface (MPI). FFTW-MPI parallelizes multi-dimensional
FFTs by a mixture of serial lower-dimensional FFTs and parallel data
transpositions. However, FFTW-MPI makes use of a one-dimensional data
decomposition, which shows to be a scalability bottleneck on large
scale, parallel computers. For example, a three-dimensional FFT of size
:math:`1024^3` can be computed with at most :math:`1024` MPI processes.
In contrast, using a two-dimensional data decomposition would increase
the maximum number of MPI processes to :math:`1024^2` in this case.

The main goal of PFFT is to extend the MPI part of the FFTW software
library to multi-dimensional data decompositions, i.e.,
:math:`d`-dimensional FFTs of size :math:`N^d` can be computed in
parallel with at most :math:`N^{d-1}` MPI processes. In addition, PFFT
offers several extra features that are particular usefull for parallel,
distributed memory FFTs but are not yet present in FFTW-MPI. We refer to
the publication  for a closer look on the different data decompositions
and the underlying algorithms of the PFFT library.

The interface of PFFT is as close as possible to the FFTW-MPI interface.
In fact, we consider every difference between PFFT and FFTW that is not
explicitly mentioned within this manual as a bug that should be reported
to https://github.com/mpip/pfft.git. Therefore, porting code that uses
FFTW-MPI to PFFT is almost trivial, e.g. see Section [sec:porting].

Most features of PFFT are inherited from FFTW or similarily implemented.
These include the following:

We employ fast :math:`\mathcal{O}(N\log N)` algorithms of FFTW to
compute arbitrary-size discrete Fourier transforms of complex data, real
data, and even- or odd-symmetric real data.

The dimension of the FFT can be arbitrary. However, parallel data
decomposition must be at least one dimension smaller.

PFFT offers portable performance; e.g., it will perform well on most
platforms.

The application of PFFT is split into a time consuming planning step and
a high performance execution step.

Installing the library is easy. It is based on the common sequence of
configure, make, and make install.

The interface of PFFT is very close to the MPI interface of FFTW. In
fact, we tried to add as few extra parameters as possible.

PFFT is written in C but also offers a Fortran interface, see
Section [sec:fortran].

FFTW includes shared memory parallelism for all serial transforms. This
enables us to benefit from hybrid parallelism to a certain amount, see
Section [sec:openmp].

All steps of our parallel FFT can be performed completely in place. This
is especially remarkable for the global transposition routines.

Confirming to good MPI programming practice, all PFFT transforms can be
performed on user defined communicators. In other words, PFFT does not
enforce the user to work with ``MPI_COMM_WORLD``.

PFFT uses the same algorithm to compute the size of the local array
blocks as FFTW. This implies that the FFT size need not be divisible by
the number of processes.

PFFT supports single, double and long double precision.

PFFT supports new-array execution, i.e., a PFFT plan can be planned and
executed on different plans up to some restrictions, see
Section [sec:new-array] for details. Thanks to Yu Feng for the new-array
execute patch.

Furthermore, we added some special features to support repeated tasks
that often occur in practical application of parallel FFTs.

PFFT includes a very flexible ghost cell exchange module. A detailed
description of this module is given in Section [sec:gc].

PFFT accepts three-dimensional data decomposition even for
three-dimensional FFTs. However, the underlying parallel FFT framework
is still based on two-dimensional decomposition. A more detailed
description can be found in Section [sec:3don2d].

PFFT explicitly supports the parallel calculation of pruned FFTs.
Details are given in Section [sec:pruned].

Finally, we complete this overview with a list of features that are (not
yet) implemented in PFFT.

Parallel one-dimensional FFT based on MPI. FFTW-MPI uses another
parallelization strategy for one-dimensional FFTs, which is not
implemented in PFFT. The reason is that we can not achive a scalability
benefit due to higher dimensional data decomposition if the FFT has only
one dimension. Therefore, one can also call FFTW directly in this case.

There is no equivalent of FFTW *wisdom* in PFFT, i.e., you can not save
a PFFT plan to disk and restore it for later use.

PFFT does not have full OpenMP support. All serial FFT computations and
global communications are implemented with FFTW, which offers OpenMP
support, see Section [sec:openmp]. However, most of the PFFT-only
features, such as pruned FFT, ghost cell send and 3d decompostion of 3d
FFTs are not yet parallelized with OpenMP.

PFFT does not have full SIMD support. All serial FFT computations and
global communications are implemented with FFTW, which offers SIMD
support, see Section [sec:simd]. However, most of the PFFT-only
features, such as pruned FFT, ghost cell send and 3d decompostion of 3d
FFTs are not yet parallelized with SIMD.

PFFT does not overlap communication and computation. The code of PFFT is
build in a very modularized structure. Most of these modules consist of
FFTWs routines. Therefore, the global transposition does not support non
blocking communication.

Similar to FFTW, we do not provide any parallel IO routines. The user is
responsible of load and store of parallel data.

PFFT depends on FFTW to perform its serial transforms and does not
support different vendor FFTs (such as Intel’s MKL or IBM’s ESSL).
However, this is not assumed to be a big drawback, since FFTW seems to
perform very well on most platforms.

The global communication routines can not be called separately. However,
it should be possible to implement a user interface to our global
transposition routines.

PFFT does not support GPU parallelization.

You are welcome to propose new PFFT features at
https://github.com/mpip/pfft.git.

Alternative parallel FFT implementations
----------------------------------------

There have been several FFT implementations that aim to circumvent the
scalability bottleneck for at least three dimensional FFTs by using
two-dimensional decomposition approach. However, these implementations
are often fitted to special problems and where not published as a stand
alone software library. Remarkable exceptions are the parallel FFT
software library by S. Plimpton , the P3DFFT software library by
D. Pekurovsky  and the software library by N. Li .

Parallel nonequispaced FFT
--------------------------

If your are interested in a parallel implementation of nonequispaced
fast Fourier transforms (NFFT) for distributed memory architectures, you
should have a look at our PNFFT software library  that is also available
at https://github.com/mpip/pnfft.git.
