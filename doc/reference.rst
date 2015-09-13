[2]ifpackageloaded#1#2 [2]ifpackageloaded#1#2 [3]ifpackageloaded#1#2#3

#1

PFFT Reference
==============

Files and Data Types
--------------------

You must include the PFFT header file by

::

    #include <pfft.h>

in the preamble of each source file that calls PFFT. This header
automatically includes ``fftw.h`` and ``fftw3-mpi.h``. Therefore, PFFT
can use the ``fftw_complex`` data type defined in ``fftw.h``, see . Note
that ``fftw_complex`` is defined to be the C99 native complex whenever
``<complex.h>`` is included *before* ``<fftw.h>``, ``<fftw-mpi.h>`` and
``<pfft.h>``. Otherwise it is defined as

::

    typedef double fftw_complex[2];

For the sake of a clean namespace we define the wrapper data type
``pfft_complex`` as

::

    typedef fftw_complex pfft_complex;

that can be used equivallently to ``fftw_complex``. Futhermore, we
define the wrapper functions

::

    void *pfft_malloc(size_t n);
    double *pfft_alloc_real(size_t n);
    pfft_complex *pfft_alloc_complex(size_t n);
    void pfft_free(void *p);

as substitues for their corresponding FFTW equivalents, see . Note that
memory allocated by one of these functions must be freed with
``pfft_free`` (or its equivalent ``fftw_free``). Because of the
performance reasons given in  we recommend to use one of the ``pfft_``
(or its equivalent ``fftw_``) allocation functions for all arrays
containing FFT inputs and outputs. However, PFFT will also work
(possibly slower) with any other memory allocation method.

Different precisions are handled as in FFTW: That is ``pfft_`` functions
and datatypes become ``pfftf_`` (single precision) or ``pfftl_`` (long
double precision) prefixed. Quadruple precision is not yet supported.
The main problem is that we do not know about a suitable MPI datatype to
represent ``__float128``.

MPI Initialization
------------------

Initialization and cleanup of PFFT in done in the same way as for
FFTW-MPI, see . In order to keep a clean name space, PFFT offers the
wrapper functions

::

    void pfft_init(void);
    void pfft_cleanup(void);

that can be used as substitutes for ``fftw_mpi_init`` and
``fftw_mpi_cleanup``, respectively.

Using PFFT Plans
----------------

PFFT follows exactly the same workflow as FFTW-MPI. A plan created by
one of the functions given in Section [sec:create-plan] is executed with

::

    void pfft_execute(const pfft_plan plan);

and freed with

::

    void pfft_destroy_plan(const pfft_plan plan);

Note, that you can *not* apply ``fftw_mpi_execute`` or ``fftw_destroy``
on PFFT plans.

The new array execute functions are given by

::

    void pfft_execute_dft(const pfft_plan plan, pfft_complex *in, pfft_complex *out);
    void pfft_execute_dft_r2c(const pfft_plan plan, double *in, pfft_complex *out);
    void pfft_execute_dft_c2r(const pfft_plan plan, pfft_complex *in, double *out);
    void pfft_execute_r2r(const pfft_plan plan, double *in, double *out);

The arrays given by ``in`` and ``out`` must have the correct size and
the same alignement as the array that were used to create the plan, just
as it is the case for FFTW, see [fftw-new-array].

Data Distribution Functions
---------------------------

Complex-to-Complex FFT
~~~~~~~~~~~~~~~~~~~~~~

::

    ptrdiff_t pfft_local_size_dft_3d(
        const ptrdiff_t *n, MPI_Comm comm_cart, unsigned pfft_flags,
        ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
        ptrdiff_t *local_no, ptrdiff_t *local_o_start);
    ptrdiff_t pfft_local_size_dft(
        int rnk_n, const ptrdiff_t *n,
        MPI_Comm comm_cart, unsigned pfft_flags,
        ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
        ptrdiff_t *local_no, ptrdiff_t *local_o_start);
    ptrdiff_t pfft_local_size_many_dft(
        int rnk_n, const ptrdiff_t *n, const ptrdiff_t *ni, const ptrdiff_t *no,
        ptrdiff_t howmany, const ptrdiff_t *iblock, const ptrdiff_t *oblock,
        MPI_Comm comm_cart, unsigned pfft_flags,
        ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
        ptrdiff_t *local_no, ptrdiff_t *local_o_start);

Compute the data distribution of a parallel, complex input/output
discrete Fourier transform (DFT) in two or more dimensions, returning
the number of *complex* numbers that must be allocated to hold the
parallel transform.

Arguments:

``rnk_n`` is the rank of the transform (typically the size of the arrays
``n``, ``ni``, ``no``) that can be any integer :math:`\ge 2`. The
``_3d`` planner corresponds to a ``rnk_n`` of 3.

The array ``n`` of size ``rnk_n`` specifies the transform dimensions.
They can be any positive integer.

The array ``ni`` of size ``rnk_n`` specifies the input array dimensions.
They can be any positive integer with ``ni[t] <= n[t]`` for all
dimensions ``t=0,...,rnk_n-1``. For ``ni[t]<n[t]`` the inputs will be
padded with zeros up to size ``n[t]`` along the ``t``-th dimension
before the transform, see Section [sec:pruned-fft].

The array ``no`` of size ``rnk_n`` specifies the output array
dimensions. They can be any positive integer with ``no[t] <= n[t]`` for
all dimensions ``t=0,...,rnk_n-1``. For ``no[t]<n[t]`` the outputs will
be pruned to size ``no[t]`` along the ``t``-th dimension after the
transform, see Section [sec:pruned-fft].

``howmany`` is the number of transforms to compute. The resulting plan
computes howmany transforms, where the input of the k-th transform is at
location in+k (in C pointer arithmetic) with stride ``howmany``, and its
output is at location out+k with stride ``howmany``. The basic
``pfft_plan_dft`` interface corresponds to howmany=1.

``comm_cart`` is a Cartesian communicator of dimension ``rnk_pm`` that
specifies the parallel data decomposition, see
Section [sec:data-decomp]. Most of the time, PFFT requires
``rnk_pm < rnk_n``. The only exception is the case
``rnk_pm == rnk_n == 3``, see Section [sec:3don3d]. If an ordinary (i.e.
non-Cartesian) communicator is passed, PFFT internally converts it into
a one-dimensional Cartesian communicator while retaining the MPI ranks
(this results in the FFTW-MPI data decomposition).

The arrays ``iblock`` and ``oblock`` of size ``rnk_pm``\ 1+ specify the
block sizes for the first ``rnk_pm``\ 1+ dimensions of the input and
output data, respectively. These must be the same block sizes as were
passed to the corresponding ``local_size`` function. You can pass
``PFFT_DEFAULT_BLOCKS`` to use PFFT’s default block sizes. Furthermore,
you can use ``PFFT_DEFAULT_BLOCK`` to set the default block size in
separate dimensions, e.g., ``iblock[t]=PFFT_DEFAULT_BLOCK``.

``pfft_flags`` is a bitwise OR (’``|``\ ’) of zero or more planner
flags, as defined in Section [sec:flags].

The array ``local_ni`` of size ``rnk_n`` returns the size of the local
input array block in every dimension (counted in units of complex
numbers).

The array ``local_i_start`` of size ``rnk_n`` returns the offset of the
local input array block in every dimension (counted in units of complex
numbers).

The array ``local_no`` of size ``rnk_n`` returns the size of the local
output array block in every dimension (counted in units of complex
numbers).

The array ``local_o_start`` of size ``rnk_n`` returns the offset of the
local output array block in every dimension (counted in units of complex
numbers).

In addition, the following ``local_block`` functions compute the local
data distribution of the process with MPI rank ``pid``. The
``local_size`` interface can be understood as a call of ``local_block``
where ``pid`` is given by ``MPI_Comm_rank(comm_cart, &pid)``, i.e., each
MPI process computes its own data block. However, ``local_block``
functions have a ``void`` return type, i.e., they omit the computation
of the local array size that is necessary to hold the parallel
transform. This makes ``local_block`` functions substantially faster in
exectuion.

::

    void pfft_local_block_dft_3d(
        const ptrdiff_t *n, MPI_Comm comm_cart, int pid, unsigned pfft_flags,
        ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
        ptrdiff_t *local_no, ptrdiff_t *local_o_start);
    void pfft_local_block_dft(
        int rnk_n, const ptrdiff_t *n,
        MPI_Comm comm_cart, int pid, unsigned pfft_flags,
        ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
        ptrdiff_t *local_no, ptrdiff_t *local_o_start);
    void pfft_local_block_many_dft(
        int rnk_n, const ptrdiff_t *ni, const ptrdiff_t *no,
        const ptrdiff_t *iblock, const ptrdiff_t *oblock,
        MPI_Comm comm_cart, int pid, unsigned pfft_flags,
        ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
        ptrdiff_t *local_no, ptrdiff_t *local_o_start);

Real-to-Complex FFT
~~~~~~~~~~~~~~~~~~~

::

    ptrdiff_t pfft_local_size_dft_r2c_3d(
        const ptrdiff_t *n, MPI_Comm comm_cart, unsigned pfft_flags,
        ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
        ptrdiff_t *local_no, ptrdiff_t *local_o_start);
    ptrdiff_t pfft_local_size_dft_r2c(
        int rnk_n, const ptrdiff_t *n,
        MPI_Comm comm_cart, unsigned pfft_flags,
        ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
        ptrdiff_t *local_no, ptrdiff_t *local_o_start);
    ptrdiff_t pfft_local_size_many_dft_r2c(
        int rnk_n, const ptrdiff_t *n, const ptrdiff_t *ni, const ptrdiff_t *no,
        ptrdiff_t howmany, const ptrdiff_t *iblock, const ptrdiff_t *oblock,
        MPI_Comm comm_cart, unsigned pfft_flags,
        ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
        ptrdiff_t *local_no, ptrdiff_t *local_o_start);

Compute the data distribution of a parallel, real-input/complex-output
discrete Fourier transform (DFT) in two or more dimensions, returning
the number of *complex* numbers that must be allocated to hold the
parallel transform.

Arguments are the same as for c2c transforms (see
Section [sec:local-size-c2c]) with the following exceptions:

The logical input array size ``ni`` will differ from the physical array
size of the real inputs if the flag ``PFFT_PADDED_R2C`` is included in
``pfft_flags``. This results from the padding at the end of the last
dimension that is necessary to align the real valued inputs and complex
valued outputs for inplace transforms, see . In contrast to FFTW-MPI,
PFFT does not pad the r2c inputs per default.

``local_ni`` is counted in units of real numbers. It will include
padding

``local_i_start`` is counted in units of real numbers.

The corresponding ``local_block`` functions compute the local data
distribution of the process with MPI rank ``pid``.

::

    void pfft_local_block_dft_r2c_3d(
        const ptrdiff_t *n, MPI_Comm comm_cart, int pid, unsigned pfft_flags,
        ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
        ptrdiff_t *local_no, ptrdiff_t *local_o_start);
    void pfft_local_block_dft_r2c(
        int rnk_n, const ptrdiff_t *n,
        MPI_Comm comm_cart, int pid, unsigned pfft_flags,
        ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
        ptrdiff_t *local_no, ptrdiff_t *local_o_start);
    void pfft_local_block_many_dft_r2c(
        int rnk_n, const ptrdiff_t *ni, const ptrdiff_t *no,
        const ptrdiff_t *iblock, const ptrdiff_t *oblock,
        MPI_Comm comm_cart, int pid, unsigned pfft_flags,
        ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
        ptrdiff_t *local_no, ptrdiff_t *local_o_start);

Complex-to-Real FFT
~~~~~~~~~~~~~~~~~~~

::

    ptrdiff_t pfft_local_size_dft_c2r_3d(
        const ptrdiff_t *n, MPI_Comm comm_cart, unsigned pfft_flags,
        ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
        ptrdiff_t *local_no, ptrdiff_t *local_o_start);
    ptrdiff_t pfft_local_size_dft_c2r(
        int rnk_n, const ptrdiff_t *n,
        MPI_Comm comm_cart, unsigned pfft_flags,
        ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
        ptrdiff_t *local_no, ptrdiff_t *local_o_start);
    ptrdiff_t pfft_local_size_many_dft_c2r(
        int rnk_n, const ptrdiff_t *n, const ptrdiff_t *ni, const ptrdiff_t *no,
        ptrdiff_t howmany, const ptrdiff_t *iblock, const ptrdiff_t *oblock,
        MPI_Comm comm_cart, unsigned pfft_flags,
        ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
        ptrdiff_t *local_no, ptrdiff_t *local_o_start);

Compute the data distribution of a parallel, complex-input/real-output
discrete Fourier transform (DFT) in two or more dimensions, returning
the number of *complex* numbers that must be allocated to hold the
parallel transform.

Arguments are the same as for c2c transforms (see
Section [sec:local-size-c2c]) with the following exceptions:

The logical output array size ``no`` will differ from the physical array
size of the real outputs if the flag ``PFFT_PADDED_C2R`` is included in
``pfft_flags``. This results from the padding at the end of the last
dimension that is necessary to align the real valued outputs and complex
valued inputs for inplace transforms, see . In contrast to FFTW-MPI,
PFFT does not pad the c2r outputs per default.

``local_no`` is counted in units of real numbers.

``local_o_start`` is counted in units of real numbers.

The corresponding ``local_block`` functions compute the local data
distribution of the process with MPI rank ``pid``.

::

    void pfft_local_block_dft_c2r_3d(
        const ptrdiff_t *n, MPI_Comm comm_cart, int pid, unsigned pfft_flags,
        ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
        ptrdiff_t *local_no, ptrdiff_t *local_o_start);
    void pfft_local_block_dft_c2r(
        int rnk_n, const ptrdiff_t *n,
        MPI_Comm comm_cart, int pid, unsigned pfft_flags,
        ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
        ptrdiff_t *local_no, ptrdiff_t *local_o_start);
    void pfft_local_block_many_dft_c2r(
        int rnk_n, const ptrdiff_t *ni, const ptrdiff_t *no,
        const ptrdiff_t *iblock, const ptrdiff_t *oblock,
        MPI_Comm comm_cart, int pid, unsigned pfft_flags,
        ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
        ptrdiff_t *local_no, ptrdiff_t *local_o_start);

Real-to-Real FFT
~~~~~~~~~~~~~~~~

::

    ptrdiff_t pfft_local_size_r2r_3d(
        const ptrdiff_t *n, MPI_Comm comm_cart, unsigned pfft_flags,
        ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
        ptrdiff_t *local_no, ptrdiff_t *local_o_start);
    ptrdiff_t pfft_local_size_r2r(
        int rnk_n, const ptrdiff_t *n,
        MPI_Comm comm_cart, unsigned pfft_flags,
        ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
        ptrdiff_t *local_no, ptrdiff_t *local_o_start);
    ptrdiff_t pfft_local_size_many_r2r(
        int rnk_n, const ptrdiff_t *n, const ptrdiff_t *ni, const ptrdiff_t *no,
        ptrdiff_t howmany, const ptrdiff_t *iblock, const ptrdiff_t *oblock,
        MPI_Comm comm_cart, unsigned pfft_flags,
        ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
        ptrdiff_t *local_no, ptrdiff_t *local_o_start);

Compute the data distribution of a parallel, complex input/output
discrete Fourier transform (DFT) in two or more dimensions, returning
the number of *real* numbers that must be allocated to hold the parallel
transform.

Arguments are the same as for c2c transforms (see
Section [sec:local-size-c2c]) with the following exceptions:

``local_ni`` is counted in units of real numbers.

``local_i_start`` is counted in units of real numbers.

``local_no`` is counted in units of real numbers.

``local_o_start`` is counted in units of real numbers.

The corresponding ``local_block`` functions compute the local data
distribution of the process with MPI rank ``pid``.

::

    void pfft_local_block_r2r_3d(
        const ptrdiff_t *n, MPI_Comm comm_cart, int pid, unsigned pfft_flags,
        ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
        ptrdiff_t *local_no, ptrdiff_t *local_o_start);
    void pfft_local_block_r2r(
        int rnk_n, const ptrdiff_t *n,
        MPI_Comm comm_cart, int pid, unsigned pfft_flags,
        ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
        ptrdiff_t *local_no, ptrdiff_t *local_o_start);
    void pfft_local_block_many_r2r(
        int rnk_n, const ptrdiff_t *ni, const ptrdiff_t *no,
        const ptrdiff_t *iblock, const ptrdiff_t *oblock,
        MPI_Comm comm_cart, int pid, unsigned pfft_flags,
        ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
        ptrdiff_t *local_no, ptrdiff_t *local_o_start);

Plan Creation
-------------

Complex-to-Complex FFT
~~~~~~~~~~~~~~~~~~~~~~

::

    pfft_plan pfft_plan_dft_3d(
        const ptrdiff_t *n, pfft_complex *in, pfft_complex *out, MPI_Comm comm_cart,
        int sign, unsigned pfft_flags);
    pfft_plan pfft_plan_dft(
        int rnk_n, const ptrdiff_t *n, pfft_complex *in, pfft_complex *out, MPI_Comm comm_cart,
        int sign, unsigned pfft_flags);
    pfft_plan pfft_plan_many_dft(
        int rnk_n, const ptrdiff_t *n, const ptrdiff_t *ni, const ptrdiff_t *no,
        ptrdiff_t howmany, const ptrdiff_t *iblock, const ptrdiff_t *oblock,
        pfft_complex *in, pfft_complex *out, MPI_Comm comm_cart,
        int sign, unsigned pfft_flags);
    pfft_plan pfft_plan_many_dft_skipped(
        int rnk_n, const ptrdiff_t *n, const ptrdiff_t *ni, const ptrdiff_t *no,
        ptrdiff_t howmany, const ptrdiff_t *iblock, const ptrdiff_t *oblock,
        const int *skip_trafos, pfft_complex *in, pfft_complex *out, MPI_Comm comm_cart,
        int sign, unsigned pfft_flags);

Plan a parallel, complex input/output discrete Fourier transform (DFT)
in two or more dimensions, returning an ``pfft_plan``. The planner
returns NULL if the plan cannot be created.

Arguments:

``rnk_n``, ``n``, ``ni``, ``no``, ``howmany``, ``iblock``, ``oblock``,
``comm_cart`` must be the same as passed to the corresponding
``pfft_local_size_dft`` function, see Section [sec:local-size-c2c].

The array ``skip_trafos`` of size ``rnk_pm``\ 1+ specifies the serial
transforms that will be omitted. For ``t=0,...,rnk_pm`` set
``skip_trafos[t]=1`` if the ``t``-th serial transformation should be
computed, otherwise set ``skip_trafos[t]=0``, see
Section [sec:skip-trafo] for more details.

``in`` and ``out`` point to the complex valued input and output arrays
of the transform, which may be the same (yielding an in-place
transform). These arrays are overwritten during planning, unless
``PFFT_ESTIMATE | PFFT_NO_TUNE`` is used in the flags. (The arrays need
not be initialized, but they must be allocated.)

``sign`` is the sign of the exponent in the formula that defines the
Fourier transform. It can be -1 (= ``PFFT_FORWARD``) or +1 (=
``PFFT_BACKWARD``).

``pfft_flags`` is a bitwise OR (’``|``\ ’) of zero or more planner
flags, as defined in Section [sec:flags].

PFFT computes an unnormalized transform: computing a forward followed by
a backward transform (or vice versa) will result in the original data
multiplied by the size of the transform (the product of the dimensions
``n[t]``).

Real-to-Complex FFT
~~~~~~~~~~~~~~~~~~~

::

    pfft_plan pfft_plan_dft_r2c_3d(
        const ptrdiff_t *n, double *in, pfft_complex *out, MPI_Comm comm_cart,
        int sign, unsigned pfft_flags);
    pfft_plan pfft_plan_dft_r2c(
        int rnk_n, const ptrdiff_t *n, double *in, pfft_complex *out, MPI_Comm comm_cart,
        int sign, unsigned pfft_flags);
    pfft_plan pfft_plan_many_dft_r2c(
        int rnk_n, const ptrdiff_t *n, const ptrdiff_t *ni, const ptrdiff_t *no,
        ptrdiff_t howmany, const ptrdiff_t *iblock, const ptrdiff_t *oblock,
        double *in, pfft_complex *out, MPI_Comm comm_cart,
        int sign, unsigned pfft_flags);
    pfft_plan pfft_plan_many_dft_r2c_skipped(
        int rnk_n, const ptrdiff_t *n, const ptrdiff_t *ni, const ptrdiff_t *no,
        ptrdiff_t howmany, const ptrdiff_t *iblock, const ptrdiff_t *oblock,
        const int *skip_trafos, double *in, pfft_complex *out, MPI_Comm comm_cart,
        int sign, unsigned pfft_flags);

Plan a parallel, real-input/complex-output discrete Fourier transform
(DFT) in two or more dimensions, returning an ``pfft_plan``. The planner
returns NULL if the plan cannot be created.

Arguments:

``rnk_n``, ``n``, ``ni``, ``no``, ``howmany``, ``iblock``, ``oblock``,
``comm_cart`` must be the same as passed to the corresponding
``pfft_local_size_dft_r2c`` function, see Section [sec:local-size-r2c].

``in`` and ``out`` point to the real valued input and complex valued
output arrays of the transform, which may be the same (yielding an
in-place transform). These arrays are overwritten during planning,
unless ``PFFT_ESTIMATE | PFFT_NO_TUNE`` is used in the flags. (The
arrays need not be initialized, but they must be allocated.)

``sign`` is the sign of the exponent in the formula that defines the
Fourier transform. It can be -1 (= ``PFFT_FORWARD``) or +1 (=
``PFFT_BACKWARD``). Note that this parameter is not part of the FFTW-MPI
interface, where r2c transforms are defined to be forward transforms.
However, the backward transform can be easily realized by an additional
conjugation of the complex outputs as done by PFFT.

Complex-to-Real FFT
~~~~~~~~~~~~~~~~~~~

::

    pfft_plan pfft_plan_dft_c2r_3d(
        const ptrdiff_t *n, pfft_complex *in, double *out, MPI_Comm comm_cart,
        int sign, unsigned pfft_flags);
    pfft_plan pfft_plan_dft_c2r(
        int rnk_n, const ptrdiff_t *n, pfft_complex *in, double *out, MPI_Comm comm_cart,
        int sign, unsigned pfft_flags);
    pfft_plan pfft_plan_many_dft_c2r(
        int rnk_n, const ptrdiff_t *n, const ptrdiff_t *ni, const ptrdiff_t *no,
        ptrdiff_t howmany, const ptrdiff_t *iblock, const ptrdiff_t *oblock,
        pfft_complex *in, double *out, MPI_Comm comm_cart,
        int sign, unsigned pfft_flags);
    pfft_plan pfft_plan_many_dft_c2r_skipped(
        int rnk_n, const ptrdiff_t *n, const ptrdiff_t *ni, const ptrdiff_t *no,
        ptrdiff_t howmany, const ptrdiff_t *iblock, const ptrdiff_t *oblock,
        const int *skip_trafos, pfft_complex *in, double *out, MPI_Comm comm_cart,
        int sign, unsigned pfft_flags);

Plan a parallel, complex-input/real-output discrete Fourier transform
(DFT) in two or more dimensions, returning an ``pfft_plan``. The planner
returns NULL if the plan cannot be created.

Arguments:

``rnk_n``, ``n``, ``ni``, ``no``, ``howmany``, ``iblock``, ``oblock``,
``comm_cart`` must be the same as passed to the corresponding
``pfft_local_size_dft_c2r`` function, see Section [sec:local-size-c2r].

``in`` and ``out`` point to the complex valued input and real valued
output arrays of the transform, which may be the same (yielding an
in-place transform). These arrays are overwritten during planning,
unless ``PFFT_ESTIMATE | PFFT_NO_TUNE`` is used in the flags. (The
arrays need not be initialized, but they must be allocated.)

``sign`` is the sign of the exponent in the formula that defines the
Fourier transform. It can be -1 (= ``PFFT_FORWARD``) or +1 (=
``PFFT_BACKWARD``). Note that this parameter is not part of the FFTW-MPI
interface, where c2r transforms are defined to be backward transforms.
However, the forward transform can be easily realized by an additional
conjugation of the complex inputs as done by PFFT.

Real-to-Real FFT
~~~~~~~~~~~~~~~~

::

    pfft_plan pfft_plan_r2r_3d(
        const ptrdiff_t *n, double *in, double *out, MPI_Comm comm_cart,
        const pfft_r2r_kind *kinds, unsigned pfft_flags);
    pfft_plan pfft_plan_r2r(
        int rnk_n, const ptrdiff_t *n, double *in, double *out, MPI_Comm comm_cart,
        const pfft_r2r_kind *kinds, unsigned pfft_flags);
    pfft_plan pfft_plan_many_r2r(
        int rnk_n, const ptrdiff_t *n, const ptrdiff_t *ni, const ptrdiff_t *no,
        ptrdiff_t howmany, const ptrdiff_t *iblock, const ptrdiff_t *oblock,
        double *in, double *out, MPI_Comm comm_cart,
        const pfft_r2r_kind *kinds, unsigned pfft_flags);
    pfft_plan pfft_plan_many_r2r_skipped(
        int rnk_n, const ptrdiff_t *n, const ptrdiff_t *ni, const ptrdiff_t *no,
        ptrdiff_t howmany, const ptrdiff_t *iblock, const ptrdiff_t *oblock,
        const int *skip_trafos, double *in, double *out, MPI_Comm comm_cart,
        const pfft_r2r_kind *kinds, unsigned pfft_flags);

Plan a parallel, real input/output (r2r) transform in two or more
dimensions, returning an ``pfft_plan``. The planner returns NULL if the
plan cannot be created.

Arguments:

``rnk_n``, ``n``, ``ni``, ``no``, ``howmany``, ``iblock``, ``oblock``,
``comm_cart`` must be the same as passed to the corresponding
``pfft_local_size_r2r`` function, see Section [sec:local-size-r2r].

``in`` and ``out`` point to the real valued input and output arrays of
the transform, which may be the same (yielding an in-place transform).
These arrays are overwritten during planning, unless
``PFFT_ESTIMATE | PFFT_NO_TUNE`` is used in the flags. (The arrays need
not be initialized, but they must be allocated.)

The array ``kinds`` of length ``rnk_n`` specifies the kind of r2r
transform that is computed in the corresponding dimensions. Just like
FFTW-MPI we compute the separable product formed by taking each
transform kind along the corresponding dimension, one dimension after
another.

FFT Execution Timer
-------------------

PFFT offers an easy way to perform run time measurements and print/write
the results.

Basis Run Time Measurements
~~~~~~~~~~~~~~~~~~~~~~~~~~~

PFFT-plans automatically accumulate the local run times of every call to
``pfft_execute``. For most applications it is sufficient to print run
time of a plan ``ths`` averaged over all runs with

::

    void pfft_print_average_timer(
        const pfft_plan ths, MPI_Comm comm);

Note, that for each timer the maximum time over all processes is reduced
to rank ``0`` of communicator ``comm``, i.e., a call to ``MPI_Reduce``
is performed and the output is only printed on this process. The
following function works in the same way but prints more verbose output

::

    void pfft_print_average_timer_adv(
        const pfft_plan ths, MPI_Comm comm);

To write the averaged run time of plan ``ths`` into a file called
``name`` use

::

    void pfft_write_average_timer(
        const pfft_plan ths, const char *name, MPI_Comm comm);
    void pfft_write_average_timer_adv(
        const pfft_plan ths, const char *name, MPI_Comm comm);

Again, the output is only written on rank ``0`` of communicator
``comm``.

Discard all the recorded run times with

::

    void pfft_reset_timer(
        pfft_plan ths);

This function is called per default at the end of every PFFT plan
creation function.

Advanced Timer Manipulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to access the run times directly a new typedef ``pfft_timer``
is introduced. The following function returns a copy of the timer
corresponding to PFFT plan ``ths``

::

    pfft_timer pfft_get_timer(
        const pfft_plan ths);

Note that the memory of the returned ``pfft_timer`` must be released
with

::

    void pfft_destroy_timer(
        pfft_timer ths);

as soon as the timer is not needed anymore.

In the following we introduce some routines to perform basic operations
on timers. For all functions with a ``pfft_timer`` return value you must
use ``pfft_destroy_timer`` in order to release the allocated memory of
the timer. Create a copy of a PFFT-timer ``orig`` with

::

    pfft_timer pfft_copy_timer(
        const pfft_timer orig);

Compute the average, local time over all runs of a timer ``ths`` with

::

    void pfft_average_timer(
        pfft_timer ths);

Create a new timer that contains the sum of two timers ``sum1`` and
``sum2`` with

::

    pfft_timer pfft_add_timers(
        const pfft_timer sum1, const pfft_timer sum2);

Create a timer that contains the maximum times of all the timers ``ths``
from all processes belonging to communicator ``comm`` with

::

    pfft_timer pfft_reduce_max_timer(
        const pfft_timer ths, MPI_Comm comm);

Since this function calls ``MPI_Reduce``, only the first process (rank
0) of ``comm`` will get the desired data while all the other processes
have timers with undefined values.

Note, that you can not access the elements of a timer directly, since it
is only a pointer to a ``struct``. However, PFFT offers a routine that
creates an array and copies all the entries of the timer into it

::

    double* pfft_convert_timer2vec(
        const pfft_timer ths);

Remember to use ``free`` in order to release the allocated memory of the
returned array at the moment it is not needed anymore. The entries of
the returned array are ordered as follows:

dimension of the process mesh ``rnk_pm``

number of serial trafos ``rnk_trafo``

number of global remaps ``rnk_remap``

number of ``pfft_execute`` runs ``iter``

local run time of all runs

``rnk_n`` local times of the serial trafos

``rnk_remap`` local times of the global remaps

2 times of the global remaps that are only necessary for
three-dimensional FFTs on three-dimensional process meshes

time for computing twiddled input (as needed for ``PFFT_SHIFTED_OUT``)

time for computing twiddled output (as needed for ``PFFT_SHIFTED_IN``)

The complementary function

::

    pfft_timer pfft_convert_vec2timer(
        const double *times);

creates a timer and fills it’s entries with the data from array
``times``. Thereby, the entries of ``times`` must be in the same order
as above.

Ghost Cell Communication
------------------------

In the following we describe the PFFT ghost cell communication module.
At the moment, PFFT ghost cell communication is restricted to
three-dimensional arrays.

Assume a three-dimensional array ``data`` of size ``n[3]`` that is
distributed in blocks such that each process has a local copy of
``data[k[0],k[1],k[2]]`` with

::

    local_start[t] <= k[t] < local_start[t] + local_n[t]

Here and in the following, we assume ``t=0,1,2``. The “classical” ghost
cell exchange communicates all the necessary data between neighboring
processes, such that each process gets a local copy of
``data[k[0],k[1],k[2]]`` with

::

    local_gc_start[t] <= k[t] < local_gc_start[t] + local_ngc[t]

where

::

    local_gc_start[t] = local_start[t] - gc_below[t];
    local_ngc[t] = local_n[t] + gc_below[t] + gc_above[t];

I.e., the local array block is increased in every dimension by
``gc_below`` elements below and ``gc_above`` elements above. Hereby, the
``data`` is wrapped periodically whenever ``k[t]`` exceeds the array
dimensions. The number of ghost cells in every dimension can be chosen
independently and can be arbitrary large, i.e., PFFT ghost cell
communication also handles the case where the requested data exceeds
next neighbor communication. The number of ghost cells can even be
bigger than the array size, which results in multiple local copies of
the same data elements at every process. However, the arrays
``gc_below`` and ``gc_above`` must be equal among all MPI processes.

PFFT ghost cell communication can work on both, the input and output
array distributions. Substitute ``local_n`` and ``local_start`` by
``local_ni`` and ``local_i_start`` if you are interested in ghost cell
communication of the input array. For ghost cell communication of the
output array, substitute ``local_n`` and ``local_start`` by ``local_no``
and ``local_o_start``.

Using Ghost Cell Plans
~~~~~~~~~~~~~~~~~~~~~~

We introduce a new datatype ``pfft_gcplan`` that stores all the
necessary information for ghost cell communication. Using a ghost cell
plan follows the typical workflow: At first, determine the parallel data
distribution; cf. Section [sec:gc:local-size]. Next, create a ghost cell
plan; cf. Section [sec:gc:plan-cdata] and Section [sec:gc:plan-rdata].
Execute the ghost cell communication with one of the following two
collective functions

::

    void pfft_exchange(
        pfft_gcplan ths);
    void pfft_reduce(
        pfft_gcplan ths);

Hereby, a ghost cell exchange creates duplicates of local data elements
on next neighboring processes, while a ghost cell reduce is the adjoint
counter part of the exchange, i.e., it adds the sum of all the
duplicates of a local data element to the original data element.
Finally, free the allocated memory with

::

    void pfft_destroy_gcplan(
        pfft_gcplan ths);

if the plan is not needed anymore. Passing a freed plan to
``pfft_exchange`` or ``pfft_reduce`` results in undefined behavior.

Data Distribution
~~~~~~~~~~~~~~~~~

Corresponding to the three interface layers for FFT planning, there are
the following three layers for computing the ghost cell data
distribution:

::

    ptrdiff_t pfft_local_size_gc_3d(
        const ptrdiff_t *local_n, const ptrdiff_t *local_start,
        const ptrdiff_t *gc_below, const ptrdiff_t *gc_above,
        ptrdiff_t *local_ngc, ptrdiff_t *local_gc_start);
    ptrdiff_t pfft_local_size_gc(
        int rnk_n, 
        const ptrdiff_t *local_n, const ptrdiff_t *local_start,
        const ptrdiff_t *gc_below, const ptrdiff_t *gc_above,
        ptrdiff_t *local_ngc, ptrdiff_t *local_gc_start);
    ptrdiff_t pfft_local_size_many_gc(
        int rnk_n,
        const ptrdiff_t *local_n, const ptrdiff_t *local_start,
        ptrdiff_t howmany,
        const ptrdiff_t *gc_below, const ptrdiff_t *gc_above,
        ptrdiff_t *local_ngc, ptrdiff_t *local_gc_start);

Hereby, ``rnk_n`` and ``howmany`` must be the exactly same variables
that were used for the PFFT plan creation. However, only the case
``rnk_n==3`` is completely implemented at the moment. The local array
size ``local_n`` must be equal to ``local_ni`` or ``local_no`` (computed
by an appropriate call of ``pfft_local_size``; cf.
Section [sec:local-size]) depending on whether the ghost cell plan works
on the FFT input or output array. Analogously, ``local_start`` becomes
``local_i_start`` or ``local_o_start``. The number of ghost cells is
given by the two arrays ``gc_below`` and ``gc_above`` that must be equal
among all MPI processes. All the ghost cell data distribution functions
return the local array plus ghost cell size ``local_ngc`` and the
corresponding offset ``local_gc_start`` as two arrays of length
``rnk_n``. In addition, the ``ptrdiff_t`` return value gives the number
of data elements that are necessary in order to store the array plus
ghost cells.

Note, that the array distribution functions do not distinguish between
real and complex valued data. That is because ``local_n`` and
``local_start`` count array elements in units of complex or real
depending on the transform. In addition, it does not matter if the local
array is transposed or not, i.e., it is not necessary to pass the flags
``PFFT_TRANSPOSED_IN`` and ``PFFT_TRANSPOSED_OUT`` to the ghost cell
distribution function. In constrast, the ghost cell plan creation
depends on the transform type as well as the transposition flags.

Memory Allocation
~~~~~~~~~~~~~~~~~

In most applications we must ensure that the data array is large enough
to suit the memory requirements of a parallel FFT and the ghost cell
communication. The following two code snippets illustrate the correct
allocation of memory in for complex valued and real valued arrays.

::

    /* Get parameters of data distribution */
    /* alloc_local, local_no, local_o_start are given in complex units */
    /* local_ni, local_i_start are given in real units */
    alloc_local = pfft_local_size_dft_r2c_3d(n, comm_cart_2d, PFFT_TRANSPOSED_NONE,
        local_ni, local_i_start, local_no, local_o_start);

    /* alloc_local_gc, local_ngc, local_gc_start are given in complex units */
    alloc_local_gc = pfft_local_size_gc_3d(
        local_no, local_o_start, gc_below, gc_above,
        local_ngc, local_gc_start);

    /* Allocate enough memory for FFT and ghost cells */
    pfft_complex *cdata = pfft_alloc_complex(alloc_local_gc > alloc_local ? alloc_local_gc : alloc_local);

Here, ``alloc_local`` gives the number of data elements that are
necessary to hold all steps of the parallel FFT, while
``alloc_local_gc`` gives the number of data elements that are necessary
to hold all steps of the ghost cell communication. Note that we took the
maximum of these both numbers as argument for ``pfft_alloc_complex``.
The code snippet for real valued arrays looks very similar.

::

    /* Get parameters of data distribution */
    /* alloc_local, local_no, local_o_start are given in complex units */
    /* local_ni, local_i_start are given in real units */
    alloc_local = pfft_local_size_dft_r2c_3d(n, comm_cart_2d, PFFT_TRANSPOSED_NONE,
        local_ni, local_i_start, local_no, local_o_start);

    /* alloc_local_gc, local_ngc, local_gc_start are given in real units */
    alloc_local_gc = pfft_local_size_gc_3d(
        local_ni, local_i_start, gc_below, gc_above,
        local_ngc, local_gc_start);

    /* Allocate enough memory for FFT and ghost cells */
    double *rdata = pfft_alloc_real(alloc_local_gc > 2*alloc_local ? alloc_local_gc : 2*alloc_local);

Note that the number of real valued data elements is given by two times
``alloc_local`` for r2c transforms, whereas the last line would change
into

::

    double *rdata = pfft_alloc_real(alloc_local_gc > alloc_local ? alloc_local_gc : alloc_local);

for r2r transforms.

Plan Creation for Complex Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following functions create ghost cell plans that operate on complex
valued arrays, i.e.,

c2c inputs,

c2c outputs,

r2c outputs (use flag ``PFFT_GC_C2R``), and

c2r inputs (use flag ``PFFT_GC_R2C``).

Corresponding to the three interface layers for FFT planning, there are
the following three layers for creating a complex valued ghost cell
plan:

::

    pfft_gcplan pfft_plan_cgc_3d(
        const ptrdiff_t *n,
        const ptrdiff_t *gc_below, const ptrdiff_t *gc_above,
        pfft_complex *data, MPI_Comm comm_cart, unsigned gc_flags);
    pfft_gcplan pfft_plan_cgc(
        int rnk_n, const ptrdiff_t *n,
        const ptrdiff_t *gc_below, const ptrdiff_t *gc_above,
        pfft_complex *data, MPI_Comm comm_cart, unsigned gc_flags);
    pfft_gcplan pfft_plan_many_cgc(
        int rnk_n, const ptrdiff_t *n,
        ptrdiff_t howmany, const ptrdiff_t *block,
        const ptrdiff_t *gc_below, const ptrdiff_t *gc_above,
        pfft_complex *data, MPI_Comm comm_cart, unsigned gc_flags);

Hereby, ``rnk_n``, ``n``, ``howmany`` and ``comm_cart`` must be the
variables that were used for the PFFT plan creation. However, only the
case ``rnk_n==3`` is completely implemented at the moment. Remember that
``n`` is the logical FFT size just as it is the case for FFT planning.
The block size ``block`` must be equal to ``iblock`` or ``oblock``
depending on whether the ghost cell plan works on the FFT input or
output array. Analogously, ``data`` becomes ``in`` or ``out``. Set the
number of ghost cells by ``gc_below`` and ``gc_above`` as described in
Section [sec:gc]. The flags ``gc_flags`` must be set appropriately to
the flags that were passed to the FFT planner. Table [tab:map-cgcflags]
shows the ghost cell planner flags that must be set in dependence on the
listed FFT planner flags.

[h]

+----------------------------+-------------------------------+
| FFT flag                   | ghost cell flag               |
+============================+===============================+
| ``PFFT_TRANSPOSED_NONE``   | ``PFFT_GC_TRANSPOSED_NONE``   |
+----------------------------+-------------------------------+
| ``PFFT_TRANSPOSED_IN``     | ``PFFT_GC_TRANSPOSED``        |
+----------------------------+-------------------------------+
| ``PFFT_TRANSPOSED_OUT``    | ``PFFT_GC_TRANSPOSED``        |
+----------------------------+-------------------------------+

[tab:map-cgcflags]

In addition, we introduce the flag ``PFFT_GC_R2C`` (and its equivalent
``PFFT_GC_C2R``) to handle the complex array storage format of r2c and
c2r transforms. In fact, these two flags imply an ordinary complex
valued ghost cell communication on an array of size
``n[0] x ... x n[rnk_n-2] x (n[rnk_n-1]/2``\ 1)+. Please note that we
wrongly assume periodic boundary conditions in this case. Therefore, you
should ignore the data elements with the last index behind
``n[rnk_n-1]/2``.

Plan Creation for Real Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following functions create ghost cell plans that operate on real
valued arrays, i.e.,

r2r inputs,

r2r outputs,

r2c inputs, and

c2r outputs.

Corresponding to the three interface layers for FFT planning, there are
the following three layers for creating a real valued ghost cell plan:

::

    pfft_gcplan pfft_plan_rgc_3d(
        const ptrdiff_t *n,
        const ptrdiff_t *gc_below, const ptrdiff_t *gc_above,
        double *data, MPI_Comm comm_cart, unsigned gc_flags);
    pfft_gcplan pfft_plan_rgc(
        int rnk_n, const ptrdiff_t *n,
        const ptrdiff_t *gc_below, const ptrdiff_t *gc_above,
        double *data, MPI_Comm comm_cart, unsigned gc_flags);
    pfft_gcplan pfft_plan_many_rgc(
        int rnk_n, const ptrdiff_t *n,
        ptrdiff_t howmany, const ptrdiff_t *block,
        const ptrdiff_t *gc_below, const ptrdiff_t *gc_above,
        double *data, MPI_Comm comm_cart, unsigned gc_flags);

Hereby, ``rnk_n``, ``n``, ``howmany`` and ``comm_cart`` must be the
variables that were used for the PFFT plan creation. Remember that ``n``
is the logical FFT size just as it is the case for FFT planning. The
block size ``block`` must be equal to ``iblock`` or ``oblock`` depending
on whether the ghost cell plan works on the FFT input or output array.
Analogously, ``data`` becomes ``in`` or ``out``. Set the number of ghost
cells by ``gc_below`` and ``gc_above`` as described in
Section [sec:gc:local-size]. The flags ``gc_flags`` must be set
appropriately to the flags that were passed to the FFT planner.
Table [tab:map-rgcflags] shows the ghost cell planner flags that must be
set in dependence on the listed FFT planner flags.

[h]

+----------------------------+-------------------------------+
| FFT flag                   | ghost cell flag               |
+============================+===============================+
| ``PFFT_TRANSPOSED_NONE``   | ``PFFT_GC_TRANSPOSED_NONE``   |
+----------------------------+-------------------------------+
| ``PFFT_TRANSPOSED_IN``     | ``PFFT_GC_TRANSPOSED``        |
+----------------------------+-------------------------------+
| ``PFFT_TRANSPOSED_OUT``    | ``PFFT_GC_TRANSPOSED``        |
+----------------------------+-------------------------------+
| ``PFFT_PADDED_R2C``        | ``PFFT_GC_PADDED_R2C``        |
+----------------------------+-------------------------------+
| ``PFFT_PADDED_C2R``        | ``PFFT_GC_PADDED_C2R``        |
+----------------------------+-------------------------------+

[tab:map-rgcflags]

Note that the flag ``PFFT_GC_PADDED_R2C`` (or its equivalent
``PFFT_GC_PADDED_C2R``) implies an ordinary real valued ghost cell
communication on an array of size
``n[0] x ... x n[rnk_n-2] x 2*(n[rnk_n-1]/2``\ 1)+. Especially, the
padding elements will be handles as normal data points, i.e., you must
we aware that the numbers of ghost cells ``gc_below[rnk_n-1]`` and
``gc_above[rnk_n-1]`` include the number of padding elements.

Inofficial Flags
~~~~~~~~~~~~~~~~

Ghost Cell Execution Timer
~~~~~~~~~~~~~~~~~~~~~~~~~~

PFFT ghost cell plans automatically accumulate the local run times of
every call to ``pfft_exchange`` and ``pfft_reduce``. For most
applications it is sufficient to print run time of a plan ``ths``
averaged over all runs with

::

    void pfft_print_average_gctimer(
        const pfft_gcplan ths, MPI_Comm comm);

Note, that for each timer the maximum time over all processes is reduced
to rank ``0`` of communicator ``comm``, i.e., a call to ``MPI_Reduce``
is performed and the output is only printed on this process. The
following function works in the same way but prints more verbose output

::

    void pfft_print_average_gctimer_adv(
        const pfft_gcplan ths, MPI_Comm comm);

To write the averaged run time of a ghost cell plan ``ths`` into a file
called ``name`` use

::

    void pfft_write_average_gctimer(
        const pfft_gcplan ths, const char *name, MPI_Comm comm);
    void pfft_write_average_gctimer_adv(
        const pfft_gcplan ths, const char *name, MPI_Comm comm);

Again, the output is only written on rank ``0`` of communicator
``comm``.

Discard all the recorded run times with

::

    void pfft_reset_gctimers(
        pfft_gcplan ths);

This function is called per default at the end of every ghost cell plan
creation function.

In order to access the run times directly a new typedef ``pfft_timer``
is introduced. The following functions return a copy of the timer
corresponding to ghost cell plan ``ths`` that accumulated the time for
ghost cell exchange or ghost cell reduce, respectively:

::

    pfft_gctimer pfft_get_gctimer_exg(
        const pfft_gcplan ths);
    pfft_gctimer pfft_get_gctimer_red(
        const pfft_gcplan ths);

Note that the memory of the returned ``pfft_gctimer`` must be released
with

::

    void pfft_destroy_gctimer(
        pfft_gctimer ths);

as soon as the timer is not needed anymore.

In the following we introduce some routines to perform basic operations
on timers. For all functions with a ``pfft_gctimer`` return value you
must use ``pfft_destroy_gctimer`` in order to release the allocated
memory of the timer. Create a copy of a ghost cell timer ``orig`` with

::

    pfft_gctimer pfft_copy_gctimer(
        const pfft_gctimer orig);

Compute the average, local time over all runs of a timer ``ths`` with

::

    void pfft_average_gctimer(
        pfft_gctimer ths);

Create a new timer that contains the sum of two timers ``sum1`` and
``sum2`` with

::

    pfft_gctimer pfft_add_gctimers(
        const pfft_gctimer sum1, const pfft_gctimer sum2);

Create a timer that contains the maximum times of all the timers ``ths``
from all processes belonging to communicator ``comm`` with

::

    pfft_gctimer pfft_reduce_max_gctimer(
        const pfft_gctimer ths, MPI_Comm comm);

Since this function calls ``MPI_Reduce``, only the first process (rank
0) of ``comm`` will get the desired data while all the other processes
have timers with undefined values.

Note, that you can not access the elements of a timer directly, since it
is only a pointer to a ``struct``. However, PFFT offers a routine that
creates an array and copies all the entries of the timer into it

::

    void pfft_convert_gctimer2vec(
        const pfft_gctimer ths, double *times);

Remember to use ``free`` in order to release the allocated memory of the
returned array at the moment it is not needed anymore. The entries of
the returned array are ordered as follows:

number of ``pfft_execute`` runs ``iter``

local run time of all runs

local run time of zero padding (make room for incoming ghost cells and
init with zeros)

local run time of the ghost cell exchange or reduce (depending on the
timer)

The complementary function

::

    pfft_gctimer pfft_convert_vec2gctimer(
        const double *times);

creates a timer and fills it’s entries with the data from array
``times``. Thereby, the entries of ``times`` must be in the same order
as above.

Useful Tools
------------

The following functions are useful tools but are not necessarily needed
to perform parallel FFTs.

Initializing Complex Inputs and Checking Outputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To fill a complex array ``data`` with reproducible, complex values you
can use one of the functions

::

    void pfft_init_input_complex_3d(
        const ptrdiff_t *n,
        const ptrdiff_t *local_n, const ptrdiff_t *local_n_start,
        pfft_complex *data);
    void pfft_init_input_complex(
        int rnk_n, const ptrdiff_t *n,
        const ptrdiff_t *local_n, const ptrdiff_t *local_start,
        pfft_complex *data);

Hereby, the arrays ``n``, ``local_n`` and ``local_n_start`` of length
``rnk_n`` (``rnk_n==3`` for ``_3d``) give the size of the FFT, the local
array size and the local array offset as computed by the array
distribution functions described in Section [sec:local-size] The
functions

::

    double pfft_check_output_complex_3d(
        const ptrdiff_t *n, 
        const ptrdiff_t *local_n, const ptrdiff_t *local_n_start,
        const pfft_complex *data, MPI_Comm comm);
    double pfft_check_output_complex(
        int rnk_n, const ptrdiff_t *n,
        const ptrdiff_t *local_n, const ptrdiff_t *local_start,
        const pfft_complex *data, MPI_Comm comm);

compute the :math:`l_1`-norm between the elements of array ``data`` and
values produced by ``pfft_init_input_complex_3d``,
``pfft_init_input_complex``. In addition, we supply the following
functions for setting all the input data to zero at once

::

    void pfft_clear_input_complex_3d(
        const ptrdiff_t *n,
        const ptrdiff_t *local_n, const ptrdiff_t *local_n_start,
        pfft_complex *data);
    void pfft_clear_input_complex(
        int rnk_n, const ptrdiff_t *n,
        const ptrdiff_t *local_n, const ptrdiff_t *local_start,
        pfft_complex *data);

Note, that these functions can be combined for a quick consistency check
of the FFT. Since a forward FFT followed by a backward FFT reproduces
the inputs up to a scaling factor, the following code snippet should
give a result equal to zero up to machine precision.

::

    /* Initialize input with random numbers */
    pfft_init_input_complex_3d(n, local_ni, local_i_start,
        in);

    /* execute parallel forward FFT */
    pfft_execute(plan_forw);

    /* clear the old input */
    if(in != out) 
      pfft_clear_input_complex_3d(n, local_ni, local_i_start, in);

    /* execute parallel backward FFT */
    pfft_execute(plan_back);

    /* Scale data */
    for(ptrdiff_t l=0; l < local_ni[0] * local_ni[1] * local_ni[2]; l++)
      in[l] /= (n[0]*n[1]*n[2]);

    /* Print error of back transformed data */
    err = pfft_check_output_complex_3d(n, local_ni, local_i_start, in, comm_cart_2d);
    pfft_printf(comm_cart_2d, "Error after one forward and backward trafo of size n=(%td, %td, %td):\n", n[0], n[1], n[2]);
    pfft_printf(comm_cart_2d, "maxerror = %6.2e;\n", err);

Hereby, we set all inputs equal to zero after the forward FFT in order
to be sure that all the final results are actually computed by the
backward FFT instead of being a buggy relict of the forward transform.

Initializing Real Inputs and Checking Outputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To fill a real array ``data`` with reproducible, real values use one of
the functions

::

    void pfft_init_input_real_3d(
        const ptrdiff_t *n,
        const ptrdiff_t *local_n, const ptrdiff_t *local_n_start,
        double *data);
    void pfft_init_input_real(
        int rnk_n, const ptrdiff_t *n,
        const ptrdiff_t *local_n, const ptrdiff_t *local_start,
        double *data);

Hereby, the arrays ``n``, ``local_n`` and ``local_n_start`` give the
size of the FFT, the local array size and the local array offset as
computed by the array distribution functions described in
Section [sec:local-size] The functions

::

    double pfft_check_output_real_3d(
        const ptrdiff_t *n,
        const ptrdiff_t *local_n, const ptrdiff_t *local_n_start,
        const pfft_complex *data, MPI_Comm comm);
    double pfft_check_output_real(
        int rnk_n, const ptrdiff_t *n,
        const ptrdiff_t *local_n, const ptrdiff_t *local_start,
        const pfft_complex *data, MPI_Comm comm);

compute the :math:`l_1`-norm between the elements of array ``data`` and
values produced by ``pfft_init_input_real_3d``,
``pfft_init_input_real``. In addition, we supply the following functions
for setting all the input data to zero at once

::

    void pfft_clear_input_real_3d(
        const ptrdiff_t *n,
        const ptrdiff_t *local_n, const ptrdiff_t *local_n_start,
        double *data);
    void pfft_clear_input_real(
        int rnk_n, const ptrdiff_t *n,
        const ptrdiff_t *local_n, const ptrdiff_t *local_start,
        double *data);

Note, that both ``pfft_init_input_real*`` functions will set all array
elements to zero were ``local_n `` local\ :sub:`ns`\ tart >= n+. In
addition, both ``pfft_check_output_real*`` function will ignore all the
errors resulting from these elements. Therefore, it is safe to use all
these functions for a consistency check of a r2c transform followed by a
c2r transform since all padding elements will be ignored.

Initializing r2c/c2r Inputs and Checking Outputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The real inputs of a r2c transform can be initialized with the functions
decribed in Section [sec:init-data-3d-r2r]. However, generating suitable
inputs for a c2r transform requires more caution. In order to get real
valued results of a DFT the complex input coefficients need to satisfy
an radial Hermitian symmetry, i.e.,
:math:`X[{\ensuremath{\boldsymbol{k}}}] = {X^*[-{\ensuremath{\boldsymbol{k}}}]}`.
We use the following trick to generate the complex input values for c2r
transforms. Assume any :math:`{\ensuremath{\boldsymbol{N}}}`-periodic
complex valued function :math:`f`. It can be easily shown that the
values
:math:`X[{\ensuremath{\boldsymbol{k}}}] := \frac{1}{2}\left(f({\ensuremath{\boldsymbol{k}}})+f^*(-{\ensuremath{\boldsymbol{k}}})\right)`
satisfy the radial Hermitian symmetry.

To fill a complex array ``data`` with reproducible, complex values that
fulfill the radial Hermitian symmetry use one of the functions

::

    void pfft_init_input_complex_hermitian_3d(
        const ptrdiff_t *n,
        const ptrdiff_t *local_n, const ptrdiff_t *local_n_start,
        double *data);
    void pfft_init_input_complex_hermitian(
        int rnk_n, const ptrdiff_t *n,
        const ptrdiff_t *local_n, const ptrdiff_t *local_start,
        double *data);

Hereby, the arrays ``n``, ``local_n`` and ``local_n_start`` give the
size of the FFT, the local array size and the local array offset as
computed by the array distribution functions described in
Section [sec:local-size] The functions

::

    double pfft_check_output_complex_hermitian_3d(
        const ptrdiff_t *n,
        const ptrdiff_t *local_n, const ptrdiff_t *local_n_start,
        const pfft_complex *data, MPI_Comm comm);
    double pfft_check_output_complex_hermitian(
        int rnk_n, const ptrdiff_t *n,
        const ptrdiff_t *local_n, const ptrdiff_t *local_start,
        const pfft_complex *data, MPI_Comm comm);

compute the :math:`l_1`-norm between the elements of array ``data`` and
values produced by ``pfft_init_input_complex_hermitian_3d``,
``pfft_init_input_complex_hermitian``. In addition, we supply the
following functions for setting all the input data to zero at once

::

    void pfft_clear_input_complex_hermitian_3d(
        const ptrdiff_t *n,
        const ptrdiff_t *local_n, const ptrdiff_t *local_n_start,
        pfft_complex *data);
    void pfft_clear_input_complex_hermitian(
        int rnk_n, const ptrdiff_t *n,
        const ptrdiff_t *local_n, const ptrdiff_t *local_start,
        pfft_complex *data);

Note, that these functions can also be used in order to generate complex
inputs with radial Hermitian symmetry for ordinary c2c transforms. Of
course the results of such a c2c DFT will have all imaginary parts equal
to zero up to machine precision.

Operations on Arrays of Type ``ptrdiff_t``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following routines are shortcuts for the elementwise manipulation of
``ptrdiff_t`` valued arrays. In the following, all arrays ``vec``,
``vec1``, and ``vec2`` are of length ``d`` and type ``ptrdiff_t``.

::

    ptrdiff_t pfft_prod_INT(
        int d, const ptrdiff_t *vec);

Returns the product over all elements of ``vec``.

::

    ptrdiff_t pfft_sum_INT(
        int d, const ptrdiff_t *vec);

Returns the sum over all elements of ``vec``.

::

    int pfft_equal_INT(
        int d, const ptrdiff_t *vec1, const ptrdiff_t *vec2);

Returns 1 if both arrays have equal entries, 0 otherwise.

::

    void pfft_vcopy_INT(
        int d, const ptrdiff_t *vec1,
        ptrdiff_t *vec2);

Copies the elements of ``vec1`` into ``vec2``.

::

    void pfft_vadd_INT(
        int d, const ptrdiff_t *vec1, const ptrdiff_t *vec2,
        ptrdiff_t *sum);

Fills ``sum`` with the componentwise sum of ``vec1`` and ``vec2``.

::

    void pfft_vsub_INT(
        int d, const ptrdiff_t *vec1, const ptrdiff_t *vec2,
        ptrdiff_t *sum);

Fills ``sum`` with the componentwise difference of ``vec1`` and
``vec2``.

Print Three-Dimensional Arrays in Parallel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the following routine to print the elements of a block decomposed
three-dimensional (real or complex valued) array ``data`` in a nicely
formatted way.

::

    void pfft_apr_real_3d(
        const double *data,
        const ptrdiff_t *local_n, const ptrdiff_t *local_start,
        const char *name, MPI_Comm comm);
    void pfft_apr_complex_3d(
        const pfft_complex *data,
        const ptrdiff_t *local_n, const ptrdiff_t *local_start,
        const char *name, MPI_Comm comm);

Obviously, this makes only sense for arrays of moderate size. The block
decomposition is given by ``local_n``, ``local_n_start`` as returned by
the array distribution function decribed in Section [sec:local-size].
Furthermore, some arbitrary string ``name`` can be added at the
beginning of each output - typically this will be the name of the array.
Communicator ``comm`` must be suitable to the block decomposition and is
used to synchronize the outputs over all processes.

Generalizations for the case where the dimensions of the local arrays
are permuted are given by

::

    void pfft_apr_real_permuted_3d(
        const double *data,
        const ptrdiff_t *local_n, const ptrdiff_t *local_start,
        int perm0, int perm1, int perm2,
        const char *name, MPI_Comm comm);
    void pfft_apr_complex_permuted_3d(
        const pfft_complex *data,
        const ptrdiff_t *local_n, const ptrdiff_t *local_start,
        int perm0, int perm1, int perm2,
        const char *name, MPI_Comm comm);

Hereby, ``perm0``, ``perm1``, and ``perm2`` give the array’s permutation
of dimension.

Reading Command Line Arguments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following function offers a simple way to read command line
arguments into an array ``parameter``.

::

    void pfft_get_args(
        int argc, char **argv, const char *name,
        int neededArgs, unsigned type,
        void *parameter);

Hereby, ``argc`` and ``argv`` are the standard argument of the ``main``
routine. Furthermore, ``name``, ``neededAgrs``, and ``type`` give the
name, number of entries and the type of the command line argument.
Supported types are ``PFFT_INT``, ``PFFT_PTRDIFF_T``, ``PFFT_FLOAT``,
``PFFT_DOUBLE``, and ``PFFT_UNSIGNED``, which denote the standard C type
that is used for typecasting. In addition, you can use the special type
``PFFT_SWITCH`` that is an integer type equal to one if the
corresponding command line argument is given. The array ``parameter``
must be of sufficient size to hold ``neededArgs`` elements of the given
data type. Special attention is given

For example, a program containing the following code snippet

::

    double x=0.1;
    pfft_get_args(argc, argv, "-pfft_x", 1, PFFT_DOUBLE, &x);
    int np[2]={2,1};
    pfft_get_args(argc, argv, "-pfft_np", 2, PFFT_INT, np);
    ptrdiff_t n[3]={32,32,32};
    pfft_get_args(argc, argv, "-pfft_n", 3, PFFT_PTRDIFF_T, n);
    int switch=0;
    pfft_get_args(argc, argv, "-pfft_on", 0, PFFT_SWITCH, switch);

that is executed via

::

    ./test -pfft_x 3.1 -pfft_np 2 3 -pfft_n 8 16 32 -pfft_on

will read ``x=3.1``, ``np[2] = (2,3)``, ``n[3]= (8,16,32)``, and turn on
the ``switch=1``. Note the address operator ``&`` in front of ``x`` in
the second line! Furthermore, note that the initialization of all
variables with default values before the call of ``pfft_get_args``
avoids trouble if the user does not provide all the command line
arguments.

Parallel Substitutes for ``vprintf``, ``fprintf``, and ``printf``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following functions are similar to the standard C function
``vfprintf``, ``fprintf`` and ``printf`` with the exception, that only
rank ``0`` within the given communicator ``comm`` will produce output.
The intension is to avoid the flood of messages that is produced when
simple ``printf`` statement are run in parallel.

::

    void pfft_vfprintf(
        MPI_Comm comm, FILE *stream, const char *format, va_list ap);
    void pfft_fprintf(
        MPI_Comm comm, FILE *stream, const char *format, ...);
    void pfft_printf(
        MPI_Comm comm, const char *format, ...);

Generating Periodic Cartesian Communicators
-------------------------------------------

Based on the processes that are part of the given communicator ``comm``
the following routine

::

    int pfft_create_procmesh_1d(
        MPI_Comm comm, int np0,
        MPI_Comm *comm_cart_1d);

allocates and creates a one-dimensional, periodic, Cartesian
communicator ``comm_cart_1d`` of size ``np0``. Thereby, a non-zero error
code is returned whenever ``np0`` does not fit the size of ``comm``. The
memory of the generated communicator should be released with
``MPI_Comm_free`` after usage. Analogously, use

::

    int pfft_create_procmesh_2d(
        MPI_Comm comm, int np0, int np1,
        MPI_Comm *comm_cart_2d);

in order to allocate and create two-dimensional, periodic, Cartesian
communicator ``comm_cart_2d`` of size ``np0*np1`` or

::

    int pfft_create_procmesh(
        int rnk_np, MPI_Comm comm, const int *np,
        MPI_Comm *comm_cart);

in order to allocate and create a ``rnk_np``-dimensional, periodic,
Cartesian communicator of size ``np[0]*np[1]*...*np[rnk_np-1]``. Hereby,
``np`` is an array of length ``rnk_np``. Again, the memory of the
generated communicator should be released with ``MPI_Comm_free`` after
usage.
