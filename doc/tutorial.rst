Tutorial
========

The following chapter describes the usage of the PFFT library at the
example of a simple test file in the first section, followed by the more
advanced features of PFFT in the next sections.

A first parallel transform - Three-dimensional FFT with two-dimensional data decomposition
------------------------------------------------------------------------------------------

We explain the basic steps for computing a parallel FFT with the PFFT
library at the example of the short test program given by
Listing [lst:man\ :sub:`c`\ 2c]. This test computes a three-dimensional
c2c-FFT on a two-dimensional process mesh. The source code can be found
in directory of the library’s source code tree.

After initializing MPI with and before calling any other PFFT routine
initialize the parallel FFT computations via

::

    void pfft_init(void);

MPI introduces the concept of communicators to store all the topological
information of the physical process layout. PFFT requires to be called
on a process mesh that corresponds to a periodic, Cartesian
communicator. We assist the user in creating such a communicator with
the following routine

::

    int pfft_create_procmesh_2d(
        MPI_Comm comm, int np0, int np1,
        MPI_Comm *comm_cart_2d);

This routine uses the processes within the communicator to create a
two-dimensional process grid of size x and stores it into the Cartesian
communicator . Note that is allocated by the routine and must be freed
with after usage. The input parameter is a communicator, indicating
which processes will participate in the transform. Choosing as implies
that the FFT is computed on all available processes.

At the next step we need to know the data decomposition of the input and
output array, that depends on the array sizes, the process grid and the
chosen parallel algorithm. Therefore, we call

::

    ptrdiff_t pfft_local_size_3d(
        ptrdiff_t *n, MPI_Comm comm_cart_2d, unsigned pfft_flags,
        ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
        ptrdiff_t *local_no, ptrdiff_t *local_o_start);

Hereby, , , , , are arrays of length :math:`3` that must be allocated.
The return value of this function equals the size of the local complex
array that needs to be allocated by every process. In most cases, this
coincides with the product of the local array sizes – but may be bigger,
whenever the parallel algorithm needs some extra storage. The input
value gives the three-dimensional FFT size and the flag serves to adjust
some details of the parallel execution. For the sake of simplicity, we
restrict our self to the case for a while and explain the more
sophisticated flags at a later point. The output arrays and give the
size and the offset of the local input array that result from the
parallel block distribution of the global input array, i.e., every
process owns the input data with for . Analogously, the output
parameters and contain the size and the offset of the local output
array.

Afterward, the input and output arrays must be allocated. Hereby,

::

    pfft_complex* pfft_alloc_complex(size_t size);

is a simple wrapper of , which in turn allocates the memory via to
ensure proper alignment for SIMD. Have a look at the FFTW user manual 
for more details on SIMD memory alignment and . Nevertheless, you can
also use any other dynamic memory allocation.

The planning of a single three-dimensional parallel FFT of size x x is
done by the function

::

    pfft_plan pfft_plan_dft_3d(
        ptrdiff_t *n, pfft_complex *in, pfft_complex *out,
        MPI_Comm comm_cart_2d, int sign, unsigned pfft_flags);

We provide the address of the input and output array by the pointers and
, respectively. An inplace transform is assumed if these pointers are
equal. The integer gives the sign in the exponential of the FFT.
Possible values are (:math:`-1`) and (:math:`+1`). Flags passed to the
planner via must coincide with the flags that were passed to . Otherwise
the data layout of the parallel execution may not match calculated local
array sizes. As return value we get a PFFT plan, some structure that
stores all the information needed to perform a parallel FFT.

Once the plan is generated, we are allowed to fill the input array .
Note, that per default the planning step will overwrite input array .
Therefore, you should not write any sensitive data into until the plan
was generated. For simplicity, our test program makes use of the library
function

::

    void pfft_init_input_complex_3d(
        ptrdiff_t *n, ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
        pfft_complex *in);

to fill the input array with some numbers. Alternatively, one can fill
the array with a function of choice and the following loop that takes
account of the parallel data layout

::

    ptrdiff_t m=0;
    for(ptrdiff_t k0=0; k0 < local_ni[0]; k0++)
      for(ptrdiff_t k1=0; k1 < local_ni[1]; k1++)
        for(ptrdiff_t k2=0; k2 < local_ni[2]; k2++)
          in[m++] = func(k0 + local_i_start[0],
                         k1 + local_i_start[1],
                         k2 + local_i_start[2]);

The parallel FFT is computed when we execute the generated plan via

::

    void pfft_execute(const pfft_plan plan);

Now, the results can be read from with an analogous three-dimensional
loop. If we do not want to execute another parallel FFT of the same
type, we free the allocated memory of the plan with

::

    void pfft_destroy_plan(pfft_plan plan);

Additionally, we use

::

    int MPI_Comm_free(MPI_Comm *comm);  

to free the communicator allocated by and

::

    void pfft_free(void *ptr);

to free memory allocated by . Finally, we exit MPI via

::

    int MPI_Finalize(void);

Porting FFTW-MPI based code to PFFT
-----------------------------------

We illustrate the close connection between FFTW-MPI and PFFT at a
three-dimensional MPI example analogous to the example given in the FFTW
manual .

Exactly the same task can be performed with PFFT as given in
Listing [lst:pfft\ :sub:`3`\ don1d].

::

    #include <pfft.h>
         
    int main(int argc, char **argv)
    {
        const ptrdiff_t n[3] = {..., ..., ...};
        pfft_plan plan;
        pfft_complex *data;
        ptrdiff_t alloc_local, local_ni[3], local_i_start[3], local_no[3], local_o_start[3], i, j, k;
        unsigned pfft_flags = 0;

        MPI_Init(&argc, &argv);
        pfft_init();

        /* get local data size and allocate */
        alloc_local = pfft_local_size_dft_3d(n, MPI_COMM_WORLD, pfft_flags,
                             local_ni, local_i_start,
                             local_no, local_o_start);
        data = pfft_alloc_complex(alloc_local);

        /* create plan for in-place forward DFT */
        plan = pfft_plan_dft_3d(n, data, data, MPI_COMM_WORLD,
                    PFFT_FORWARD, PFFT_ESTIMATE);

        /* initialize data to some function my_function(x,y,z) */
        for (i = 0; i < local_n[0]; ++i) 
          for (j = 0; j < n[1]; ++j) 
            for (k = 0; k < n[2]; ++k)
              data[i*n[1]*n[2] + j*n[2] + k] = my_function(local_i_start[0] + i, j, k);

        /* compute transforms, in-place, as many times as desired */
        pfft_execute(plan);

        pfft_destroy_plan(plan);

        MPI_Finalize();
    }

substitute by

substitute all prefixes and by

substitute all prefixes by

the integers , , become arrays of length 3

in

has additional input and additional outputs ,

The loop that inits becomes splitted along all three dimensions. We
could also use

First, All prefixes are substituted by

Now, the changes in order to use a two-dimensional process mesh are
marginal as can be seen in Listing [lst:pfft\ :sub:`3`\ don2d].

::

    #include <pfft.h>
         
    int main(int argc, char **argv)
    {
        const ptrdiff_t n[3] = {..., ..., ...};
        (red@const int np0 = ..., np1 = ...;@*)
        pfft_plan plan;
        pfft_complex *data;
        ptrdiff_t alloc_local, local_ni[3], local_i_start[3], local_no[3], local_o_start[3], i, j, k;
        unsigned pfft_flags = 0;
        (red@MPI_Comm comm_cart_2d;@*)

        MPI_Init(&argc, &argv);
        pfft_init();

        (red@/* create two-dimensional process grid of size np0 x np1 */@*)
        (red@pfft_create_procmesh_2d(MPI_COMM_WORLD, np0, np1,@*)
            (red@&comm_cart_2d);@*)
        
        /* get local data size and allocate */
        alloc_local = pfft_local_size_dft_3d(n, (red@comm_cart_2d@*), pfft_flags,
                             local_ni, local_i_start,
                             local_no, local_o_start);
        data = pfft_alloc_complex(alloc_local);

        /* create plan for in-place forward DFT */
        plan = pfft_plan_dft_3d(n, data, data, MPI_COMM_WORLD,
                    PFFT_FORWARD, PFFT_ESTIMATE);

        /* initialize data to some function my_function(x,y,z) */
        for (i = 0; i < local_n[0]; ++i) 
          for (j = 0; j < (red@local_n[1]@*); ++j) 
            for (k = 0; k < (red@local_n[2]@*); ++k)
              data[i*(red@local_n[1]*local_n[2]@*) + j*(red@local_n[2]@*) + k] =
                  my_function(local_i_start[0] + i,
                      (red@local_i_start[1] +@*) j,
                      (red@local_i_start[2] +@*) k);

        /* compute transforms, in-place, as many times as desired */
        pfft_execute(plan);

        pfft_destroy_plan(plan);

        MPI_Finalize();
    }

Errorcode for communicator creation
-----------------------------------

As we have seen the function

::

    int pfft_create_procmesh_2d(
        MPI_Comm comm, int np0, int np1,
        MPI_Comm *comm_cart_2d);

creates a two-dimensional, periodic, Cartesian communicator. The return
value (not used in Listing [lst:man\ :sub:`c`\ 2c]) is the forwarded
error code of . It is equal to zero if the communicator was created
successfully. The most common error is that the number of processes
within the input communicator does not fit . In this case the Cartesian
communicator is not generated and the return value is unequal to zero.
Therefore, a typical sanity check might look like

::

    /* Create two-dimensional process grid of size np[0] x np[1],
       if possible */
    if( pfft_create_procmesh_2d(MPI_COMM_WORLD, np[0], np[1],
            &comm_cart_2d) )
    {
      pfft_fprintf(MPI_COMM_WORLD, stderr,
          "Error: This test file only works with %d processes.\n",
          np[0]*np[1]);
      MPI_Finalize();
      return 1;
    }

Hereby, we use the PFFT library function

::

    void pfft_fprintf(
        MPI_Comm comm, FILE *stream, const char *format, ...);

to print the error message. This function is similar to the standard C
function with the exception, that only the process with MPI rank
:math:`0` within the given communicator will produce some output; see
Section [sec:fprintf] for details.

Inplace transforms
------------------

Similar to FFTW, PFFT is able to compute parallel FFTs completely in
place, which means that beside some constant buffers, no second data
array is necessary. Especially, the global data communication can be
performed in place. As far as we know, there is no other parallel FFT
library beside FFTW and PFFT that supports this feature. This feature is
enabled as soon as the pointer to the output array is equal to the
pointer to the input array . E.g., in Listing [lst:man\ :sub:`c`\ 2c] we
would call

::

    /* Plan parallel forward FFT */
    plan = pfft_plan_dft_3d(n, in, in, comm_cart_2d,
        PFFT_FORWARD, PFFT_TRANSPOSED_NONE);

Higher dimensional data decomposition
-------------------------------------

The test program given in Listing [lst:man\ :sub:`c`\ 2c] used a
two-dimensional data decomposition of a three-dimensional data set.
Moreover, PFFT support the computation of any :math:`d`-dimensional FFT
with :math:`r`-dimensional data decomposition as long as
:math:`r\le d-1`. For example, one can use a one-dimensional data
decomposition for any two- or higher-dimensional data set, while the
data set must be at least four-dimensional to fit to a three-dimensional
data decomposition. The case :math:`r=d` is not supported efficiently,
since during the parallel computations there is always at least one
dimension that remains local, i.e., one dimensions stays non-decomposed.
The only exception from this rule is the case :math:`d=r=3` that is
supported by PFFT in a special way, see Section [sec:3don3d] for
details.

The dimensionality of the data decomposition is given by the dimension
of the Cartesian communicator that goes into the PFFT planing routines.
Therefore, we present a generalization of communicator creation function

::

    int pfft_create_procmesh(
        int rnk_np, MPI_Comm comm, const int *np,
        MPI_Comm *comm_cart);

Hereby, the array of length gives the size of the Cartesian communicator
.

Parallel data decomposition
---------------------------

In the following, we use the notation :math:`\frac{n}{P}` to symbolize
that an array of length :math:`n` is broken into disjoint blocks and
distributed on :math:`P` MPI processes. Hereby, the data is distributed
in compliance to the FFTW-MPI data decompostion , i.e., the first
(rounded down) processes get a contiguous chunk of elements, the next
process gets the remaining data elements, and all remaining processes
get nothing. Thereby, the block size defaults to (rounded down) but can
also be user defined.

Non-transposed and transposed data layout
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the following, we use the notation :math:`\frac{n}{P}` to symbolize
that an array of length :math:`n` is distributed on :math:`P` MPI
processes. The standard PFFT data decomposition of :math:`h` interleaved
:math:`d`-dimensional arrays of equal size
:math:`n_0 \times n_1\times \hdots \times n_{d-1}` on a
:math:`r`-dimensional process mesh of size
:math:`P_0\times \hdots \times P_{r-1}` is given by the blocks

.. math:: \frac{n_0}{P_0} \times \frac{n_1}{P_1} \times \hdots \times \frac{n_{r-1}}{P_{r-1}}  \times n_r \times n_{r+1} \times \hdots \times n_{d-1} \times h.

A PFFT created with planning flag requires the inputs to be decomposed
in this standard way and produces outputs that are decomposed in the
same way.

PFFT can save half of the global communication amount, if the data
reordering to standard decomposition is omitted. The transposed data
decomposition is given by

.. math:: \frac{n_1}{P_0} \times \frac{n_2}{P_1} \times \hdots \times \frac{n_{r}}{P_{r-1}}  \times n_0 \times n_{r+1} \times \hdots \times n_{d-1} \times h

A PFFT plan created with planning flag produces outputs with transposed
data decomposition. Analogously, a PFFT plan created with planning flag
requires its inputs to be decomposed in the transposed way. Typically,
one creates a forward plan with and a backward plan with planning flag .

Note that the flags and must be passed to the array distribution
function (see Section [sec:local-size]) *as well as* to the planner (see
Section [sec:create-plan]).

Three-dimensional FFTs with three-dimensional data decomposition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Many applications work with three-dimensional block decompositions of
three-dimensional arrays. PFFT supports decompositions of the kind

.. math:: \frac{n_0}{P_0} \times \frac{n_1}{P_1} \times \frac{n_2}{P_2} \times h.

However, PFFT applies a parallel algorithms that needs at least one
non-distributed transform dimension (we do not transform along
:math:`h`), Therefore, we split the number of processes along the last
dimension into two factors :math:`P_2=Q_1Q_2`, remap the data to the
two-dimensional decomposition

.. math:: \frac{n_0}{P_0Q_0} \times \frac{n_1}{P_1Q_1} \times n_2 \times h,

and compute the parallel FFT with this two-dimensional decomposition.
Note that the 3d to 2d remap implies some very special restrictions on
the block sizes for :math:`n_0` and :math:`n_1`, i.e., the blocks must
be divisible by :math:`Q_0` and :math:`Q_1`. More precisely, the default
blocks of the 2d-decomposition are given by and (both divisions rounded
down). This implies that the default blocks of the 3d-decomposition must
be , , and (all divisions rounded down).

Planning effort
---------------

Pass one of the following flags

,

,

, or,

to the PFFT planner in order to plan all internal FFTW plans with , , ,
or , respectively. The default value is .

PFFT uses FFTW plans for parallel array transposition and the serial
transforms. In fact, every serial transform is a combination of strided
lower-dimensional FFTs and a serial array transposition (necessary to
prepare the global transposition) which can be done by a single FFTW
plan. However, it turns out that FFTW sometimes performs better if the
serial transposition and the strided FFTs are executed separately.
Therefore, PFFT introduces the flag that enables extensive run time
tests in order to find the optimal sequence of serial strided FFT and
serial transposition for every serial transform. These tests are disable
on default which corresponds to the flag .

Preserving input data
---------------------

The following flags

,

, and,

only take effect for out-of-place transforms. The first one behaves
analogously to the FFTW flag and ensures that the input values are not
overwritten. In fact, this flag implies that only the first serial
transform is executed out-of-place and all successive steps are
performed in-place on the output array. In compliance to FFTW, this is
the default behaviour for out-of-place plans.

The second flag behaves analogously to the FFTW flag and tells the
planner that the input array can be used as scratch array. This may give
some speedup for out-of-place plans, because all the intermediate
transforms and transposition steps can be performed out-of-place.

Finally, the flag can be used for out-of-place plans that store its
inputs and outputs in the same array, i.e., array is used for
intermediate out-of-place transforms and transpositions but the PFFT
inputs and outputs are stored in array .

FFTs with shifted index sets
----------------------------

Pruned FFT and Shifted Index Sets
---------------------------------

Pruned FFT
~~~~~~~~~~

For pruned r2r- and c2c-FFT are defined as

.. math:: g_l = \sum_{k=0}^{n_i-1} \hat g_k \eim{kl/n}, \quad l=0,\hdots,n_o-1,

where :math:`n_i\le n` and :math:`n_o\le n`.

Shifted Index Sets
~~~~~~~~~~~~~~~~~~

For :math:`N\in 2\N` we define the FFT with shifted inputs

For :math:`K,L,N\in 2\N`, :math:`L<N`, :math:`L<N` we define

Precisions
----------

PFFT handles multiple precisions exactly in the same way as FFTW.
Therefore, we quote part  of the FFTW manual in the context of PFFT:

You can install single and long-double precision versions of PFFT, which
replace double with float and long double, respectively; see
[sec:install]. To use these interfaces, you must

Link to the single/long-double libraries; on Unix, or instead of (or in
addition to) . (You can link to the different-precision libraries
simultaneously.)

Include the same header file.

Replace all lowercase instances of ‘’ with ‘’ or ‘’ for single or
long-double precision, respectively. ( becomes , becomes , etcetera.)

Uppercase names, i.e. names beginning with ‘’, remain the same.

Replace with or for subroutine parameters.

Ghost cell communication
------------------------

Fortran interface
-----------------

