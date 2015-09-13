%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\chapter{Tutorial}\label{chap:tuto}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The following chapter describes the usage of the PFFT library at the example of a simple test file in the first section,
followed by the more advanced features of PFFT in the next sections.

\section{A first parallel transform - Three-dimensional FFT with two-dimensional data decomposition}
We explain the basic steps for computing a parallel FFT with the PFFT library at the example
of the short test program given by Listing~\ref{lst:man_c2c}. This test computes a three-dimensional c2c-FFT on
a two-dimensional process mesh. The source code \code{manual_c2c_3d.c} can be found in directory \code{tests/}
of the library's source code tree. 
\lstinputlisting[numbers=left, float, caption={Minimal parallel c2c-FFT test program.}, label=lst:man_c2c]{../tests/manual_c2c_3d.c}

After initializing MPI with \code{MPI_Init} and before calling any other PFFT routine initialize
the parallel FFT computations via
\begin{lstlisting}
void pfft_init(void);
\end{lstlisting}
MPI introduces the concept of communicators to store all the topological information of the physical process layout.
PFFT requires to be called on a process mesh that corresponds to a periodic, Cartesian communicator.
We assist the user in creating such a communicator with the following routine
\begin{lstlisting}
int pfft_create_procmesh_2d(
    MPI_Comm comm, int np0, int np1,
    MPI_Comm *comm_cart_2d);
\end{lstlisting}
This routine uses the processes within the communicator \code{comm} to create a two-dimensional process
grid of size \code{np0} x \code{np1} and stores it into the Cartesian communicator \code{comm_cart_2d}.
Note that \code{comm_cart_2d} is allocated by the routine and must be freed with \code{MPI_Comm_free} after usage.
The input parameter \code{comm} is a communicator, indicating which processes will participate in the transform.
Choosing \code{comm} as \code{MPI_COMM_WORLD} implies that the FFT is computed on all available processes.

At the next step we need to know the data decomposition of the input and output array, that depends on
the array sizes, the process grid and the chosen parallel algorithm. Therefore, we call
\begin{lstlisting}
ptrdiff_t pfft_local_size_3d(
    ptrdiff_t *n, MPI_Comm comm_cart_2d, unsigned pfft_flags,
    ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
    ptrdiff_t *local_no, ptrdiff_t *local_o_start);
\end{lstlisting}
Hereby, \code{n}, \code{local_ni}, \code{local_i_start}, \code{local_no}, \code{local_o_start} are arrays of length $3$ that must be allocated.
The return value of this function equals the size of the local complex array that needs to be allocated by every process.
In most cases, this coincides with the product of the local array sizes -- but may be bigger,
whenever the parallel algorithm needs some extra storage.
The input value \code{n} gives the three-dimensional FFT size and the flag \code{pfft_flags} serves to adjust
some details of the parallel execution. For the sake of simplicity, we restrict our self to the case
\code{pfft_flags=PFFT_TRANSPOSED_NONE} for a while and explain the more sophisticated flags at a later point.
The output arrays \code{local_ni} and \code{local_i_start} give the size and the offset of the local input array
that result from the parallel block distribution of the global input array, i.e.,
every process owns the input data \code{in[k[0],k[1],k[2]]} with \code{local_i_start[t] <= k[t] < local_i_start[t] + local_ni[t]}
for \code{t=0,1,2}. Analogously, the output parameters \code{local_o_start} and \code{local_no} contain the size
and the offset of the local output array.

Afterward, the input and output arrays must be allocated. Hereby,
\begin{lstlisting}
pfft_complex* pfft_alloc_complex(size_t size);
\end{lstlisting}
is a simple wrapper of \code{fftw_alloc_complex}, which in turn allocates the memory via \code{fftw_malloc} to ensure proper alignment for SIMD.
Have a look at the FFTW user manual~\cite{fftw-align-mem} for more details on SIMD memory alignment and \code{fftw_malloc}.
Nevertheless, you can also use any other dynamic memory allocation.

The planning of a single three-dimensional parallel FFT of size \code{n[0]} x \code{n[1]} x \code{n[2]}
is done by the function
\begin{lstlisting}
pfft_plan pfft_plan_dft_3d(
    ptrdiff_t *n, pfft_complex *in, pfft_complex *out,
    MPI_Comm comm_cart_2d, int sign, unsigned pfft_flags);
\end{lstlisting}
We provide the address of the input and output array by the pointers \code{in} and \code{out},
respectively. An inplace transform is assumed if these pointers are equal.
The integer \code{sign} gives the sign in the exponential of the FFT. Possible values are \code{PFFT_FORWARD} ($-1$)
and \code{PFFT_BACKWARD} ($+1$).
Flags passed to the planner via \code{pfft\_flags} must coincide with the flags that were passed to \code{pfft_local_size_3d}.
Otherwise the data layout of the parallel execution may not match calculated local array sizes.
As return value we get a PFFT plan, some structure that stores all the information needed to perform a parallel FFT.

Once the plan is generated, we are allowed to fill the input array \code{in}. Note, that per default the planning step
\code{pfft_plan_dft_3d} will overwrite input array \code{in}. Therefore, you should not write any sensitive data into \code{in} until the plan was generated.
For simplicity, our test program makes use of the library function
\begin{lstlisting}
void pfft_init_input_complex_3d(
    ptrdiff_t *n, ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
    pfft_complex *in);
\end{lstlisting}
to fill the input array with some numbers. Alternatively, one can fill the array with a function \code{func} of choice
and the following loop that takes account of the parallel data layout
\begin{lstlisting}
ptrdiff_t m=0;
for(ptrdiff_t k0=0; k0 < local_ni[0]; k0++)
  for(ptrdiff_t k1=0; k1 < local_ni[1]; k1++)
    for(ptrdiff_t k2=0; k2 < local_ni[2]; k2++)
      in[m++] = func(k0 + local_i_start[0],
                     k1 + local_i_start[1],
                     k2 + local_i_start[2]);
\end{lstlisting}
The parallel FFT is computed when we execute the generated plan via
\begin{lstlisting}
void pfft_execute(const pfft_plan plan);
\end{lstlisting}
Now, the results can be read from \code{out} with an analogous three-dimensional loop.
If we do not want to execute another parallel FFT of the same type, we free the allocated memory of the plan with
\begin{lstlisting}
void pfft_destroy_plan(pfft_plan plan);
\end{lstlisting}
Additionally, we use
\begin{lstlisting}
int MPI_Comm_free(MPI_Comm *comm);  
\end{lstlisting}
to free the communicator allocated by \code{pfft_create_procmesh_2d} and
\begin{lstlisting}
void pfft_free(void *ptr);
\end{lstlisting}
to free memory allocated by \code{pfft_alloc_complex}.
Finally, we exit MPI via
\begin{lstlisting}
int MPI_Finalize(void);
\end{lstlisting}


\section{Porting FFTW-MPI based code to PFFT}\label{sec:porting}
\todo[inline]{finish FFTW2PFFT porting example}
We illustrate the close connection between FFTW-MPI and PFFT at a three-dimensional MPI example analogous to the example given in the FFTW manual~\cite{fftw-2dmpi}.
\lstinputlisting[numbers=left, float, caption={Minimal parallel c2c-FFT test program.}, label=lst:fftw_3don1d]{man_fftw_3don1d.tex}

Exactly the same task can be performed with PFFT as given in Listing~\ref{lst:pfft_3don1d}.
\begin{lstlisting}
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
\end{lstlisting}



\begin{compactitem}
  \item substitute \code{fftw3-mpi.h} by \code{pfft.h}
  \item substitute all prefixes \code{fftw_} and \code{fftw_mpi_} by \code{pfft_}
  \item substitute all prefixes \code{FFTW_} by \code{PFFT_}
  \item the integers \code{N}, \code{local_n0}, \code{local_0_start} become arrays of length 3
  \item \code{dft_} in \code{pfft_local_size_dft_3d}
  \item \code{pfft_local_size_dft_3d} has additional input \code{pfft_flags} and additional outputs \code{local_no}, \code{local_o_start}
  \item The loop that inits \code{data} becomes splitted along all three dimensions. We could also use 
  
  
\end{compactitem}


First, All prefixes \code{fftw_} are substituted by \code{pfft_}

Now, the changes in order to use a two-dimensional process mesh are marginal as can be seen in Listing~\ref{lst:pfft_3don2d}.
\begin{lstlisting}
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
\end{lstlisting}







\section{Errorcode for communicator creation}
As we have seen the function
\begin{lstlisting}
int pfft_create_procmesh_2d(
    MPI_Comm comm, int np0, int np1,
    MPI_Comm *comm_cart_2d);
\end{lstlisting}
creates a two-dimensional, periodic, Cartesian communicator. The \code{int} return value
(not used in Listing~\ref{lst:man_c2c}) is the forwarded error code of \code{MPI_Cart_create}.
It is equal to zero if the communicator was created successfully.
The most common error is that the number of processes within the input
communicator \code{comm} does not fit \code{np0 x np1}. In this case the Cartesian communicator
is not generated and the return value is unequal to zero. Therefore, a typical sanity check might look like
\begin{lstlisting}
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
\end{lstlisting}
Hereby, we use the PFFT library function
\begin{lstlisting}
void pfft_fprintf(
    MPI_Comm comm, FILE *stream, const char *format, ...);
\end{lstlisting}
to print the error message.
This function is similar to the standard C function \code{fprintf} with the exception, that only the process with MPI rank $0$
within the given communicator \code{comm} will produce some output; see Section~\ref{sec:fprintf} for details.

\section{Inplace transforms}
Similar to FFTW, PFFT is able to compute parallel FFTs completely in place, which means that beside some
constant buffers, no second data array is necessary. Especially, the global data communication
can be performed in place. As far as we know, there is no other parallel FFT library beside FFTW and PFFT that
supports this feature.
This feature is enabled as soon as the pointer to the output array \code{out} is equal to the pointer to the input array \code{in}.
E.g., in Listing~\ref{lst:man_c2c} we would call
\begin{lstlisting}[firstnumber=34]
/* Plan parallel forward FFT */
plan = pfft_plan_dft_3d(n, in, in, comm_cart_2d,
    PFFT_FORWARD, PFFT_TRANSPOSED_NONE);
\end{lstlisting}

\section{Higher dimensional data decomposition}
The test program given in Listing~\ref{lst:man_c2c} used a two-dimensional data decomposition of a three-dimensional data set.
Moreover, PFFT support the computation of any $d$-dimensional FFT with $r$-dimensional data decomposition
as long as $r\le d-1$. For example, one can use a one-dimensional data decomposition for any two- or higher-dimensional data set,
while the data set must be at least four-dimensional to fit to a three-dimensional data decomposition.
The case $r=d$ is not supported efficiently, since during the parallel computations
there is always at least one dimension that remains local, i.e., one dimensions stays non-decomposed.
The only exception from this rule is the case $d=r=3$ that is supported by PFFT in a special way, see Section~\ref{sec:3don3d} for details.

The dimensionality of the data decomposition is given by the dimension of the Cartesian communicator that
goes into the PFFT planing routines. Therefore, we present a generalization of communicator creation function
\begin{lstlisting}
int pfft_create_procmesh(
    int rnk_np, MPI_Comm comm, const int *np,
    MPI_Comm *comm_cart);
\end{lstlisting}
Hereby, the array \code{np} of length \code{rnk_np} gives the size of the Cartesian communicator \code{cart_comm}.

\section{Parallel data decomposition}\label{sec:par-data-decomp}
In the following, we use the notation $\frac{n}{P}$ to symbolize that an array of length $n$ is broken into disjoint blocks and distributed on $P$ MPI processes.
Hereby, the data is distributed in compliance to the FFTW-MPI data decompostion~\cite{fftw-mpi-data-distribution},
i.e., the first \code{P/block} (rounded down) processes get a contiguous chunk of \code{block} elements,
the next process gets the remaining \code{n - block * (n/block)} data elements, and all remaining processes get nothing.
Thereby, the block size \code{block} defaults to \code{n/P} (rounded down) but can also be user defined.

\subsection{Non-transposed and transposed data layout}
In the following, we use the notation $\frac{n}{P}$ to symbolize that an array of length $n$ is distributed on $P$ MPI processes.
The standard PFFT data decomposition of $h$ interleaved $d$-dimensional arrays of equal size $n_0 \times n_1\times \hdots \times n_{d-1}$
on a $r$-dimensional process mesh of size $P_0\times \hdots \times P_{r-1}$ is given by the blocks
\begin{equation*}
  \frac{n_0}{P_0} \times \frac{n_1}{P_1} \times \hdots \times \frac{n_{r-1}}{P_{r-1}}  \times n_r \times n_{r+1} \times \hdots \times n_{d-1} \times h.
\end{equation*}
A PFFT created with planning flag \code{PFFT_TRANSPOSED_NONE} requires the inputs to be decomposed in this standard way and produces
outputs that are decomposed in the same way.

PFFT can save half of the global communication amount, if the data reordering to standard decomposition is omitted. 
The transposed data decomposition is given by
\begin{equation*}
  \frac{n_1}{P_0} \times \frac{n_2}{P_1} \times \hdots \times \frac{n_{r}}{P_{r-1}}  \times n_0 \times n_{r+1} \times \hdots \times n_{d-1} \times h
\end{equation*}
A PFFT plan created with planning flag \code{PFFT_TRANSPOSED_OUT} produces outputs with transposed data decomposition.
Analogously, a PFFT plan created with planning flag \code{PFFT_TRANSPOSED_IN} requires its inputs to be decomposed in the transposed way.
Typically, one creates a forward plan with \code{PFFT_TRANSPOSED_OUT} and a backward plan with planning flag \code{PFFT_TRANSPOSED_IN}.

Note that the flags \code{PFFT_TRANSPOSED_OUT} and \code{PFFT_TRANSPOSED_IN} must be passed to the array distribution function (see Section~\ref{sec:local-size})
\emph{as well as} to the planner (see Section~\ref{sec:create-plan}).


\subsection{Three-dimensional FFTs with three-dimensional data decomposition}\label{sec:3don3d}
Many applications work with three-dimensional block decompositions of three-dimensional arrays.
PFFT supports decompositions of the kind
\begin{equation*}
  \frac{n_0}{P_0} \times \frac{n_1}{P_1} \times \frac{n_2}{P_2} \times h.
\end{equation*}
However, PFFT applies a parallel algorithms that needs at least one non-distributed transform dimension (we do not transform along $h$),
Therefore, we split the number of processes along the last dimension into two factors $P_2=Q_1Q_2$, remap
the data to the two-dimensional decomposition
\begin{equation*}
  \frac{n_0}{P_0Q_0} \times \frac{n_1}{P_1Q_1} \times n_2 \times h,
\end{equation*}
and compute the parallel FFT with this two-dimensional decomposition.
Note that the 3d to 2d remap implies some very special restrictions on the block sizes for $n_0$ and $n_1$, i.e.,
the blocks must be divisible by $Q_0$ and $Q_1$. More precisely, the default blocks of the 2d-decomposition
are given by \code{n0/(P0*Q0)} and \code{n1/(P1*Q1)} (both divisions rounded down).
This implies that the default blocks of the 3d-decomposition must be \code{n0/(P0*Q0) * Q0},
\code{n1/(P1*Q1) * Q1}, and \code{n2/(Q0*Q1)} (all divisions rounded down).


\section{Planning effort}
Pass one of the following flags
\begin{compactitem}
  \item \code{PFFT_ESTIMATE},
  \item \code{PFFT_MEASURE},
  \item \code{PFFT_PATIENT}, or,
  \item \code{PFFT_EXHAUSIVE}
\end{compactitem}
to the PFFT planner in order to plan all internal FFTW plans with \code{FFTW_ESTIMATE}, \code{FFTW_MEASURE}, \code{FFTW_PATIENT}, or \code{FFTW_EXHAUSIVE},
respectively. The default value is \code{PFFT_MEASURE}.

PFFT uses FFTW plans for parallel array transposition and the serial transforms. In fact, every serial transform is a combination of
strided lower-dimensional FFTs and a serial array transposition (necessary to prepare the global transposition) which can be done by a single FFTW plan.
However, it turns out that FFTW sometimes performs better if the serial transposition and the strided FFTs are executed separately.
Therefore, PFFT introduces the flag \code{PFFT_TUNE} that enables extensive run time tests in order to find the optimal sequence of
serial strided FFT and serial transposition for every serial transform. These tests are disable on default which corresponds to the flag \code{PFFT_NO_TUNE}.

\section{Preserving input data}
The following flags
\begin{compactitem}
  \item \code{PFFT_PRESERVE_INPUT},
  \item \code{PFFT_DESTROY_INPUT}, and,
  \item \code{PFFT_BUFFERED_INPLACE}
\end{compactitem}
only take effect for out-of-place transforms.
The first one behaves analogously to the FFTW flag \code{FFTW_PRESERVE_INPUT} and ensures that the input values are not overwritten.
In fact, this flag implies that only the first serial transform is executed out-of-place and all
successive steps are performed in-place on the output array.
In compliance to FFTW, this is the default behaviour for out-of-place plans.

The second flag behaves analogously to the FFTW flag \code{FFTW_DESTROY_INPUT} and tells the planner that
the input array can be used as scratch array. This may give some speedup for out-of-place plans,
because all the intermediate transforms and transposition steps can be performed out-of-place.

Finally, the flag \code{PFFT_BUFFERED_INPLACE} can be used for out-of-place plans that store its inputs and outputs in the same array,
i.e., array \code{out} is used for intermediate out-of-place transforms and transpositions but the PFFT inputs and outputs are stored in array \code{in}.


\section{FFTs with shifted index sets}
\todo[inline]{Describe shifted input and output}
\begin{compactitem}
  \item \code{PFFT_SHIFTED_IN}
  \item \code{PFFT_SHIFTED_OUT}
\end{compactitem}

\section{Pruned FFT and Shifted Index Sets}
\todo[inline]{Describe pruned FFT with shifted input and output}
\subsection{Pruned FFT}
For pruned r2r- and c2c-FFT are defined as
\begin{equation*}
  g_l = \sum_{k=0}^{n_i-1} \hat g_k \eim{kl/n}, \quad l=0,\hdots,n_o-1,
\end{equation*}
where $n_i\le n$ and $n_o\le n$.

\subsection{Shifted Index Sets}
For $N\in 2\N$ we define the FFT with shifted inputs


For $K,L,N\in 2\N$, $L<N$, $L<N$ we define





\section{Precisions}\label{sec:prec}
PFFT handles multiple precisions exactly in the same way as FFTW. Therefore, we quote part~\cite{fftw-prec} of the FFTW manual in the context of PFFT:

You can install single and long-double precision versions of PFFT, which replace double with float and long double, respectively; see \ref{sec:install}.
To use these interfaces, you must
\begin{compactitem}
  \item Link to the single/long-double libraries; on Unix, \code{-lpfftf} or \code{-lpfftl} instead of (or in addition to) \code{-lpfft}.
        (You can link to the different-precision libraries simultaneously.)
  \item Include the same \code{<pfft.h>} header file.
  \item Replace all lowercase instances of ‘\code{pfft_}’ with ‘\code{pfftf_}’ or ‘\code{pfftl_}’ for single or long-double precision, respectively.
        (\code{pfft_complex} becomes \code{pfftf_complex}, \code{pfft_execute} becomes \code{pfftf_execute}, etcetera.)
  \item Uppercase names, i.e. names beginning with ‘\code{PFFT_}’, remain the same.
  \item Replace \code{double} with \code{float} or \code{long double} for subroutine parameters.
\end{compactitem}

\section{Ghost cell communication}
\todo[inline]{explain ghost cell communication with a test file}

\section{Fortran interface}
\todo[inline]{explain F03 interface with a test file}
