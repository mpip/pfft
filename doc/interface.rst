%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\chapter{Interface Layers of the PFFT Library}\label{chap:api}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
We give a quick overview of the PFFT interface layers in the order of increasing flexibility at the example of c2c-FFTs.
For r2c-, c2r-, and r2r-FFT similar interface layer specifications apply. A full reference list of all PFFT functions is given in Chapter~\ref{chap:ref}. 
\section{Basic Interface}
The \code{_3d} interface is the simplest interface layer. It is suitable for the planning of three-dimensional FFTs.
\begin{lstlisting}
ptrdiff_t pfft_local_size_dft_3d(
    const ptrdiff_t *n, MPI_Comm comm_cart, unsigned pfft_flags,
    ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
    ptrdiff_t *local_no, ptrdiff_t *local_o_start);
void pfft_local_block_dft_3d(
    const ptrdiff_t *n, MPI_Comm comm_cart,
    int pid, unsigned pfft_flags,
    ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
    ptrdiff_t *local_no, ptrdiff_t *local_o_start);
pfft_plan pfft_plan_dft_3d(
    const ptrdiff_t *n,
    pfft_complex *in, pfft_complex *out, MPI_Comm comm_cart,
    int sign, unsigned pfft_flags);
\end{lstlisting}
Hereby, \code{n}, \code{local_ni}, \code{local_i_start}, \code{local_no}, and \code{local_o_start} are
\code{ptrdiff_t} arrays of length \code{3}.

The basic interface generalizes the \code{_3d} interface to FFTs of arbitrary dimension \code{rnk_n}.
\begin{lstlisting}
ptrdiff_t pfft_local_size_dft(
    int rnk_n, const ptrdiff_t *n,
    MPI_Comm comm_cart, unsigned pfft_flags,
    ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
    ptrdiff_t *local_no, ptrdiff_t *local_o_start);
void pfft_local_block_dft(
    int rnk_n, const ptrdiff_t *n,
    MPI_Comm comm_cart, int pid, unsigned pfft_flags,
    ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
    ptrdiff_t *local_no, ptrdiff_t *local_o_start);
pfft_plan pfft_plan_dft(
    int rnk_n, const ptrdiff_t *n,
    pfft_complex *in, pfft_complex *out, MPI_Comm comm_cart,
    int sign, unsigned pfft_flags);
\end{lstlisting}
Therefore, \code{n}, \code{local_ni}, \code{local_i_start}, \code{local_no}, and \code{local_o_start} become
arrays of length \code{rnk_n}.

\section{Advanced Interface}
The advanced interface introduces the arrays \code{ni} and \code{no} of length \code{rnk_n}
that give the pruned FFT input and output size.
Furthermore, the arrays \code{iblock} and \code{oblock} of length \code{rnk_pm} (\code{rnk_pm} being the dimension of the process mesh)
serve to adjust the block size of the input and output block decomposition.
The additional parameter \code{howmany} gives the number of transforms that will be computed simultaneously.
\begin{lstlisting}
ptrdiff_t pfft_local_size_many_dft(
    int rnk_n, const ptrdiff_t *n,
    const ptrdiff_t *ni, const ptrdiff_t *no, ptrdiff_t howmany,
    const ptrdiff_t *iblock, const ptrdiff_t *oblock,
    MPI_Comm comm_cart, unsigned pfft_flags,
    ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
    ptrdiff_t *local_no, ptrdiff_t *local_o_start);
void pfft_local_block_many_dft(
    int rnk_n, const ptrdiff_t *ni, const ptrdiff_t *no,
    const ptrdiff_t *iblock, const ptrdiff_t *oblock,
    MPI_Comm comm_cart, int pid, unsigned pfft_flags,
    ptrdiff_t *local_ni, ptrdiff_t *local_i_start,
    ptrdiff_t *local_no, ptrdiff_t *local_o_start);
pfft_plan pfft_plan_many_dft(
    int rnk_n, const ptrdiff_t *n,
    const ptrdiff_t *ni, const ptrdiff_t *no, ptrdiff_t howmany,
    const ptrdiff_t *iblock, const ptrdiff_t *oblock,
    pfft_complex *in, pfft_complex *out, MPI_Comm comm_cart,
    int sign, unsigned pfft_flags);
\end{lstlisting}


\section{Preliminary: Skip Serial Transformations}\label{sec:skip-trafo}
The \code{_skipped} interface extends the \code{_many} interface by adding the possibility to skip some of the serial FFTs.
\begin{lstlisting}
pfft_plan pfft_plan_many_dft_skipped(
    int rnk_n, const ptrdiff_t *n,
    const ptrdiff_t *ni, const ptrdiff_t *no, ptrdiff_t howmany,
    const ptrdiff_t *iblock, const ptrdiff_t *oblock,
    (red@const int *skip_trafos,@*)
    pfft_complex *in, pfft_complex *out, MPI_Comm comm_cart,
    int sign, unsigned pfft_flags);
\end{lstlisting}
Hereby, \code{skip_trafos} is an \code{int} array of length \code{rnk_pm+1} (\code{rnk_pm} being the mesh dimension of the communicator \code{comm_cart}).
For \code{t=0,...,rnk_pm} set \code{skip_trafos[t]=1} if the \code{t}-th serial transformation should be computed, otherwise set \code{skip_trafos[t]=0}.
Note that the local transpositions are always performed, since they are a prerequisite for the global communication to work.
At the moment it is only possible to skip the whole serial transform along the last \code{rnk_n-rnk_pm-1} dimensions.
However, this behaviour can be realized by a call of a \code{(rnk_pm+1)}-dimensional PFFT with
\begin{lstlisting}
for(int t=rnk_pm+1; t<rnk_n; t++)
  howmany *= n[t];
\end{lstlisting}
and manual computation of the desired serial transforms along the last \code{rnk_n-rnk_pm-1} dimensions.
