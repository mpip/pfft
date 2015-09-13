Interface Layers of the PFFT Library
====================================

We give a quick overview of the PFFT interface layers in the order of
increasing flexibility at the example of c2c-FFTs. For r2c-, c2r-, and
r2r-FFT similar interface layer specifications apply. A full reference
list of all PFFT functions is given in Chapter [chap:ref].

Basic Interface
---------------

The interface is the simplest interface layer. It is suitable for the
planning of three-dimensional FFTs.

::

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

Hereby, , , , , and are arrays of length .

The basic interface generalizes the interface to FFTs of arbitrary
dimension .

::

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

Therefore, , , , , and become arrays of length .

Advanced Interface
------------------

The advanced interface introduces the arrays and of length that give the
pruned FFT input and output size. Furthermore, the arrays and of length
( being the dimension of the process mesh) serve to adjust the block
size of the input and output block decomposition. The additional
parameter gives the number of transforms that will be computed
simultaneously.

::

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

Preliminary: Skip Serial Transformations
----------------------------------------

The interface extends the interface by adding the possibility to skip
some of the serial FFTs.

::

    pfft_plan pfft_plan_many_dft_skipped(
        int rnk_n, const ptrdiff_t *n,
        const ptrdiff_t *ni, const ptrdiff_t *no, ptrdiff_t howmany,
        const ptrdiff_t *iblock, const ptrdiff_t *oblock,
        (red@const int *skip_trafos,@*)
        pfft_complex *in, pfft_complex *out, MPI_Comm comm_cart,
        int sign, unsigned pfft_flags);

Hereby, is an array of length ( being the mesh dimension of the
communicator ). For set if the -th serial transformation should be
computed, otherwise set . Note that the local transpositions are always
performed, since they are a prerequisite for the global communication to
work. At the moment it is only possible to skip the whole serial
transform along the last dimensions. However, this behaviour can be
realized by a call of a -dimensional PFFT with

::

    for(int t=rnk_pm+1; t<rnk_n; t++)
      howmany *= n[t];

and manual computation of the desired serial transforms along the last
dimensions.
