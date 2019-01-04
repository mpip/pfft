[2]ifpackageloaded#1#2 [2]ifpackageloaded#1#2 [3]ifpackageloaded#1#2#3

#1

Interface Layers of the PFFT Library
====================================

We give a quick overview of the PFFT interface layers in the order of
increasing flexibility at the example of c2c-FFTs. For r2c-, c2r-, and
r2r-FFT similar interface layer specifications apply. A full reference
list of all PFFT functions is given in Chapter [chap:ref].

Basic Interface
---------------

The ``_3d`` interface is the simplest interface layer. It is suitable
for the planning of three-dimensional FFTs.

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

Hereby, ``n``, ``local_ni``, ``local_i_start``, ``local_no``, and
``local_o_start`` are ``ptrdiff_t`` arrays of length ``3``.

The basic interface generalizes the ``_3d`` interface to FFTs of
arbitrary dimension ``rnk_n``.

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

Therefore, ``n``, ``local_ni``, ``local_i_start``, ``local_no``, and
``local_o_start`` become arrays of length ``rnk_n``.

Advanced Interface
------------------

The advanced interface introduces the arrays ``ni`` and ``no`` of length
``rnk_n`` that give the pruned FFT input and output size. Furthermore,
the arrays ``iblock`` and ``oblock`` of length ``rnk_pm`` (``rnk_pm``
being the dimension of the process mesh) serve to adjust the block size
of the input and output block decomposition. The additional parameter
``howmany`` gives the number of transforms that will be computed
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

The ``_skipped`` interface extends the ``_many`` interface by adding the
possibility to skip some of the serial FFTs.

::

    pfft_plan pfft_plan_many_dft_skipped(
        int rnk_n, const ptrdiff_t *n,
        const ptrdiff_t *ni, const ptrdiff_t *no, ptrdiff_t howmany,
        const ptrdiff_t *iblock, const ptrdiff_t *oblock,
        (red@const int *skip_trafos,@*)
        pfft_complex *in, pfft_complex *out, MPI_Comm comm_cart,
        int sign, unsigned pfft_flags);

Hereby, ``skip_trafos`` is an ``int`` array of length ``rnk_pm``\ 1+
(``rnk_pm`` being the mesh dimension of the communicator ``comm_cart``).
For ``t=0,...,rnk_pm`` set ``skip_trafos[t]=1`` if the ``t``-th serial
transformation should be computed, otherwise set ``skip_trafos[t]=0``.
Note that the local transpositions are always performed, since they are
a prerequisite for the global communication to work. At the moment it is
only possible to skip the whole serial transform along the last
``rnk_n-rnk_pm-1`` dimensions. However, this behaviour can be realized
by a call of a ``(rnk_pm``\ 1)+-dimensional PFFT with

::

    for(int t=rnk_pm+1; t<rnk_n; t++)
      howmany *= n[t];

and manual computation of the desired serial transforms along the last
``rnk_n-rnk_pm-1`` dimensions.
