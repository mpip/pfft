[2]ifpackageloaded#1#2 [2]ifpackageloaded#1#2 [3]ifpackageloaded#1#2#3

#1

Advanced Features
=================

How to Deal with FFT Index Shifts in Parallel
---------------------------------------------

Let :math:`n\in2{\ensuremath{\mathbb{N}}}`. A common problem is that the
index of the FFT input and/or output array runs between
:math:`-\nicefrac n2,\dots,\nicefrac n2-1`, but the FFT library requires
them to run between :math:`0,\dots,n-1`. With serial program execution
one can easily remap the input data :math:`\hat g_k` in a way that is
suitable for the library, i.e.,

.. math:: \hat f_k := \hat g_{(k-\nicefrac n2\bmod n)}, \quad k = 0,\dots,n-1.

Similarly, one could remap the outputs of the library :math:`f_l`,
:math:`l=0,\cdots,n-1` in the opposite direction in order to get the
required outputs, i.e.,

.. math:: g_l := f_{l \bmod n}, \quad l = -\nicefrac n2,\dots,\nicefrac n2-1.

These shifts are also known as ``fftshift`` in Matlab.

However, with distributed memory these ``fftshift`` operations require
more complex data movements and result in a global communication. For
example, the first index of the array moves to the middle and,
therefore, the corresponding data move to another MPI process.
Fortunately, this communication can be avoided at the cost of little
extra computation. At the end of the section we present two PFFT library
functions that perform the necessary pre- and postprocessing for shifted
input and output index sets.

Shift with half the FFT size
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The special case of input shift :math:`k_s=-\nicefrac n2` and/or output
shift :math:`l_s=-\nicefrac n2` is supported by PFFT. User can choose to
shift the input (``PFFT_SHIFTED_IN``) and/or to shift the output
(``PFFT_SHIFTED_OUT``).

Here, we are interested in the computation of

.. math:: g_l = \sum_{k=-\nicefrac{n_i}{2}}^{\nicefrac{n_i}{2}-1} \hat g_k {\ensuremath{\mathrm{e}^{-2\pi{{\ensuremath{\text{\scriptsize{i}}}}} kl/n}}}, \quad l=-\nicefrac{n_o}{2},\dots,\nicefrac{n_o}{2}-1

with :math:`n, n_i, n_o \in 2{\ensuremath{\mathbb{N}}}` and
:math:`n>n_i`, :math:`n>n_o`.

With an index shift of :math:`\nicefrac n2` both in :math:`k` and
:math:`l` this equivalent to the computation of

.. math::

   \begin{aligned}
     g_{(l-\nicefrac{n}{2})}
     &= \sum_{k=\nicefrac{n}{2}-\nicefrac{n_i}{2}}^{\nicefrac{n}{2}+\nicefrac{n_i}{2}-1}
        \hat g_{(k-\nicefrac{n}{2})} {\ensuremath{\mathrm{e}^{-2\pi{{\ensuremath{\text{\scriptsize{i}}}}} (k-\nicefrac n2)(l-\nicefrac n2)/n}}} \\
     &= {{\ensuremath{\mathrm{e}}}}^{+\pi{\ensuremath{\text{\scriptsize{i}}}}l} 
          \sum_{k=\nicefrac{n}{2}-\nicefrac{n_i}{2}}^{\nicefrac{n}{2}+\nicefrac{n_i}{2}-1}
          \left(\hat g_{(k-\nicefrac{n}{2})}{{\ensuremath{\mathrm{e}}}}^{+\pi{\ensuremath{\text{\scriptsize{i}}}}(k-\nicefrac n2)}\right) {\ensuremath{\mathrm{e}^{-2\pi{{\ensuremath{\text{\scriptsize{i}}}}} kl/n}}} \\
     &= {{\ensuremath{\mathrm{e}}}}^{+\pi{\ensuremath{\text{\scriptsize{i}}}}(l-\nicefrac n2)} 
        \underbrace{
          \sum_{k=\nicefrac{n}{2}-\nicefrac{n_i}{2}}^{\nicefrac{n}{2}+\nicefrac{n_i}{2}-1}
          \underbrace{\left(\hat g_{(k-\nicefrac{n}{2})}{{\ensuremath{\mathrm{e}}}}^{+\pi{\ensuremath{\text{\scriptsize{i}}}}k}\right)}_{\hat f_k} {\ensuremath{\mathrm{e}^{-2\pi{{\ensuremath{\text{\scriptsize{i}}}}} kl/n}}}
        }_{f_l}\end{aligned}

for
:math:` l=\nicefrac n2-\nicefrac{n_o}{2},\dots,\nicefrac n2 +\nicefrac{n_o}{2}-1`.
Therefore, we get the following algorithm

.. math:: f_l = \sum_{k=0}^n \hat g_k {\ensuremath{\mathrm{e}^{-2\pi{{\ensuremath{\text{\scriptsize{i}}}}} kl/n}}}, \quad l=-\nicefrac{n_o}{2},\dots,\nicefrac{n_o}{2}-1

The special case :math:`k_s=-\frac{n_i}{2}, l_s=-\frac{n_o}{2}`
corresponds to the shifts the arrays ()

[1] =1.1ex For :math:`k=0,\dots,n-1` set :math:`\hat f_k = 0`. For
:math:`k=-\nicefrac{n_i}{2},\dots,\nicefrac{n_i}{2}-1` compute
:math:`\hat f_{(k+\nicefrac{n}{2})} = (-1)^{(k+\nicefrac{n}{2})} \hat g_{k}`.
For :math:`l=0,\dots,n-1` compute
:math:`f_l = \sum_{k=0}^{n} \hat f_k {\ensuremath{\mathrm{e}^{-2\pi{{\ensuremath{\text{\scriptsize{i}}}}} kl/n}}}`
using PFFT. For :math:`l=-\nicefrac{n_o}{2},\dots,\nicefrac{n_o}{2}-1`
compute :math:`g_l = (-1)^l f_{(l+n/2)} `.

Note, that this shift implies that the library deals with pruned FFTs in
a special way, i.e., half of the zeros are added at the beginning of the
inputs and the other half is added at the end.

Arbitrary shifts
~~~~~~~~~~~~~~~~

More general shifts must be done by the user.

In a more general setting, we are interested in the computation of FFTs
with shifted index sets, i.e., assume
:math:`k_s,l_s\in{\ensuremath{\mathbb{Z}}}` and compute

.. math::

   g_l = \sum_{k=k_s}^{n_i+k_s-1} \hat g_k {\ensuremath{\mathrm{e}^{-2\pi{{\ensuremath{\text{\scriptsize{i}}}}} kl/n}}},
     \quad l=l_s,\dots,n_o+l_s-1\,.

Because of the periodicity of the FFT this can be easily performed by
Alg. [alg:fftshift:sub:`t`\ ranslation].

[alg:fftshift:sub:`t`\ ranslation]

[1] =1.1ex For :math:`k=0,\dots,n_i-1` assign
:math:`\hat f_k = \hat g_{(k+k_s\bmod n_i)}`. For
:math:`l=0,\dots,n_o-1` compute
:math:`f_l = \sum_{k=0}^{n_i} \hat f_k {\ensuremath{\mathrm{e}^{-2\pi{{\ensuremath{\text{\scriptsize{i}}}}} kl/n}}}`
using PFFT. For :math:`l=0,\dots,n_o-1` assign
:math:`g_l = f_{(l-l_s\bmod n_o)}`.

However, this involves explicit data movement since the sequence of data
changes. For a our parallel data decomposition the change of data layout
requires data communication. A simple index shift results in the
computation of

.. math::

   \begin{aligned}
     g_{l+l_s}
     &=
       \sum_{k=k_s}^{n_i+k_s-1} \hat g_k {\ensuremath{\mathrm{e}^{-2\pi{{\ensuremath{\text{\scriptsize{i}}}}} k(l+l_s)/n}}}
       =
       \sum_{k=0}^{n_i-1} \hat g_{k+k_s} {\ensuremath{\mathrm{e}^{-2\pi{{\ensuremath{\text{\scriptsize{i}}}}} (k+k_s)(l+l_s)/n}}} \\
     &=
       {\ensuremath{\mathrm{e}^{-2\pi{{\ensuremath{\text{\scriptsize{i}}}}} k_sl/n}}} \sum_{k=0}^{n_i-1} \underbrace{\left(\hat g_{k+k_s}{\ensuremath{\mathrm{e}^{-2\pi{{\ensuremath{\text{\scriptsize{i}}}}} (k+k_s)l_s/n}}}\right)}_{=: \hat f_k} {\ensuremath{\mathrm{e}^{-2\pi{{\ensuremath{\text{\scriptsize{i}}}}} kl/n}}}\end{aligned}

for all :math:`l=0,\dots,n_o-1`. The resulting
Alg. [alg:fftshift:sub:`m`\ odulation] preserves the sequence of data at
the price of some extra computation.

[alg:fftshift:sub:`m`\ odulation]

[1] =1.1ex For :math:`k=0,\dots,n_i-1` compute
:math:`\hat f_k = \hat g_{(k+k_s)} {\ensuremath{\mathrm{e}^{-2\pi{{\ensuremath{\text{\scriptsize{i}}}}} (k+k_s)l_s/n}}}`.
For :math:`l=0,\dots,n_o-1` compute
:math:`f_l = \sum_{k=0}^{n_i} \hat f_k {\ensuremath{\mathrm{e}^{-2\pi{{\ensuremath{\text{\scriptsize{i}}}}} kl/n}}}`
using PFFT. For :math:`l=0,\dots,n_o-1` compute
:math:`g_{(l+l_s)} = f_l {\ensuremath{\mathrm{e}^{-2\pi{{\ensuremath{\text{\scriptsize{i}}}}} k_sl/n}}}`.

The special case :math:`k_s=-\frac{n_i}{2}, l_s=-\frac{n_o}{2}`
corresponds to the shifts the arrays ()

[1] =1.1ex For :math:`k=0,\dots,n_i-1` compute
:math:`\hat f_k = \hat g_{(k-\nicefrac{n_i}{2})} {{\ensuremath{\mathrm{e}}}}^{+\pi{\ensuremath{\text{\scriptsize{i}}}}(k-\nicefrac{n_i}{2})n_o/n}`.
For :math:`l=0,\dots,n_o-1` compute
:math:`f_l = \sum_{k=0}^{n_i} \hat f_k {\ensuremath{\mathrm{e}^{-2\pi{{\ensuremath{\text{\scriptsize{i}}}}} kl/n}}}`
using PFFT. For :math:`l=0,\dots,n_o-1` compute
:math:`g_{(l-\nicefrac{n_o}{2})} = f_l {{\ensuremath{\mathrm{e}}}}^{+\pi{\ensuremath{\text{\scriptsize{i}}}}n_i l/n}`.

Parallel pruned FFT
-------------------

Within PFFT we define a pruned FFT as

.. math:: g_l = \sum_{k=0}^{n_i-1} \hat g_{k} {\ensuremath{\mathrm{e}^{-2\pi{{\ensuremath{\text{\scriptsize{i}}}}} kl/n}}}, \quad l=0,\dots,n_o-1.

Formally, this is equivallent to the following regular size :math:`n`
FFT

.. math:: f_l = \sum_{k=0}^{n-1} \hat f_{k} {\ensuremath{\mathrm{e}^{-2\pi{{\ensuremath{\text{\scriptsize{i}}}}} kl/n}}}, \quad l=0,\dots,n,

with

.. math::

   \hat g_k := 
     \begin{cases}
     \hat f_k, &: k=0,\dots,n_1-1, \\
     0         &: k=n_i,\dots,n-1,    
     \end{cases}

and :math:`f_l := g_l`, :math:`k=0,\dots,n_o-1`. I.e., we add
:math:`n-n_i` zeros at the end of the input array and throw away
:math:`n-n_o` entries at the end of the output array.

The definition of pruned FFT changes for ``PFFT_SHIFTED_IN`` and
``PFFT_SHIFTED_OUT``.
