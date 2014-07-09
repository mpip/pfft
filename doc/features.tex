%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\chapter{Advanced Features}\label{chap:feat}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------
\section{How to Deal with FFT Index Shifts in Parallel}
%------------------------------------------------------------------------------
Let $n\in2\N$. A common problem is that the index of the FFT input and/or output array runs between $-\nicefrac n2,\hdots,\nicefrac n2-1$,
but the FFT library requires them to run between $0,\hdots,n-1$. With serial program execution one can easily remap the input data $\hat g_k$
in a way that is suitable for the library, i.e.,
\begin{equation*}
  \hat f_k := \hat g_{(k-\nicefrac n2\bmod n)}, \quad k = 0,\hdots,n-1.
\end{equation*}
Similarly, one could remap the outputs of the library $f_l$, $l=0,\cdots,n-1$ in the opposite direction in order to get the
required outputs, i.e.,
\begin{equation*}
  g_l := f_{l \bmod n}, \quad l = -\nicefrac n2,\hdots,\nicefrac n2-1.
\end{equation*}
These shifts are also known as \code{fftshift} in Matlab.

However, with distributed memory these \code{fftshift} operations require more complex data movements and result in a global communication.
For example, the first index of the array moves to the middle and, therefore, the corresponding data move to another MPI process.
Fortunately, this communication can be avoided at the cost of little extra computation.
At the end of the section we present two PFFT library functions that perform the necessary pre- and postprocessing
for shifted input and output index sets.

\subsection{Shift with half the FFT size}

The special case of input shift $k_s=-\nicefrac n2$ and/or output shift $l_s=-\nicefrac n2$ is supported by PFFT.
User can choose to shift the input (\verb+PFFT_SHIFTED_IN+) and/or to shift the output (\verb+PFFT_SHIFTED_OUT+).
\todo{this flag can be used for \code{local_size} and planning}

Here, we are interested in the computation of
\begin{equation*}
  g_l = \sum_{k=-\nicefrac{n_i}{2}}^{\nicefrac{n_i}{2}-1} \hat g_k \eim{kl/n}, \quad l=-\nicefrac{n_o}{2},\hdots,\nicefrac{n_o}{2}-1  
\end{equation*}
with $n, n_i, n_o \in 2\N$ and $n>n_i$, $n>n_o$.

With an index shift of $\nicefrac n2$ both in $k$ and $l$ this equivalent to the computation of
\begin{align*}
  g_{(l-\nicefrac{n}{2})}
  &= \sum_{k=\nicefrac{n}{2}-\nicefrac{n_i}{2}}^{\nicefrac{n}{2}+\nicefrac{n_i}{2}-1}
     \hat g_{(k-\nicefrac{n}{2})} \eim{(k-\nicefrac n2)(l-\nicefrac n2)/n} \\
  &= \e^{+\pi\ti l} 
       \sum_{k=\nicefrac{n}{2}-\nicefrac{n_i}{2}}^{\nicefrac{n}{2}+\nicefrac{n_i}{2}-1}
       \left(\hat g_{(k-\nicefrac{n}{2})}\e^{+\pi\ti (k-\nicefrac n2)}\right) \eim{kl/n} \\
  &= \e^{+\pi\ti(l-\nicefrac n2)} 
     \underbrace{
       \sum_{k=\nicefrac{n}{2}-\nicefrac{n_i}{2}}^{\nicefrac{n}{2}+\nicefrac{n_i}{2}-1}
       \underbrace{\left(\hat g_{(k-\nicefrac{n}{2})}\e^{+\pi\ti k}\right)}_{\hat f_k} \eim{kl/n}
     }_{f_l}
\end{align*}
for $ l=\nicefrac n2-\nicefrac{n_o}{2},\hdots,\nicefrac n2 +\nicefrac{n_o}{2}-1$.
Therefore, we get the following algorithm

\begin{equation*}
  f_l = \sum_{k=0}^n \hat g_k \eim{kl/n}, \quad l=-\nicefrac{n_o}{2},\hdots,\nicefrac{n_o}{2}-1  
\end{equation*}

The special case $k_s=-\frac{n_i}{2}, l_s=-\frac{n_o}{2}$ corresponds to the shifts the arrays (\textsf{FFTSHIFT})
\begin{algorithm}
  \begin{algorithmic}[1]
    \itemsep=1.1ex
    \State For $k=0,\hdots,n-1$ set $\hat f_k = 0$.
    \State For $k=-\nicefrac{n_i}{2},\hdots,\nicefrac{n_i}{2}-1$ compute $\hat f_{(k+\nicefrac{n}{2})} = (-1)^{(k+\nicefrac{n}{2})} \hat g_{k}$.
    \State For $l=0,\hdots,n-1$ compute $f_l = \sum_{k=0}^{n} \hat f_k \eim{kl/n}$ using PFFT.
    \State For $l=-\nicefrac{n_o}{2},\hdots,\nicefrac{n_o}{2}-1$ compute $g_l = (-1)^l f_{(l+n/2)} $.
  \end{algorithmic}
\end{algorithm}


Note, that this shift implies that the library deals with pruned FFTs in a special way, i.e., half of the zeros are added
at the beginning of the inputs and the other half is added at the end.






\subsection{Arbitrary shifts}
More general shifts must be done by the user.


In a more general setting, we are interested in the computation of FFTs with shifted index sets, i.e., assume $k_s,l_s\in\Z$ and compute
\begin{equation*}
  g_l = \sum_{k=k_s}^{n_i+k_s-1} \hat g_k \eim{kl/n},
  \quad l=l_s,\hdots,n_o+l_s-1\,.
\end{equation*}
Because of the periodicity of the FFT this can be easily performed by \algname~\ref{alg:fftshift_translation}.
\begin{algorithm}\label{alg:fftshift_translation}
  \begin{algorithmic}[1]
    \itemsep=1.1ex
    \State For $k=0,\hdots,n_i-1$ assign $\hat f_k = \hat g_{(k+k_s\bmod n_i)}$.
    \State For $l=0,\hdots,n_o-1$ compute $f_l = \sum_{k=0}^{n_i} \hat f_k \eim{kl/n}$ using PFFT.
    \State For $l=0,\hdots,n_o-1$ assign $g_l = f_{(l-l_s\bmod n_o)}$.
  \end{algorithmic}
  \caption{Shifted FFT with explicit data movement.}
\end{algorithm}
However, this involves explicit data movement since the sequence of data changes.
For a our parallel data decomposition the change of data layout requires data communication.
A simple index shift results in the computation of
\begin{align*}
  g_{l+l_s}
  &=
    \sum_{k=k_s}^{n_i+k_s-1} \hat g_k \eim{k(l+l_s)/n}
    =
    \sum_{k=0}^{n_i-1} \hat g_{k+k_s} \eim{(k+k_s)(l+l_s)/n} \\
  &=
    \eim{k_sl/n} \sum_{k=0}^{n_i-1} \underbrace{\left(\hat g_{k+k_s}\eim{(k+k_s)l_s/n}\right)}_{=: \hat f_k} \eim{kl/n}
\end{align*}
for all $l=0,\hdots,n_o-1$. The resulting \algname~\ref{alg:fftshift_modulation} preserves the sequence of
data at the price of some extra computation.
\begin{algorithm}\label{alg:fftshift_modulation}
  \begin{algorithmic}[1]
    \itemsep=1.1ex
    \State For $k=0,\hdots,n_i-1$ compute $\hat f_k = \hat g_{(k+k_s)} \eim{(k+k_s)l_s/n}$.
    \State For $l=0,\hdots,n_o-1$ compute $f_l = \sum_{k=0}^{n_i} \hat f_k \eim{kl/n}$ using PFFT.
    \State For $l=0,\hdots,n_o-1$ compute $g_{(l+l_s)} = f_l \eim{k_sl/n}$.
  \end{algorithmic}
  \caption{Shifted FFT without explicit data movement.}
\end{algorithm}

The special case $k_s=-\frac{n_i}{2}, l_s=-\frac{n_o}{2}$ corresponds to the shifts the arrays (\textsf{FFTSHIFT})
\begin{algorithm}
  \begin{algorithmic}[1]
    \itemsep=1.1ex
    \State For $k=0,\hdots,n_i-1$ compute $\hat f_k = \hat g_{(k-\nicefrac{n_i}{2})} \e^{+\pi\ti (k-\nicefrac{n_i}{2})n_o/n}$.
    \State For $l=0,\hdots,n_o-1$ compute $f_l = \sum_{k=0}^{n_i} \hat f_k \eim{kl/n}$ using PFFT.
    \State For $l=0,\hdots,n_o-1$ compute $g_{(l-\nicefrac{n_o}{2})} = f_l \e^{+\pi\ti n_i l/n}$.
  \end{algorithmic}
\end{algorithm}




%------------------------------------------------------------------------------
\section{Parallel pruned FFT}
%------------------------------------------------------------------------------
Within PFFT we define a pruned FFT as
\begin{equation*}
  g_l = \sum_{k=0}^{n_i-1} \hat g_{k} \eim{kl/n}, \quad l=0,\hdots,n_o-1.
\end{equation*}
Formally, this is equivallent to the following regular size $n$ FFT
\begin{equation*}
  f_l = \sum_{k=0}^{n-1} \hat f_{k} \eim{kl/n}, \quad l=0,\hdots,n,
\end{equation*}
with 
\begin{equation*}
  \hat g_k := 
  \begin{cases}
  \hat f_k, &: k=0,\hdots,n_1-1, \\
  0         &: k=n_i,\hdots,n-1,    
  \end{cases}
\end{equation*}
and $f_l := g_l$, $k=0,\hdots,n_o-1$. I.e., we add $n-n_i$ zeros at the end of the input array and throw away $n-n_o$ entries at the end of the output array.


The definition of pruned FFT changes for \code{PFFT_SHIFTED_IN} and \code{PFFT_SHIFTED_OUT}.


