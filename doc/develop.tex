%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\chapter{Developers Guide}\label{chap:develop}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



\section{Search and replace patterns}

Correct alignment of pfft.h header
\begin{lstlisting}
%s/^\(    [^ ]\+[^\\]*\)  \\/  \1\\/g  
\end{lstlisting}

Expand most macros of pfft.h to generate the function reference of this manual:
\begin{lstlisting}[language=bash,prebreak=\textbackslash,]
sed -e 's/ *\\$//g' -e 's/PFFT_EXTERN //g' \
    -e 's/PX(\([^)]*\))/pfft_\1/g' -e 's/ INT/ ptrdiff_t/g' \
    -e 's/ R/ double/g' -e 's/ C/ pfft_complex/g' \
    -e 's/^  //g' pfft.h > pfft.h.expanded
\end{lstlisting}








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\chapter{ToDo}\label{chap:todo}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\begin{itemize}
  \item \code{PFFT_FORWARD} is defined as \code{FFTW_FORWARD}
  \item \code{FFTW_FORWARD} is defined as $-1$
  \item PFFT allows to chose between \code{FFTW_FORWARD} and \code{FFTW_BACKWARD}, which is not implemented by FFTW.
  \item Matlab uses the same sign convention, i.e., $-1$ for \code{fft} and $+1$ for \code{ifftn}
\end{itemize}

\section{Measuring parallel run times}
Use \code{MPI_Barrier} in front of every call to \code{pfft_} function to avoid unbalanced run times.

