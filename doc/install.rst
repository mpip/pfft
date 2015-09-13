%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\chapter{Installation and linking}\label{chap:inst}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The install of PFFT is based on the Autotools and follows the typical workflow
\begin{lstlisting}[escapechar=§]
./configure
make
make install
\end{lstlisting}


\section{Install of the latest official FFTW release}\label{sec:fftw_inst}
PFFT depends on Release~\fftwversion{} of the FFTW library~\cite{fftw}.
For the sake of completeness, we show the command line based install procedure in the following.
However, note that we provide install scripts on \websoft that simplify the install a lot.
We highly recommend to use these install scripts, since they additionally apply several
performance patches and bugfixes that have been submitted to the FFTW developers but
are not yet included in the official FFTW releases.
\begin{lstlisting}[escapechar=§]
wget http://www.fftw.org/fftw-§\fftwversionsl§.tar.gz
tar xzvf fftw-§\fftwversion§.tar.gz
cd fftw-§\fftwversion§
./configure --enable-mpi --prefix=$HOME/local/fftw3_mpi §\label{lst:fftw:conf}§
make
make install
\end{lstlisting}
The MPI algorithms of FFTW must be build with a MPI C compiler. Add the statement \code{MPICC=\$MPICCOMP}
at the end of line~\ref{lst:fftw:conf} if the \code{configure} script fails to determine the right
MPI C compiler \code{\$MPICCOMP}. Similarly, the MPI Fortran compiler \code{\$MPIFCOMP} is set by \code{MPIFC=\$MPIFCOMP}.

\section{Install of the PFFT library}\label{sec:pfft-inst}
In the simplest case, the hardware platform and the \fftw-\fftwversion{} library are
recognized by the PFFT configure script automatically, so all we have to do is
\begin{lstlisting}[escapechar=§]
wget http://www.tu-chemnitz.de/~mpip/software/pfft-§\pfftversionsl§.tar.gz
tar xzvf pfft-§\pfftversion§.tar.gz
cd pfft-§\pfftversion§
./configure
make
make check
make install
\end{lstlisting}
Hereby, the optional call \code{make check} builds the test programs.
If the \fftw-\fftwversion{} software library is already installed on your system but not found by the configure script,
you can provide the FFTW installation directory \code{\$FFTWDIR} to configure by
\begin{lstlisting}[language=bash]
./configure --with-fftw3=$FFTWDIR
\end{lstlisting}
This call implies that the FFTW header files are located in \code{\$FFTWDIR/include} and the FFTW library files are located
in \code{\$FFTWDIR/lib}. Otherwise, one should specify the FFTW include path \code{\$FFTWINC} and the FFTW library path
\code{\$FFTWLIB} separately by
\begin{lstlisting}[prebreak = {\textbackslash}]
./configure --with-fftw3-includedir=$FFTWINC --with-fftw3-libdir=$FFTWLIB
\end{lstlisting}
At the end, this is equivalent to
\begin{lstlisting}[prebreak = {\textbackslash}]
./configure CPPFLAGS=-I$FFTWINC LDFLAGS=-L$FFTWLIB
\end{lstlisting}
which is more common to experienced users of the Autotools.
To install PFFT in a user specified directory \code{\$PFFTINSTDIR} call configure with the option
\begin{lstlisting}
./configure --prefix=$PFFTINSTDIR
\end{lstlisting}
However, this option is mandatory whenever you do not have root permissions on your machine, since the default install paths of 
\code{configure} are not accessible by standard users.
The PFFT library must be built with a MPI compiler. In Section~\ref{sec:fftw_inst} we already described how to hand the right compilers to the \code{configure} script.
Some more options are
\begin{compactitem}
  \item \code[keywords=]{--enable-float}: Produces a single-precision version of PFFT (float) instead of the default double-precision (double); see \ref{sec:prec}.
  \item \code[keywords=]{--enable-long-double}: Produces a long-double precision version of PFFT (long double) instead of the default double-precision (double); see \ref{sec:prec}.
  \item \code{--disable-fortran}: Disables inclusion of Fortran wrapper routines in the standard PFFT libraries.
  \item \code{--disable-tests}: Disables build of test programs.
\end{compactitem}
For more details on the options of the \code{configure} script call
\begin{lstlisting}
./configure --help
\end{lstlisting}


\section{How to include PFFT in your program}
All programs using PFFT should include its header file
\begin{lstlisting}
#include <pfft.h>
\end{lstlisting}
This header includes the FFTW headers \code{fftw.h}, \code{fftw-mpi.h} automatically. Make sure that the compiler can find them by setting
the include flags appropriately.
You must also link to the PFFT, FFTW and FFTW-MPI libraries. On Unix, this means adding \code{-lpfft -lfftw3_mpi -lfftw3 -lm} at the end of the link command.
For example, to build \code{pfft_test.c} use the following compiler invocation
\begin{lstlisting}[prebreak = {\textbackslash}]
mpicc pfft_test.c -I$PFFTINC -I$FFTWINC -L$PFFTLIB -L$FFTWLIB -lpfft -lfftw3_mpi -lfftw3 -lm
\end{lstlisting}
Substitute \code{mpicc} by any other MPI C compiler if you like.
\code{\$PFFTINC}, \code{\$FFTWINC}, \code{\$PFFTLIB}, and \code{\$FFTWLIB} denote the PFFT and FFTW include and library paths, respectively.
If you use the install scripts mentioned in Sect.~\ref{sec:pfft-inst}, these paths will be
\begin{lstlisting}[escapechar=§,numbers=none]
PFFTINC = $HOME/local/pfft-§\pfftversion§/include
FFTWINC = $HOME/local/fftw-§\fftwversion§/include
PFFTINC = $HOME/local/pfft-§\pfftversion§/lib
FFTWINC = $HOME/local/fftw-§\fftwversion§/lib
\end{lstlisting}


