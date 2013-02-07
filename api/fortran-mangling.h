/*
 * Copyright (c) 2003, 2007-8 Matteo Frigo
 * Copyright (c) 2003, 2007-8 Massachusetts Institute of Technology
 * Copyright (c) 2010-2013 Michael Pippig
 *
 * This file is part of PFFT.
 *
 * PFFT is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PFFT is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PFFT.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef FORTRAN_MANGLING_H
#define FORTRAN_MANGLING_H 1

/* Fortran-like (e.g. as in BLAS) type prefixes for Fortran interface */
#if defined(PFFT_SINGLE)
#  define pxf(name) CONCAT(spfft_, name)
#  define PXF(NAME) CONCAT(SPFFT_, NAME)
#elif defined(PFFT_LDOUBLE)
/* FIXME: what is best?  BLAS uses D..._X, apparently.  Ugh. */
#  define pxf(name) CONCAT(lpfft_, name)
#  define PXF(NAME) CONCAT(LPFFT_, NAME)
#else
#  define pxf(name) CONCAT(dpfft_, name)
#  define PXF(NAME) CONCAT(DPFFT_, NAME)
#endif

/* If FC_FUNC is not defined and the user didn't explicitly specify
   --disable-fortran, then make our best guess at default wrappers
   (since FC_FUNC_EQUIV should not be defined in this case, we
    will use both double-underscored g77 wrappers and single- or
    non-underscored wrappers).  This saves us from dealing with
    complaints in the cases where the user failed to specify
    an Fortran compiler or wrapper detection failed for some reason. */
#if !defined(FC_FUNC) && !defined(DISABLE_FORTRAN)
#  if (defined(_WIN32) || defined(__WIN32__)) && !defined(WINDOWS_FC_MANGLING)
#    define WINDOWS_FC_MANGLING 1
#  endif
#  if defined(_AIX) || defined(__hpux) || defined(hpux)
#    define FC_FUNC(a, A) a
#  elif defined(CRAY) || defined(_CRAY) || defined(_UNICOS)
#    define FC_FUNC(a, A) A
#  else
#    define FC_FUNC(a, A) a ## _
#  endif
#  define FC_FUNC_(a, A) a ## __
#endif

#if defined(WITH_G77_WRAPPERS) && !defined(DISABLE_FORTRAN)
#  undef FC_FUNC_
#  define FC_FUNC_(a, A) a ## __
#  undef FC_FUNC_EQUIV
#endif

/* annoying Windows syntax for shared-library declarations */
#if defined(PFFT_DLL) && (defined(_WIN32) || defined(__WIN32__))
#  define PFFT_VOIDFUNC __declspec(dllexport) void
#else
#  define PFFT_VOIDFUNC void
#endif

#endif /* !FORTRAN_MANGLING_H */
