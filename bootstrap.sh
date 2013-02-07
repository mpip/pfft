#!/bin/sh
#
# $Id$
#
# Copyright (c) 2003, 2006 Matteo Frigo
# Copyright (c) 2003, 2006 Massachusetts Institute of Technology
# Copyright (c) 2007 Jens Keiner
# Copyright (c) 2010 Michael Pippig
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 51
# Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
################################################################################
# NOTE: If you just want to build PFFT, do not use this file. Just follow the 
# installation instructions as described in the tutorial found under 
# doc/tutorial.
#
# This file is based on the bootstrap.sh script from NFFT 3.2 by Jens Keiner,
# which in turn is based on the bootstrap.sh script from FFTW 3.1.2 by
# M. Frigo and S. G. Johnson
################################################################################

# alias to allow for systems having glibtoolize
alias libtoolize=$(type -p glibtoolize libtoolize | head -1)

touch ChangeLog

echo "PLEASE IGNORE WARNINGS AND ERRORS"

rm -rf autom4te.cache
libtoolize
autoreconf --verbose --install --force

rm -f config.cache

# Add dependency tracking support for IBM C/C++ xlc/xlC Compilers
if grep 'xlc' build-aux/depcomp >/dev/null 2>&1
then :
else
  patch -b -p0 2>/dev/null <<\EOF
--- build-aux/depcomp.orig
+++ build-aux/depcomp
@@ -102,6 +102,12 @@
    depmode=msvc7
 fi
 
+if test "$depmode" = xlc; then
+   # IBM C/C++ Compilers xlc/xlC can output gcc-like dependency informations.
+   gccflag=-qmakedep=gcc,-MF
+   depmode=gcc
+fi
+
 case "$depmode" in
 gcc3)
 ## gcc 3 implements dependency tracking that does exactly what
@@ -226,6 +232,13 @@
   rm -f "$tmpdepfile"
   ;;
 
+xlc)
+  # This case exists only to let depend.m4 do its work.  It works by
+  # looking at the text of this script.  This case will never be run,
+  # since it is checked for above.
+  exit 1
+  ;;
+
 aix)
   # The C for AIX Compiler uses -M and outputs the dependencies
   # in a .u file.  In older versions, this file always lives in the
EOF
fi
