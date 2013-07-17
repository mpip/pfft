#!/bin/sh 

# Copyright (C) 2011 The ScaFaCoS project
#  
# This file is part of ScaFaCoS.
#  
# ScaFaCoS is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  
#  ScaFaCoS is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser Public License for more details.
#  
#  You should have received a copy of the GNU Lesser Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
#

# Default values for pdflatex, bibtex and makeindex
PDFLATEX=pdflatex
BIBTEX=bibtex
MAKEINDEX=makeindex

# If the incfile exists, source it
INCFILE=latexinc.sh
test -f ./$INCFILE && . ./$INCFILE

TEXFILE=manual

echo "TEXINPUTS=$TEXINPUTS"

PDFLATEX_CALL="$PDFLATEX -halt-on-error -interaction=batchmode"

echo "Running LaTeX stage 1..."

$PDFLATEX_CALL $TEXFILE
EC=$?
if test $EC -ne 0; then
    echo "ERROR: LaTeX stage 1 failed."
    echo "These are the last 20 lines of $TEXFILE.log:"
    echo "--------------------------------------------"
    tail -n 20 $TEXFILE.log
    echo "--------------------------------------------"
    test -f $TEXFILE.pdf && rm $TEXFILE.pdf
    exit $EC
fi

echo "Running bibtex..."
# Don't fail on bibtex error: it fails if no .bib-file is there!
$BIBTEX $TEXFILE

echo "Running LaTeX stage 2..."
$PDFLATEX_CALL $TEXFILE 
EC=$?
if test $EC -ne 0; then
    echo "ERROR: LaTeX stage 2 failed."
    echo "These are the last 20 lines of $TEXFILE.log:"
    echo "--------------------------------------------"
    tail -n 20 $TEXFILE.log
    echo "--------------------------------------------"
    test -f $TEXFILE.pdf && rm $TEXFILE.pdf
    exit $EC
fi

if test -e $TEXFILE.idx; then
    echo "Running makeindex..."
    $MAKEINDEX $TEXFILE || exit $?
fi

echo "Running LaTeX stage 3..."
$PDFLATEX_CALL $TEXFILE 
EC=$?
if test $EC -ne 0; then
    echo "ERROR: LaTeX stage 3 failed."
    echo "These are the last 20 lines of $TEXFILE.log:"
    echo "--------------------------------------------"
    tail -n 20 $TEXFILE.log
    echo "--------------------------------------------"
    test -f $TEXFILE.pdf && rm $TEXFILE.pdf
    exit $EC
fi

echo "See LaTeX output in $TEXFILE.log."
