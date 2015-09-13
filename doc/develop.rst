Developers Guide
================

Search and replace patterns
---------------------------

Correct alignment of pfft.h header

::

    %s/^\(    [^ ]\+[^\\]*\)  \\/  \1\\/g  

Expand most macros of pfft.h to generate the function reference of this
manual:

::

    sed -e 's/ *\\$//g' -e 's/PFFT_EXTERN //g' \
        -e 's/PX(\([^)]*\))/pfft_\1/g' -e 's/ INT/ ptrdiff_t/g' \
        -e 's/ R/ double/g' -e 's/ C/ pfft_complex/g' \
        -e 's/^  //g' pfft.h > pfft.h.expanded

ToDo
====

-  is defined as

-  is defined as :math:`-1`

-  PFFT allows to chose between and , which is not implemented by FFTW.

-  Matlab uses the same sign convention, i.e., :math:`-1` for and
   :math:`+1` for

Measuring parallel run times
----------------------------

Use in front of every call to function to avoid unbalanced run times.
