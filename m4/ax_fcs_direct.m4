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

AC_DEFUN_ONCE([AX_FCS_DIRECT_ARGS],[

])


AC_DEFUN([AX_FCS_DIRECT_TOP],[

AX_FCS_DIRECT_ARGS

#if $1 ; then
#  AC_MSG_NOTICE([do direct solver stuff in top-level configure])
#fi

])


AC_DEFUN([AX_FCS_DIRECT_SOLVER],[

AX_FCS_DIRECT_ARGS

])
