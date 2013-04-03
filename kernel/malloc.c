/*
 * Copyright (c) 2003, 2007-8 Matteo Frigo
 * Copyright (c) 2003, 2007-8 Massachusetts Institute of Technology
 * Copyright (c) 2002, 2009 Jens Keiner, Stefan Kunis, Daniel Potts
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

#include "pfft.h"
#include "ipfft.h"


void *PX(malloc)(size_t n){
  return X(malloc)(n);  
}

R *PX(alloc_real)(size_t n){
  return X(alloc_real)(n);
}

C *PX(alloc_complex)(size_t n){
  return X(alloc_complex)(n);
}
  
void PX(free)(void *p){
  X(free)(p);
}


void PX(die)(
    const char *s, MPI_Comm comm
    )
{
  fflush(stdout);
  PX(fprintf)(comm, stderr, s);
  PX(cleanup)();
  MPI_Finalize();
  exit(EXIT_FAILURE);
}
