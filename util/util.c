/*
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
#include "util.h"
#include <string.h>     /* strcmp */


/* calculate block sizes from physical array size
 * the 'rnk_pm' distributed dimensions */
/*  */
void PX(decompose)(
    const INT *pn, const INT* block,
    int rnk_pm, const int *coords_pm,
    INT *local_n, INT *local_start
    )
{
  for(int t=0; t<rnk_pm; t++)
    PX(decompose_1d)(pn[t], block[t], coords_pm[t], &local_n[t], &local_start[t]);
}

void PX(decompose_1d)(
    INT pn, INT block_size, int which_block,
    INT *local_n, INT *local_n_start
    )
{
  PX(local_block_size_and_offset)(pn, block_size, which_block, local_n, local_n_start);
}

int PX(pos_mod)(int dividend, int divisor)
{
  /* wrapper for C modulo function with positive result */

  /* assure -divisor <= dividend < divisor */
  dividend %= divisor;
  return (dividend < 0) ? dividend + ABS(divisor) : dividend;
}



void PX(evaluate_user_block_size)(
    int rnk_pm, const INT *pn, const INT *block, const MPI_Comm *comms_pm,
    INT *block_intern
    )
{
  int num_procs;
  
  for(int t=0; t<rnk_pm; t++){
    MPI_Comm_size(comms_pm[t], &num_procs);
    
    if(block == PFFT_DEFAULT_BLOCKS)
      block_intern[t] = PX(global_block_size)(pn[t], PFFT_DEFAULT_BLOCK, num_procs);
    else
      block_intern[t] = PX(global_block_size)(pn[t], block[t], num_procs);
  }
}

void PX(evaluate_user_gcells)(
    int rnk_n, const INT *gc_below, const INT *gc_above,
    INT *gc_below_intern, INT *gc_above_intern
    )
{
  for(int t=0; t<rnk_n; t++){
    gc_below_intern[t] = (gc_below == PFFT_NO_GCELLS) ? 0 : gc_below[t];
    gc_above_intern[t] = (gc_above == PFFT_NO_GCELLS) ? 0 : gc_above[t];
  }
}


int PX(flag_active)(
    unsigned flags, unsigned search_flag
    )
{
  return (int) ( flags & search_flag );
}


int PX(flag_not_active)(
    unsigned flags, unsigned search_flag
    )
{
  return (int) ( (~flags) & search_flag );
}

INT PX(prod_INT)(
    int d, const INT *vec
    )
{
  INT prod, t;

  for(prod=1, t=0; t<d; t++)
    prod *= vec[t];

  return prod;
}

INT PX(sum_INT)(
    int d, const INT *vec
    )
{
  INT sum, t;

  for(sum=0, t=0; t<d; t++)
    sum += vec[t];

  return sum;
}

int PX(equal_INT)(
    int d, const INT *vec1, const INT *vec2
    )
{
  int equal, t;
  for(equal=1, t=0; t<d; t++)
    equal *= (vec1[t] == vec2[t]);

  return equal;
}

void PX(vcopy_INT)(
    int d, const INT *vec1,
    INT *vec2
    )
{
  for(int t=0; t<d; t++)
    vec2[t] = vec1[t];
}


void PX(vadd_INT)(
    int d, const INT *vec1, const INT *vec2,
    INT *sum
    )
{
  for(int t=0; t<d; t++)
    sum[t] = vec1[t] + vec2[t];
}

void PX(vsub_INT)(
    int d, const INT *vec1, const INT *vec2,
    INT *sum
    )
{
  for(int t=0; t<d; t++)
    sum[t] = vec1[t] - vec2[t];
}

void PX(vfprintf)(
    MPI_Comm comm, FILE *stream,
    const char *format, va_list ap
    )
{
  int myrank;
  MPI_Comm_rank(comm, &myrank);
  if(!myrank)
    vfprintf(stream, format, ap);
}

void PX(fprintf)(
    MPI_Comm comm, FILE *stream, const char *format, ...
    )
{
  va_list ap;
  va_start(ap,format);
  PX(vfprintf)(comm, stream, format, ap);
  va_end(ap);
}

void PX(printf)(
    MPI_Comm comm, const char *format, ...
    )
{
  va_list ap;
  va_start(ap,format);
  PX(vfprintf)(comm, stdout, format, ap);
  va_end(ap);
}


void PX(physical_dft_size)(
    int rnk_n, const INT *n, unsigned trafo_flag,
    INT *pn
    )
{
  /* phy. size of first dims equals input size */
  for(int t=0; t<rnk_n-1; t++)
    pn[t] = n[t];

  /* physical size of last dim changes for r2c or c2r */
  pn[rnk_n-1] = PX(physical_dft_size_1d)(
    n[rnk_n-1], trafo_flag);
}


INT PX(physical_dft_size_1d)(
    INT n, unsigned trafo_flag
    )
{
  return ( trafo_flag & PFFTI_TRAFO_RDFT ) ? n/2 + 1 : n;
}


INT* PX(malloc_INT)(
    size_t size
    )
{
  return (INT*) malloc(sizeof(INT) * size);
}

int* PX(malloc_int)(
    size_t size
    )
{
  return (int*) malloc(sizeof(int) * size);
}

unsigned* PX(malloc_unsigned)(
    size_t size
    )
{
  return (unsigned*) malloc(sizeof(unsigned) * size);
}

int PX(needs_3dto2d_remap)(
    int rnk_n, MPI_Comm comm_cart
    )
{
  int rnk_pm;

  MPI_Cartdim_get(comm_cart, &rnk_pm);

  return (rnk_n == 3) && (rnk_pm == 3);
}

