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

#include "ipfft.h"
#include <util.h>

static INT default_block_size(INT global_array_size, int num_procs);


void PX(local_block_size_and_offset)(INT global_array_size, INT global_block_size, int which_block,
    INT *local_block_size, INT *local_block_start)
{
  /* Block index runs from blockStart <= index < blockEnd */

  *local_block_start = PX(local_block_offset)(global_array_size, global_block_size, which_block);
  *local_block_size  = PX(local_block_size)(global_array_size, global_block_size, which_block);
}


INT PX(local_block_size_shifted)(
    INT global_array_size, INT global_block_size, int shift, MPI_Comm comm_cart_1d
    )
{
  int num_procs, coord;
  MPI_Comm_size(comm_cart_1d, &num_procs);
  PX(get_mpi_cart_coord_1d)(comm_cart_1d, &coord);

  return PX(local_block_size)(global_array_size, global_block_size, PX(pos_mod)(coord+shift, num_procs));
}


INT PX(local_block_offset)(INT global_array_size, INT global_block_size, int which_block)
{
  /* For a given globalBlockSize and arraySize, compute the localBlockSize
   * and the localBlockOffset on the given process. */
  INT blocks = PX(num_blocks)(global_array_size, global_block_size);
  return (which_block >= blocks) ? 0 : which_block * global_block_size;
}


INT PX(local_block_size)(INT global_array_size, INT global_block_size, int which_block)
{
  INT blocks = PX(num_blocks)(global_array_size, global_block_size);

  if (which_block >= blocks)
    return 0;
  else
    return ((which_block == blocks - 1)
        ? (global_array_size - which_block * global_block_size) : global_block_size);
}


INT PX(global_block_size)(INT global_array_size, INT user_block_size, int num_procs)
{
  /* calculate default globalBlockSize if userBlockSize equals FFTW_MPI_DEFAULT_BLOCK */
  return (user_block_size == FFTW_MPI_DEFAULT_BLOCK) ?
      default_block_size(global_array_size, num_procs) : user_block_size;
}


INT PX(num_blocks)(INT global_array_size, INT global_block_size)
{
  return (global_array_size + global_block_size - 1) / global_block_size;
}


static INT default_block_size(INT global_array_size, int num_procs)
{
  /* Pick a default block size for dividing an array of size arrayLength among
   * numProcs processes. Divide as equally as possible, while minimizing
   * the maximum block size among the processes as well as the number of
   * processes with nonzero blocks. */
  return ((global_array_size + num_procs - 1) / num_procs);
}

