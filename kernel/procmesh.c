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

static void extract_comm_1d(
    int dim, MPI_Comm comm_cart,
    MPI_Comm *comm_1d);

int PX(create_procmesh)(
    int rnk, MPI_Comm comm, const int *np,
    MPI_Comm *comm_cart
    )
{
  int num_procs=1, size; 

  MPI_Comm_size(comm, &size);
  for(int t=0; t<rnk; t++)
    num_procs *= np[t];

  if(num_procs != size)
    return 1;

  int *periods = (int *) malloc(sizeof(int) * (size_t) rnk);
  for(int t=0; t<rnk; t++)
    periods[t] = 1;

  int *dims = (int*) malloc(sizeof(int) * (size_t) rnk);
  for(int t=0; t<rnk; t++)
    dims[t] = np[t];

  const int reorder=1;
  int ret = MPI_Cart_create(comm, rnk, dims, periods, reorder, comm_cart);

  free(periods); free(dims);

  return ret;
}

int PX(create_procmesh_1d)(
    MPI_Comm comm, int np0,
    MPI_Comm *comm_cart_1d
    )
{
  const int rnk=1;

  return PX(create_procmesh)(rnk, comm, &np0,
      comm_cart_1d);
}

int PX(create_procmesh_2d)(
    MPI_Comm comm, int np0, int np1,
    MPI_Comm *comm_cart_2d
    )
{
  const int rnk=2;
  int np[2];

  np[0] = np0; np[1] = np1;

  return PX(create_procmesh)(rnk, comm, np,
      comm_cart_2d);
}


int PX(is_cart_procmesh_2d)(
    MPI_Comm comm_cart_2d
    )
{
  int status, ndims;

  MPI_Topo_test(comm_cart_2d, &status);
  MPI_Cartdim_get(comm_cart_2d, &ndims);

  return ( (status == MPI_CART) && (ndims == 2) );
}

int PX(is_cart_procmesh)(
    MPI_Comm comm_cart
    )
{
  int status;

  MPI_Topo_test(comm_cart, &status);
  return (status == MPI_CART);
}

MPI_Comm PX(assure_cart_comm)(
    MPI_Comm comm
    )
{
  MPI_Comm comm_cart;

  if(PX(is_cart_procmesh)(comm)){
    MPI_Comm_dup(comm, &comm_cart);
  } else {
    int np;
    MPI_Comm_size(comm, &np);
    PX(create_procmesh_1d)(comm, np, &comm_cart);
  }

  return comm_cart;
}

/* allocate comms_1d before call of split_cart_procmesh */
void PX(split_cart_procmesh)(
    MPI_Comm comm_cart,
    MPI_Comm *comms_1d
    )
{
  int ndims;
  MPI_Cartdim_get(comm_cart, &ndims);
 
  for(int dim=0; dim<ndims; dim++)
    extract_comm_1d(dim, comm_cart,
        &comms_1d[dim]);
}

static void extract_comm_1d(
    int dim, MPI_Comm comm_cart,
    MPI_Comm *comm_1d
    )
{
  int ndims, *remain_dims;
  MPI_Cartdim_get(comm_cart, &ndims);
  
  remain_dims = (int *) malloc(sizeof(int) * (size_t) ndims);
  for(int t=0; t<ndims; t++)
    remain_dims[t] = (t==dim) ? 1 : 0;

  MPI_Cart_sub(comm_cart, remain_dims, comm_1d);
  
  free(remain_dims);
}


int PX(get_mpi_cart_coord_1d)(MPI_Comm comm_cart_1d, int *coord)
{
  int maxdims = 1;
  return PX(get_mpi_cart_coords)(comm_cart_1d, maxdims, coord);
}


int PX(get_mpi_cart_coords)(MPI_Comm comm_cart, int maxdims, int *coords)
{
  int rnk;
  MPI_Comm_rank(comm_cart, &rnk);
  return MPI_Cart_coords(comm_cart, rnk, maxdims, coords);
}



int PX(get_mpi_cart_dims)(MPI_Comm comm_cart, int maxdims, int *dims)
{
  int ret;
  int *periods = PX(malloc_int)(maxdims);
  int *coords  = PX(malloc_int)(maxdims);

  ret =  MPI_Cart_get(comm_cart, maxdims, dims, periods, coords);
  
  free(periods); free(coords);
  return ret;
}

