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
static void factorize(
    const INT *N, int p0, int p1, int q,
    int *ptr_q0, int *ptr_q1);


int PX(create_procmesh)(
    int rnk, MPI_Comm comm, const int *np,
    MPI_Comm *comm_cart
    )
{
  int *periods, reorder=1, num_procs=1, ret = 0, size; 
  int *dims;

  dims = (int*) malloc(sizeof(int) * (size_t) rnk);
  for(int t=0; t<rnk; t++)
    dims[t] = np[t];

  MPI_Comm_size(comm, &size);
  for(int t=0; t<rnk; t++)
    num_procs *= np[t];

  if(num_procs != size)
    return 1;

  periods = (int *) malloc(sizeof(int) * (size_t) rnk);
  *comm_cart = MPI_COMM_NULL;

  for(int t=0; t<rnk; t++)
    periods[t] = 1;

  ret = MPI_Cart_create(comm, rnk, dims, periods, reorder, comm_cart);

  free(periods); free(dims);

  return ret;
}

int PX(create_procmesh_2d)(
    MPI_Comm comm, int np0, int np1,
    MPI_Comm *comm_cart_2d
    )
{
  int np[2], rnk=2;

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


static int no_3dto2d(
    PX(plan) ths
    )
{
  if(ths->remap_3dto2d[0] != NULL)
    return 0;

  if(ths->remap_3dto2d[1] != NULL)
    return 0;

  PX(fprintf)(ths->comm_cart, stderr, "Error: PFFT plan does not use 3dto2d remap.\n");
  return 1;
}


void PX(procmesh_3dto2d)(
    PX(plan) ths, int *np_2d
    )
{
  if(no_3dto2d(ths)) return;

  np_2d[0] = ths->np[0];
  np_2d[1] = ths->np[1];
}

void PX(coords_3dto2d)(
    PX(plan) ths, int *coords_3d, int *coords_2d
    )
{
  int q0, q1;

  if(no_3dto2d(ths)) return;

  if(ths->remap_3dto2d[0] != NULL){
    q0 = ths->remap_3dto2d[0]->q0;
    q1 = ths->remap_3dto2d[0]->q1;
  } else {
    q0 = ths->remap_3dto2d[1]->q0;
    q1 = ths->remap_3dto2d[1]->q1;
  }

  coords_2d[0] = coords_3d[0]*q0 + coords_3d[2]/q1;
  coords_2d[1] = coords_3d[1]*q1 + coords_3d[2]%q1;
}

void PX(split_cart_procmesh_3dto2d_p0q0)(
    const INT *n, MPI_Comm comm_cart_3d,
    MPI_Comm *comm_1d
    )
{
  int p0, p1, q0=0, q1=0;
  int ndims, coords_3d[3];
  int dim_1d, period_1d, reorder=0;
  int color, key;
  MPI_Comm comm;

  if( !PX(is_cart_procmesh)(comm_cart_3d) )
    return;

  MPI_Cartdim_get(comm_cart_3d, &ndims);
  if(ndims != 3)
    return;

  PX(get_mpi_cart_coords)(comm_cart_3d, ndims, coords_3d);
  PX(get_procmesh_dims_2d)(n, comm_cart_3d, &p0, &p1, &q0, &q1);

  /* split into p1*q1 comms of size p0*q0 */
  key   = coords_3d[0]*q0 + coords_3d[2]/q1;
  color = coords_3d[1]*q1 + coords_3d[2]%q1;
  MPI_Comm_split(comm_cart_3d, color, key, &comm);

  dim_1d = p0*q0; period_1d = 1;
  MPI_Cart_create(comm, ndims=1, &dim_1d, &period_1d, reorder,
      comm_1d);

  MPI_Comm_free(&comm);
}


void PX(split_cart_procmesh_3dto2d_p1q1)(
    const INT *n, MPI_Comm comm_cart_3d,
    MPI_Comm *comm_1d
    )
{
  int p0, p1, q0=0, q1=0;
  int ndims, coords_3d[3];
  int dim_1d, period_1d, reorder=0;
  int color, key;
  MPI_Comm comm;

  if( !PX(is_cart_procmesh)(comm_cart_3d) )
    return;

  MPI_Cartdim_get(comm_cart_3d, &ndims);
  if(ndims != 3)
    return;

  PX(get_mpi_cart_coords)(comm_cart_3d, ndims, coords_3d);
  PX(get_procmesh_dims_2d)(n, comm_cart_3d, &p0, &p1, &q0, &q1);

  /* split into p0*q0 comms of size p1*q1 */
  color = coords_3d[0]*q0 + coords_3d[2]/q1;
  key   = coords_3d[1]*q1 + coords_3d[2]%q1;
  MPI_Comm_split(comm_cart_3d, color, key, &comm);

  dim_1d = p1*q1; period_1d = 1;
  MPI_Cart_create(comm, ndims=1, &dim_1d, &period_1d, reorder,
      comm_1d);

  MPI_Comm_free(&comm);
}


/* implement the splitting to create p0*p1*q0 comms of size q1
 * and p0*p1*q1 comms of size q0 */
void PX(split_cart_procmesh_for_3dto2d_remap_q0)(
    const INT *n, MPI_Comm comm_cart_3d,
    MPI_Comm *comm_1d
    )
{
  int p0, p1, q0=0, q1=0;
  int ndims, coords_3d[3];
  int dim_1d, period_1d, reorder=0;
  int color, key;
  MPI_Comm comm;

  if( !PX(is_cart_procmesh)(comm_cart_3d) )
    return;

  MPI_Cartdim_get(comm_cart_3d, &ndims);
  if(ndims != 3)
    return;

  PX(get_mpi_cart_coords)(comm_cart_3d, ndims, coords_3d);
  PX(get_procmesh_dims_2d)(n, comm_cart_3d, &p0, &p1, &q0, &q1);

  /* split into p0*p1*q1 comms of size q0 */
  color = coords_3d[0]*p1*q1 + coords_3d[1]*q1 + coords_3d[2]%q1;
  key = coords_3d[2]/q1;
  MPI_Comm_split(comm_cart_3d, color, key, &comm);

  dim_1d = q0; period_1d = 1;
  MPI_Cart_create(comm, ndims=1, &dim_1d, &period_1d, reorder,
      comm_1d);

  MPI_Comm_free(&comm);
}


void PX(split_cart_procmesh_for_3dto2d_remap_q1)(
    const INT *n, MPI_Comm comm_cart_3d,
    MPI_Comm *comm_1d
    )
{
  int p0, p1, q0=0, q1=0;
  int ndims, coords_3d[3];
  int dim_1d, period_1d, reorder=0;
  int color, key;
  MPI_Comm comm;

  if( !PX(is_cart_procmesh)(comm_cart_3d) )
    return;

  MPI_Cartdim_get(comm_cart_3d, &ndims);
  if(ndims != 3)
    return;

  PX(get_mpi_cart_coords)(comm_cart_3d, ndims, coords_3d);
  PX(get_procmesh_dims_2d)(n, comm_cart_3d, &p0, &p1, &q0, &q1);

  /* split into p0*p1*q0 comms of size q1 */
  color = coords_3d[0]*p1*q0 + coords_3d[1]*q0 + coords_3d[2]/q1;
  key = coords_3d[2]%q1;
//   key = coords_3d[2]/q0; /* TODO: delete this line after several tests */
  MPI_Comm_split(comm_cart_3d, color, key, &comm);

  dim_1d = q1; period_1d = 1;
  MPI_Cart_create(comm, ndims=1, &dim_1d, &period_1d, reorder,
      comm_1d);

  MPI_Comm_free(&comm);
}

/* factorize an integer q into q0*q1 with
 * q1 <= q0 and q0-q1 -> min */
// static void factorize_old(
//     int q, 
//     int *ptr_q0, int *ptr_q1
//     )
// {
//   for(int t1 = 1; t1 <= sqrt(q); t1++)
//     if(t1 * (q/t1) == q)
//       *ptr_q1 = t1;
// 
//   *ptr_q0 = q / (*ptr_q1);
// }

void PX(get_procmesh_dims_2d)(
    const INT *n, MPI_Comm comm_cart_3d,
    int *p0, int *p1, int *q0, int *q1
    )
{
  int ndims=3, dims[3];

  PX(get_mpi_cart_dims)(comm_cart_3d, ndims, dims);
  *p0 = dims[0]; *p1 = dims[1];
  factorize(n, dims[0], dims[1], dims[2], q0, q1);
}

static INT num_used_procs(
    INT N, int p
    )
{
  INT blk = PX(global_block_size)(N, PFFT_DEFAULT_BLOCK, p);
  return PX(num_blocks)(N, blk);
}

static INT used_procs(
    const INT *n, int p0, int p1, int q0, int q1
    )
{
  INT up00 = num_used_procs(n[0], p0*q0);
  INT up10 = num_used_procs(n[1], p0*q0);
  INT up11 = num_used_procs(n[1], p1*q1);
  INT up21 = num_used_procs(n[2], p1*q1);

  INT used_procs = up00 * up11;
  used_procs = MIN(used_procs, up21 * up00);
  used_procs = MIN(used_procs, up10 * up21);
 
  return used_procs; 
}

/* Factorize an integer q into q0*q1 with
 * Min{ n0/(p0*q0) * n1/(p1*q1), n2/(p1*q1) * n0/(p0*q0), n1/(p0*q0) * n2/(p1*q1) } -> Max
 * */
static void factorize(
    const INT *n, int p0, int p1, int q,
    int *ptr_q0, int *ptr_q1
    )
{
  INT up, max_used_procs = 0;

  for(int t = 1; t <= sqrt(q); t++){
    if(t * (q/t) == q){
      
      up = used_procs(n, p0, p1, q/t, t);
      if(up > max_used_procs){
        max_used_procs = up;
        *ptr_q0 = q/t;
      }

      up = used_procs(n, p0, p1, t, q/t);
      if(up > max_used_procs){
        max_used_procs = up;
        *ptr_q0 = t;
      }
    }
  }

  *ptr_q1 = q / (*ptr_q0);

}

