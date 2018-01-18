/*
 * Copyright (c) 2011-2013 Michael Pippig
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

/* Global infos about procmesh are only enabled in debug mode
 * Otherwise we do not use any global variables. */
#if PFFT_DEBUG_GVARS
  extern MPI_Comm *gdbg_comm_cart;
  extern int gdbg_rnk_pm;
  extern MPI_Comm *gdbg_comms_pm;
#endif

static void factorize(
    int q, 
    int *ptr_q0, int *ptr_q1);
static void factorize_equal(
    int p0, int p1, int q, 
    int *ptr_q0, int *ptr_q1);


/* TODO: think about order of in and out for T_IN */
/* TODO: implement user blocksize */

static void init_blks_comms_local_size(
    const INT *n, MPI_Comm comm_cart_2d,
    INT *iblk, INT *oblk,
    MPI_Comm *icomms, MPI_Comm *ocomms,
    INT *local_ni, INT *local_no);

static void split_comms_2dto1d(
    MPI_Comm comm_cart_2d,
    MPI_Comm *icomms, MPI_Comm *ocomms);
static void get_local_blocks_by_comms(
    const INT *n,
    const INT *iblks, const MPI_Comm *icomms,
    const INT *oblks, const MPI_Comm *ocomms,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start);
static void get_local_n_2d_by_comms(
    const INT *n, const INT *blks, const MPI_Comm *comms,
    INT *local_n);
static void get_local_start_2d_by_comms(
    const INT *n, const INT *blks, const MPI_Comm *comms,
    INT *local_start);
static void get_local_blocks_by_coords(
    const INT *n,
    const INT *iblks, const int *icoords,
    const INT *oblks, const int *ocoords,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start);
static void get_local_n_2d_by_coords(
    const INT *n, const INT *blks, const int *coords,
    INT *local_n);
static void get_local_start_2d_by_coords(
    const INT *n, const INT *blks, const int *coords,
    INT *local_start);

void PX(local_block_remap_2dto1d_transposed)(
    int rnk_n, const INT *pn, 
    MPI_Comm comm_cart_2d, int pid, 
    unsigned transp_flag, unsigned trafo_flag,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start 
    )
{
  int p0, q0, rnk_pm;
  int icoords[2], ocoords[2];
  INT iblks[2], oblks[2];

  /* remap only works for 3d data on 3d procmesh */
  if(rnk_n != 2) return;

  MPI_Cartdim_get(comm_cart_2d, &rnk_pm);
  if(rnk_pm != 2) return;

  /* Handle r2c input and c2r output like r2r. For complex data we use the C2C flag. */
  if(trafo_flag & PFFTI_TRAFO_RDFT)
    trafo_flag = PFFTI_TRAFO_R2R;

  PX(get_procmesh_dims_1d)(comm_cart_2d, &p0, &q0);
  PX(default_block_size_2dto1d)(pn, p0, q0, iblks, oblks);

  MPI_Cart_coords(comm_cart_2d, pid, 2, icoords);
  PX(coords_2dto1d)(q0, icoords, ocoords);
  ocoords[1] = 0;

  /* take care of transposed data ordering */
  if(transp_flag & PFFT_TRANSPOSED_OUT){
    get_local_blocks_by_coords(pn, iblks, icoords, oblks, ocoords,
        local_ni, local_i_start, local_no, local_o_start);
  } else {
    get_local_blocks_by_coords(pn, iblks, icoords, oblks, ocoords,
        local_no, local_o_start, local_ni, local_i_start);
  }
} 

static void free_three_comms(MPI_Comm *comms)
{
  const int num_comms = 2;
  for(int t=0; t<num_comms; ++t){
    if(MPI_COMM_NULL != comms[t]){
      MPI_Comm_free(&comms[t]);
    }
  }
}

int PX(local_size_remap_2dto1d_transposed)(
    int rnk_n, const INT *pn, INT howmany, 
    MPI_Comm comm_cart_2d, 
    unsigned transp_flag, unsigned trafo_flag,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start 
    )
{
  INT mem=1, mem_tmp;
  int p0, q0, rnk_pm;
  INT nb, nt, N0, N1, h0, h1, hm, blk0, blk1;
  INT iblk[3], oblk[3];
  MPI_Comm icomms[3], ocomms[3];
  MPI_Comm comm_q0, comm_q1;

  /* remap only works for 2d data on 2d procmesh */
  if(rnk_n != 2)
    return 0;

  MPI_Cartdim_get(comm_cart_2d, &rnk_pm);
  if(rnk_pm != 2)
    return 0;

  /* Handle r2c input and c2r output like r2r. For complex data we use the C2C flag. */
  if(trafo_flag & PFFTI_TRAFO_RDFT)
    trafo_flag = PFFTI_TRAFO_R2R;

  PX(get_procmesh_dims_1d)(comm_cart_2d, &p0, &q0);
  init_blks_comms_local_size(pn, comm_cart_2d,
      iblk, oblk, icomms, ocomms,
      local_ni, local_no);

  /* n0/p0 x n1/p1 -> n1/p1 x n0/p0 */
  nb = local_ni[0];
  nt = local_ni[1];
 
  mem_tmp = PX(local_size_sertrafo)(
        nb, 1, &nt, howmany, trafo_flag);
  mem = MAX(mem, mem_tmp);

  pfft_fprintf(MPI_COMM_WORLD, stderr, "mem_tmp local1 = %td\n", mem_tmp);
  /* n1/q0 x n0/p0 -> n1 x n0/(p0*q0)
   * for each q0 ranks, a transpose of a matrix n1 x n0/p0,
   * from divided along n1, to along n0/p0
   *  
   * */
  N0 = local_ni[0]; h0 = 1; /* n0 / p0 */
  N1 = local_no[1]; h1 = 1; /* n1 */
  blk0 = oblk[0];  /* n0/(p0*q0) */
  blk1 = iblk[1];  /* n1 / q0 */
  hm = 1; /* set hm to 1 since mem will be in units of real/complex */

  pfft_fprintf(MPI_COMM_WORLD, stderr, "sizing N0 = %td N1 = %td, blk0 = %td blk1 = %td\n", N0, N1, blk0, blk1);

  PX(split_cart_procmesh_for_2dto1d_remap_q0)(comm_cart_2d, &comm_q0);
  mem_tmp = PX(local_size_global_transp)(
      N0, N1, h0, h1, hm, blk0, blk1, comm_q0);
  mem = MAX(mem, mem_tmp);
  MPI_Comm_free(&comm_q0);

  pfft_fprintf(MPI_COMM_WORLD, stderr, "mem_tmp global = %td\n", mem_tmp);

  /* n1 x n0/(p0*q0) -> n0/(p0*q0) x n1 */
  nb = local_no[1];
  nt = local_no[0];

  mem_tmp = PX(local_size_sertrafo)(
        nb, 1, &nt, howmany, trafo_flag);
  mem = MAX(mem, mem_tmp);

  pfft_fprintf(MPI_COMM_WORLD, stderr, "mem_tmp local2 = %td\n", mem_tmp);

  /* take care of transposed data ordering */
  if(transp_flag & PFFT_TRANSPOSED_OUT){
    get_local_blocks_by_comms(pn, iblk, icomms, oblk, ocomms,
        local_ni, local_i_start, local_no, local_o_start);
  } else {
    get_local_blocks_by_comms(pn, iblk, icomms, oblk, ocomms,
        local_no, local_o_start, local_ni, local_i_start);
  }

  /* free communicators */
  free_three_comms(icomms);
  free_three_comms(ocomms);

  return mem;
} 


/* ouput is written to 'in', also for outofplace */
remap_nd_plan PX(plan_remap_2dto1d_transposed)(
    remap_nd_plan ths,
    int rnk_n, const INT *pn, INT howmany, 
    MPI_Comm comm_cart_2d, R *in_user, R *out_user, 
    unsigned transp_flag, unsigned trafo_flag,
    unsigned opt_flag, unsigned io_flag, unsigned fftw_flags 
    )
{
  int p0, p1, rnk_pm;
  INT nb, nt, N0, N1, h0, h1, hm, blk0, blk1;
  INT local_ni[3], local_no[3];
  INT iblk[3],oblk[3];
  MPI_Comm icomms[3], ocomms[3];
  MPI_Comm comm_q0, comm_q1;
  R *in=in_user, *out=out_user;

  /* remap only works for 2d data on 2d procmesh */
  if(rnk_n != 2)
    return NULL;

  MPI_Cartdim_get(comm_cart_2d, &rnk_pm);
  if(rnk_pm != 2)
    return NULL;

  /* Handle r2c input and c2r output like r2r. For complex data we use the C2C flag. */
  if(trafo_flag & PFFTI_TRAFO_RDFT)
    trafo_flag = PFFTI_TRAFO_R2R;

  PX(get_procmesh_dims_1d)(comm_cart_2d, &p0, &ths->q0);
  init_blks_comms_local_size(pn, comm_cart_2d,
      iblk, oblk, icomms, ocomms,
      local_ni, local_no);

  /* n0/p0 x n1/p1 > n1/q0 x n0/p0 */
  nb = local_ni[0];
  nt = local_ni[1];

  ths->local_transp[0] = NULL;
  ths->local_transp[1] = NULL;
  ths->global_remap[0] = NULL;
  ths->global_remap[1] = NULL;
 
  /* plan TRANSPOSED_IN in opposite direction */
  if(transp_flag & PFFT_TRANSPOSED_IN){
  } else {
  }

  /* n1/q0 x n0/p0 -> n0/(p0*q0) x n1 */
  N0 = local_ni[0]; h0 = 1; /* n0 / p0 */
  N1 = local_no[1]; h1 = 1; /* n1 */
  blk0 = oblk[0];  /* n0/(p0*q0) */
  blk1 = iblk[1];  /* n1 / q0 */
  hm = howmany * (trafo_flag & PFFTI_TRAFO_C2C ? 2 : 1);

  pfft_fprintf(MPI_COMM_WORLD, stderr, "planning, local_ni = %td %td\n", local_ni[0], local_ni[1]);

  pfft_fprintf(MPI_COMM_WORLD, stderr, "planning, local_no = %td %td\n", local_no[0], local_no[1]);
  PX(split_cart_procmesh_for_2dto1d_remap_q0)(comm_cart_2d, &comm_q0);

  if(transp_flag & PFFT_TRANSPOSED_IN) {
    R * tmp;
    /* compute in-place plans on 'out' in order to preserve inputs,
     * global transp is preferable outof place */

    ths->global_remap[0] = PX(plan_global_transp)(
        N1, N0, h1, h0, hm, blk1, blk0,
        comm_q0, in, out, PFFT_TRANSPOSED_OUT, fftw_flags);

    ths->local_transp[1] = PX(plan_sertrafo)(
        nb, 1, &nt, howmany, out, out, 0, NULL,
        trafo_flag| PFFTI_TRAFO_SKIP, transp_flag, 0,
        opt_flag, fftw_flags);

  } else {
    R * tmp;
    if(io_flag & PFFT_DESTROY_INPUT)
      tmp = in; 
    else;
      /* default: compute in-place plans on 'out' in order to preserve inputs */
      tmp = out;
    ths->local_transp[0] = PX(plan_sertrafo)(
        nb, 1, &nt, howmany, in, tmp, 0, NULL,
        trafo_flag| PFFTI_TRAFO_SKIP, transp_flag, 0,
        opt_flag, fftw_flags);

    ths->global_remap[1] = PX(plan_global_transp)(
        N0, N1, h0, h1, hm, blk0, blk1,
        comm_q0, tmp, out, PFFT_TRANSPOSED_IN, fftw_flags);

  }
  MPI_Comm_free(&comm_q0);

  pfft_fprintf(MPI_COMM_WORLD, stderr, "planning, N0 = %td N1 = %td, blk0 = %td blk1 = %td hm=%td transp_flag %d\n", N0, N1, blk0, blk1, hm, transp_flag );

  /* free communicators */
  free_three_comms(icomms);
  free_three_comms(ocomms);

  return ths;
} 

static void init_blks_comms_local_size(
    const INT *n, MPI_Comm comm_cart_2d,
    INT *iblk, INT *oblk,
    MPI_Comm *icomms, MPI_Comm *ocomms,
    INT *local_ni, INT *local_no
    )
{
  int p0, p1, q0, q1;

  PX(get_procmesh_dims_1d)(comm_cart_2d, &p0, &q0);

  PX(default_block_size_2dto1d)(n, p0, q0,
      iblk, oblk);
  split_comms_2dto1d(comm_cart_2d,
    icomms, ocomms);
  get_local_n_2d_by_comms(n, iblk, icomms,
      local_ni);
  get_local_n_2d_by_comms(n, oblk, ocomms,
      local_no);
}

static void get_local_blocks_by_comms(
    const INT *n,
    const INT *iblks, const MPI_Comm *icomms,
    const INT *oblks, const MPI_Comm *ocomms,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  get_local_n_2d_by_comms(n, iblks, icomms, local_ni);
  get_local_start_2d_by_comms(n, iblks, icomms, local_i_start);
  get_local_n_2d_by_comms(n, oblks, ocomms, local_no);
  get_local_start_2d_by_comms(n, oblks, ocomms, local_o_start);
}

static void get_local_n_2d_by_comms(
    const INT *n, const INT *blks, const MPI_Comm *comms,
    INT *local_n
    )
{
  int coord;

  for(int t=0; t<2; t++){
    PX(get_mpi_cart_coord_1d)(comms[t], &coord);
    local_n[t] = PX(local_block_size)(n[t], blks[t], coord);
  }
}

static void get_local_start_2d_by_comms(
    const INT *n, const INT *blks, const MPI_Comm *comms,
    INT *local_start
    )
{
  int coord;

  for(int t=0; t<2; t++){
    PX(get_mpi_cart_coord_1d)(comms[t], &coord);
    local_start[t] = PX(local_block_offset)(n[t], blks[t], coord);
  }
}

static void get_local_blocks_by_coords(
    const INT *n,
    const INT *iblks, const int *icoords,
    const INT *oblks, const int *ocoords,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  get_local_n_2d_by_coords(n, iblks, icoords, local_ni);
  get_local_start_2d_by_coords(n, iblks, icoords, local_i_start);
  get_local_n_2d_by_coords(n, oblks, ocoords, local_no);
  get_local_start_2d_by_coords(n, oblks, ocoords, local_o_start);
}

static void get_local_n_2d_by_coords(
    const INT *n, const INT *blks, const int *coords,
    INT *local_n
    )
{
  for(int t=0; t<2; t++)
    local_n[t] = PX(local_block_size)(n[t], blks[t], coords[t]);
}

static void get_local_start_2d_by_coords(
    const INT *n, const INT *blks, const int *coords,
    INT *local_start
    )
{
  for(int t=0; t<2; t++)
    local_start[t] = PX(local_block_offset)(n[t], blks[t], coords[t]);
}



void PX(default_block_size_2dto1d)(
    const INT *n, int p0, int q0,
    INT *iblk, INT *oblk
    )
{
  /* n0/(p0*q0) x n1 */
  oblk[0] = PX(global_block_size)(n[0], PFFT_DEFAULT_BLOCK, p0*q0);
  oblk[1] = n[1];

  /* n0/p0 x n1/q0 */
  iblk[0] = oblk[0]*q0;
  iblk[1] = PX(global_block_size)(n[1], PFFT_DEFAULT_BLOCK, q0);
}


/* allocate array of length 3 for communicators */
static void split_comms_2dto1d(
    MPI_Comm comm_cart_2d,
    MPI_Comm *icomms, MPI_Comm *ocomms
    )
{
  int ndims=1, dim_1d, period_1d, reorder=0;
  
  /* n0/p0 x n1/p1 */
  PX(split_cart_procmesh)(comm_cart_2d, icomms);

  /* n0/(p0 * q0) x n1 */
  PX(split_cart_procmesh_2dto1d_p0q0)(comm_cart_2d, &ocomms[0]);
  dim_1d=1; period_1d=1; 
  MPI_Cart_create(MPI_COMM_SELF, ndims, &dim_1d, &period_1d, reorder,
      &ocomms[1]);

}

void PX(coords_2dto1d)(
    int q0, const int *coords_2d,
    int *coords_1d
    )
{
  coords_1d[0] = coords_2d[0]*q0 + coords_2d[1];
}

void PX(split_cart_procmesh_2dto1d_p0q0)(
    MPI_Comm comm_cart_2d,
    MPI_Comm *comm_1d
    )
{
  int p0, q0=0;
  int ndims, coords_2d[3];
  int dim_1d, period_1d, reorder=0;
  int color, key;
  MPI_Comm comm;

  if( !PX(is_cart_procmesh)(comm_cart_2d) )
    return;

  MPI_Cartdim_get(comm_cart_2d, &ndims);
  if(ndims != 2)
    return;

  PX(get_mpi_cart_coords)(comm_cart_2d, ndims, coords_2d);
  PX(get_procmesh_dims_1d)(comm_cart_2d, &p0, &q0);

  /* split into 1 comms of size p0*q0 */
  key   = coords_2d[0]*q0 + coords_2d[1];
  color = 0;
  MPI_Comm_split(comm_cart_2d, color, key, &comm);

  dim_1d = p0*q0; period_1d = 1;
  MPI_Cart_create(comm, ndims=1, &dim_1d, &period_1d, reorder,
      comm_1d);

  MPI_Comm_free(&comm);
}


/* implement the splitting to create p0 comms of size q0 */
void PX(split_cart_procmesh_for_2dto1d_remap_q0)(
    MPI_Comm comm_cart_2d,
    MPI_Comm *comm_1d
    )
{
  int p0, q0=0;
  int ndims, coords_2d[2];
  int dim_1d, period_1d, reorder=0;
  int color, key;
  MPI_Comm comm;

  if( !PX(is_cart_procmesh)(comm_cart_2d) )
    return;

  MPI_Cartdim_get(comm_cart_2d, &ndims);
  if(ndims != 2)
    return;

  PX(get_mpi_cart_coords)(comm_cart_2d, ndims, coords_2d);
  PX(get_procmesh_dims_1d)(comm_cart_2d, &p0, &q0);

  /* split into p0 comms of size q0 */
  color = coords_2d[0];
  key = coords_2d[1];
  MPI_Comm_split(comm_cart_2d, color, key, &comm);

  dim_1d = q0; period_1d = 1;
  MPI_Cart_create(comm, ndims=1, &dim_1d, &period_1d, reorder,
      comm_1d);

  MPI_Comm_free(&comm);
}

void PX(get_procmesh_dims_1d)(
    MPI_Comm comm_cart_2d,
    int *p0, int *q0
    )
{
  int ndims=2, dims[2];

  PX(get_mpi_cart_dims)(comm_cart_2d, ndims, dims);
  *p0 = dims[0];
  *q0 = dims[1];
}






