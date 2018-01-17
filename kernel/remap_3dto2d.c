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

/*
 * useful table of i, o, and m.
 *
 * i : n0 / p0 x n1 / p1 x n2 / p2
 * o : n0 / (p0 * q0) x n1 / (p1 * q1) x n2 / 1
 * m : n0 / p0 x n1 / (p1 * q1) x n2 / q0
 *
 * */

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
    const INT *n, MPI_Comm comm_cart_3d,
    INT *iblk, INT *mblk, INT *oblk,
    MPI_Comm *icomms, MPI_Comm *mcomms, MPI_Comm *ocomms,
    INT *local_ni, INT *local_nm, INT *local_no);
static void split_comms_3dto2d(
    MPI_Comm comm_cart_3d,
    MPI_Comm *icomms, MPI_Comm *mcomms, MPI_Comm *ocomms);
static void get_local_blocks_by_comms(
    const INT *n,
    const INT *iblks, const MPI_Comm *icomms,
    const INT *oblks, const MPI_Comm *ocomms,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start);
static void get_local_n_3d_by_comms(
    const INT *n, const INT *blks, const MPI_Comm *comms,
    INT *local_n);
static void get_local_start_3d_by_comms(
    const INT *n, const INT *blks, const MPI_Comm *comms,
    INT *local_start);
static void get_local_blocks_by_coords(
    const INT *n,
    const INT *iblks, const int *icoords,
    const INT *oblks, const int *ocoords,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start);
static void get_local_n_3d_by_coords(
    const INT *n, const INT *blks, const int *coords,
    INT *local_n);
static void get_local_start_3d_by_coords(
    const INT *n, const INT *blks, const int *coords,
    INT *local_start);

void PX(local_block_remap_3dto2d_transposed)(
    int rnk_n, const INT *pn, 
    MPI_Comm comm_cart_3d, int pid, 
    unsigned transp_flag, unsigned trafo_flag,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start 
    )
{
  int p0, p1, q0, q1, rnk_pm;
  int icoords[3], ocoords[3];
  INT iblks[3], mblks[3], oblks[3];

  /* remap only works for 3d data on 3d procmesh */
  if(rnk_n != 3) return;

  MPI_Cartdim_get(comm_cart_3d, &rnk_pm);
  if(rnk_pm != 3) return;

  /* Handle r2c input and c2r output like r2r. For complex data we use the C2C flag. */
  if(trafo_flag & PFFTI_TRAFO_RDFT)
    trafo_flag = PFFTI_TRAFO_R2R;

  PX(get_procmesh_dims_2d)(comm_cart_3d, &p0, &p1, &q0, &q1);
  PX(default_block_size_3dto2d)(pn, p0, p1, q0, q1,
      iblks, mblks, oblks);

  MPI_Cart_coords(comm_cart_3d, pid, 3, icoords);
  PX(coords_3dto2d)(q0, q1, icoords, ocoords);
  ocoords[2] = 0;

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
  const int num_comms = 3;
  for(int t=0; t<num_comms; ++t){
    if(MPI_COMM_NULL != comms[t]){
      MPI_Comm_free(&comms[t]);
    }
  }
}

int PX(local_size_remap_3dto2d_transposed)(
    int rnk_n, const INT *pn, INT howmany, 
    MPI_Comm comm_cart_3d, 
    unsigned transp_flag, unsigned trafo_flag,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start 
    )
{
  INT mem=1, mem_tmp;
  int p0, p1, q0, q1, rnk_pm;
  INT nb, nt, N0, N1, h0, h1, hm, blk0, blk1;
  INT local_nm[3];
  INT iblk[3], mblk[3], oblk[3];
  MPI_Comm icomms[3], mcomms[3], ocomms[3];
  MPI_Comm comm_q0, comm_q1;

  /* remap only works for 3d data on 3d procmesh */
  if(rnk_n != 3)
    return 0;

  MPI_Cartdim_get(comm_cart_3d, &rnk_pm);
  if(rnk_pm != 3)
    return 0;

  /* Handle r2c input and c2r output like r2r. For complex data we use the C2C flag. */
  if(trafo_flag & PFFTI_TRAFO_RDFT)
    trafo_flag = PFFTI_TRAFO_R2R;

  PX(get_procmesh_dims_2d)(comm_cart_3d, &p0, &p1, &q0, &q1);
  init_blks_comms_local_size(pn, comm_cart_3d,
      iblk, mblk, oblk, icomms, mcomms, ocomms,
      local_ni, local_nm, local_no);

  /* n0/p0 x n1/p1 x n2/(q0*q1) -> n2/(q0*q1) x n0/p0 x n1/p1 */
  nb = local_ni[0] * local_ni[1];
  nt = local_ni[2];
 
  mem_tmp = PX(local_size_sertrafo)(
        nb, 1, &nt, howmany, trafo_flag);
  mem = MAX(mem, mem_tmp);

  /* n2/(q0*q1) x n0/p0 x n1/p1 -> n2/q0 x n0/p0 x n1/(p1*q1) */
  N0 = local_ni[1]; h0 = 1;
  N1 = local_nm[2]; h1 = local_ni[0];
  blk0 = mblk[1];
  blk1 = iblk[2];
  hm = 1; /* set hm to 1 since mem will be in units of real/complex */

  PX(split_cart_procmesh_for_3dto2d_remap_q1)(comm_cart_3d, &comm_q1);
  mem_tmp = PX(local_size_global_transp)(
      N0, N1, h0, h1, hm, blk0, blk1, comm_q1);
  mem = MAX(mem, mem_tmp);
  MPI_Comm_free(&comm_q1);

  /* n2/q0 x n0/p0 x n1/(p1*q1) -> n2 x n0/(p0*q0) x n1/(p1*q1) */
  N0 = local_nm[0]; h0 = local_nm[1];
  N1 = local_no[2]; h1 = 1;
  blk0 = oblk[0];
  blk1 = mblk[2];
  hm = 1; /* set hm to 1 since mem will be in units of real/complex */

  PX(split_cart_procmesh_for_3dto2d_remap_q0)(comm_cart_3d, &comm_q0);
  mem_tmp = PX(local_size_global_transp)(
      N0, N1, h0, h1, hm, blk0, blk1, comm_q0);
  mem = MAX(mem, mem_tmp);
  MPI_Comm_free(&comm_q0);

  /* n2 x n0/(p0*q0) x n1/(p1*q1) -> n0/(p0*q0) x n1/(p1*q1) x n2 */
  nb = local_no[2];
  nt = local_no[0] * local_no[1];
  
  mem_tmp = PX(local_size_sertrafo)(
        nb, 1, &nt, howmany, trafo_flag);
  mem = MAX(mem, mem_tmp);

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
  free_three_comms(mcomms);
  free_three_comms(ocomms);

  return mem;
} 


/* ouput is written to 'in', also for outofplace */
remap_nd_plan PX(plan_remap_3dto2d_transposed)(
    remap_nd_plan ths,
    int rnk_n, const INT *pn, INT howmany, 
    MPI_Comm comm_cart_3d, R *in_user, R *out_user, 
    unsigned transp_flag, unsigned trafo_flag,
    unsigned opt_flag, unsigned io_flag, unsigned fftw_flags 
    )
{
  int p0, p1, rnk_pm;
  INT nb, nt, N0, N1, h0, h1, hm, blk0, blk1;
  INT local_ni[3], local_nm[3], local_no[3];
  INT iblk[3], mblk[3], oblk[3];
  MPI_Comm icomms[3], mcomms[3], ocomms[3];
  MPI_Comm comm_q0, comm_q1;
  R *in=in_user, *out=out_user;

  /* remap only works for 3d data on 3d procmesh */
  if(rnk_n != 3)
    return NULL;

  MPI_Cartdim_get(comm_cart_3d, &rnk_pm);
  if(rnk_pm != 3)
    return NULL;

  /* Handle r2c input and c2r output like r2r. For complex data we use the C2C flag. */
  if(trafo_flag & PFFTI_TRAFO_RDFT)
    trafo_flag = PFFTI_TRAFO_R2R;

  PX(get_procmesh_dims_2d)(comm_cart_3d, &p0, &p1, &ths->q0, &ths->q1);
  init_blks_comms_local_size(pn, comm_cart_3d,
      iblk, mblk, oblk, icomms, mcomms, ocomms,
      local_ni, local_nm, local_no);

  /* n0/p0 x n1/p1 x n2/(q0*q1) -> n2/(q0*q1) x n0/p0 x n1/p1 */
  nb = local_ni[0] * local_ni[1];
  nt = local_ni[2];
 
  /* plan TRANSPOSED_IN in opposite direction */
  if(transp_flag & PFFT_TRANSPOSED_IN){
    if(~io_flag & PFFT_DESTROY_INPUT)
      in = out; /* default: compute in-place plans on 'out' in order to preserve inputs */
    ths->local_transp[1] = PX(plan_sertrafo)(
        nb, 1, &nt, howmany, out, in, 0, NULL,
        trafo_flag| PFFTI_TRAFO_SKIP, transp_flag, 0,
        opt_flag, fftw_flags);
  } else {
    ths->local_transp[0] = PX(plan_sertrafo)(
        nb, 1, &nt, howmany, in, out, 0, NULL,
        trafo_flag| PFFTI_TRAFO_SKIP, transp_flag, 0,
        opt_flag, fftw_flags);
    if(~io_flag & PFFT_DESTROY_INPUT)
      in = out; /* default: compute in-place plans on 'out' in order to preserve inputs */
  }

  /* n2/(q0*q1) x n0/p0 x n1/p1 -> n2/q0 x n0/p0 x n1/(p1*q1) */
  /* for each q0, we are looking at a transpose of
   * local_ni[1] x (local_nm[2] x local_ni[0]),
   * The intial partition is along (local_nm[2] x local_ni[0]),
   * by size iblk[2] x local_ni[0].
   * The final partition is along local_ni[1], by size
   * mblk[1].
   *
   * The math works out by referring to the table at the beginning of the code.
   * */
  N0 = local_ni[1]; h0 = 1;
  N1 = local_nm[2]; h1 = local_ni[0];
  blk0 = mblk[1];
  blk1 = iblk[2];
  hm = howmany * (trafo_flag & PFFTI_TRAFO_C2C ? 2 : 1);

  PX(split_cart_procmesh_for_3dto2d_remap_q1)(comm_cart_3d, &comm_q1);
  if(transp_flag & PFFT_TRANSPOSED_IN)
    ths->global_remap[1] = PX(plan_global_transp)(
        N1, N0, h1, h0, hm, blk1, blk0,
        comm_q1, in, out, PFFT_TRANSPOSED_OUT, fftw_flags);
  else
    ths->global_remap[0] = PX(plan_global_transp)(
        N0, N1, h0, h1, hm, blk0, blk1,
        comm_q1, out, in, PFFT_TRANSPOSED_IN, fftw_flags);
  MPI_Comm_free(&comm_q1);

  /* n2/q0 x n0/p0 x n1/(p1*q1) -> n2 x n0/(p0*q0) x n1/(p1*q1) */
  N0 = local_nm[0]; h0 = local_nm[1];
  N1 = local_no[2]; h1 = 1;
  blk0 = oblk[0];
  blk1 = mblk[2];
  hm = howmany * (trafo_flag & PFFTI_TRAFO_C2C ? 2 : 1);

  PX(split_cart_procmesh_for_3dto2d_remap_q0)(comm_cart_3d, &comm_q0);
  if(transp_flag & PFFT_TRANSPOSED_IN)
    ths->global_remap[0] = PX(plan_global_transp)(
        N1, N0, h1, h0, hm, blk1, blk0,
        comm_q0, out, in, PFFT_TRANSPOSED_OUT, fftw_flags);
  else
    ths->global_remap[1] = PX(plan_global_transp)(
        N0, N1, h0, h1, hm, blk0, blk1,
        comm_q0, in, out, PFFT_TRANSPOSED_IN, fftw_flags);
  MPI_Comm_free(&comm_q0);

  /* n2 x n0/(p0*q0) x n1/(p1*q1) -> n0/(p0*q0) x n1/(p1*q1) x n2 */
  nb = local_no[2];
  nt = local_no[0] * local_no[1];
  
  if(transp_flag & PFFT_TRANSPOSED_IN){
    if(~io_flag & PFFT_DESTROY_INPUT){
      /* restore pointers to 'in_user' and 'out_user' in order to compute the first step out-of-place */
      in = in_user;
    }

    ths->local_transp[0] = PX(plan_sertrafo)(
        nb, 1, &nt, howmany, in, out, 0, NULL,
        trafo_flag| PFFTI_TRAFO_SKIP, transp_flag, 0,
        opt_flag, fftw_flags);
  } else {
    ths->local_transp[1] = PX(plan_sertrafo)(
        nb, 1, &nt, howmany, out, in, 0, NULL,
        trafo_flag| PFFTI_TRAFO_SKIP, transp_flag, 0,
        opt_flag, fftw_flags);
  }

  /* free communicators */
  free_three_comms(icomms);
  free_three_comms(mcomms);
  free_three_comms(ocomms);

  return ths;
} 

static void init_blks_comms_local_size(
    const INT *n, MPI_Comm comm_cart_3d,
    INT *iblk, INT *mblk, INT *oblk,
    MPI_Comm *icomms, MPI_Comm *mcomms, MPI_Comm *ocomms,
    INT *local_ni, INT *local_nm, INT *local_no
    )
{
  int p0, p1, q0, q1;

  PX(get_procmesh_dims_2d)(comm_cart_3d, &p0, &p1, &q0, &q1);

  PX(default_block_size_3dto2d)(n, p0, p1, q0, q1,
      iblk, mblk, oblk);
  split_comms_3dto2d(comm_cart_3d,
    icomms, mcomms, ocomms);
  get_local_n_3d_by_comms(n, iblk, icomms,
      local_ni);
  get_local_n_3d_by_comms(n, mblk, mcomms,
      local_nm);
  get_local_n_3d_by_comms(n, oblk, ocomms,
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
  get_local_n_3d_by_comms(n, iblks, icomms, local_ni);
  get_local_start_3d_by_comms(n, iblks, icomms, local_i_start);
  get_local_n_3d_by_comms(n, oblks, ocomms, local_no);
  get_local_start_3d_by_comms(n, oblks, ocomms, local_o_start);
}

static void get_local_n_3d_by_comms(
    const INT *n, const INT *blks, const MPI_Comm *comms,
    INT *local_n
    )
{
  int coord;

  for(int t=0; t<3; t++){
    PX(get_mpi_cart_coord_1d)(comms[t], &coord);
    local_n[t] = PX(local_block_size)(n[t], blks[t], coord);
  }
}

static void get_local_start_3d_by_comms(
    const INT *n, const INT *blks, const MPI_Comm *comms,
    INT *local_start
    )
{
  int coord;

  for(int t=0; t<3; t++){
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
  get_local_n_3d_by_coords(n, iblks, icoords, local_ni);
  get_local_start_3d_by_coords(n, iblks, icoords, local_i_start);
  get_local_n_3d_by_coords(n, oblks, ocoords, local_no);
  get_local_start_3d_by_coords(n, oblks, ocoords, local_o_start);
}

static void get_local_n_3d_by_coords(
    const INT *n, const INT *blks, const int *coords,
    INT *local_n
    )
{
  for(int t=0; t<3; t++)
    local_n[t] = PX(local_block_size)(n[t], blks[t], coords[t]);
}

static void get_local_start_3d_by_coords(
    const INT *n, const INT *blks, const int *coords,
    INT *local_start
    )
{
  for(int t=0; t<3; t++)
    local_start[t] = PX(local_block_offset)(n[t], blks[t], coords[t]);
}



void PX(default_block_size_3dto2d)(
    const INT *n, int p0, int p1, int q0, int q1,
    INT *iblk, INT *mblk, INT *oblk
    )
{
  /* n0/(p0*q0) x n1/(p1*q1) x n2 */
  oblk[0] = PX(global_block_size)(n[0], PFFT_DEFAULT_BLOCK, p0*q0);
  oblk[1] = PX(global_block_size)(n[1], PFFT_DEFAULT_BLOCK, p1*q1);
  oblk[2] = n[2];

  /* n0/p0 x n1/p1 x n2/(q0*q1) */
  iblk[0] = oblk[0]*q0;
  iblk[1] = oblk[1]*q1;
  iblk[2] = PX(global_block_size)(n[2], PFFT_DEFAULT_BLOCK, q0*q1);
  //PX(global_block_size)(blk1[2], PFFT_DEFAULT_BLOCK, q1);

  /* n0/p0 x n1/(p1*q1) x n2/q0 */
  mblk[0] = oblk[0]*q0;
  mblk[1] = oblk[1];
  mblk[2] = iblk[2]*q1;
}


/* allocate array of length 3 for communicators */
static void split_comms_3dto2d(
    MPI_Comm comm_cart_3d,
    MPI_Comm *icomms, MPI_Comm *mcomms, MPI_Comm *ocomms
    )
{
  int ndims=1, dim_1d, period_1d, reorder=0;
  
  /* n0/p0 x n1/p1 x n2/(q0*q1) */
  PX(split_cart_procmesh)(comm_cart_3d, icomms);

  /* n0/p0 x n1/p1 x n2/(q0*q1) */
  PX(split_cart_procmesh_3dto2d_p0q0)(comm_cart_3d, &ocomms[0]);
  PX(split_cart_procmesh_3dto2d_p1q1)(comm_cart_3d, &ocomms[1]);
  dim_1d=1; period_1d=1; 
  MPI_Cart_create(MPI_COMM_SELF, ndims, &dim_1d, &period_1d, reorder,
      &ocomms[2]);

  /* n0/p0 x n1/(p1*q1) x n2/q0 */
  MPI_Comm_dup(icomms[0], &mcomms[0]);
  MPI_Comm_dup(ocomms[1], &mcomms[1]);
  PX(split_cart_procmesh_for_3dto2d_remap_q0)(comm_cart_3d, &mcomms[2]);
}

void PX(coords_3dto2d)(
    int q0, int q1, const int *coords_3d,
    int *coords_2d
    )
{
  coords_2d[0] = coords_3d[0]*q0 + coords_3d[2]/q1;
  coords_2d[1] = coords_3d[1]*q1 + coords_3d[2]%q1;
}

void PX(split_cart_procmesh_3dto2d_p0q0)(
    MPI_Comm comm_cart_3d,
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
  PX(get_procmesh_dims_2d)(comm_cart_3d, &p0, &p1, &q0, &q1);

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
    MPI_Comm comm_cart_3d,
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
  PX(get_procmesh_dims_2d)(comm_cart_3d, &p0, &p1, &q0, &q1);

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
    MPI_Comm comm_cart_3d,
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
  PX(get_procmesh_dims_2d)(comm_cart_3d, &p0, &p1, &q0, &q1);

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
    MPI_Comm comm_cart_3d,
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
  PX(get_procmesh_dims_2d)(comm_cart_3d, &p0, &p1, &q0, &q1);

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

void PX(get_procmesh_dims_2d)(
    MPI_Comm comm_cart_3d,
    int *p0, int *p1, int *q0, int *q1
    )
{
  int ndims=3, dims[3];

  PX(get_mpi_cart_dims)(comm_cart_3d, ndims, dims);
  *p0 = dims[0]; *p1 = dims[1];
//   factorize(dims[2], q0, q1);
  factorize_equal(dims[0], dims[1], dims[2], q0, q1);
}

/* factorize an integer q into q0*q1 with
 * q1 <= q0 and q0-q1 -> min */
static void factorize(
    int q, 
    int *ptr_q0, int *ptr_q1
    )
{
  for(int t = 1; t <= sqrt(q); t++)
    if(t * (q/t) == q)
      *ptr_q1 = t;

  *ptr_q0 = q / (*ptr_q1);
}

/* factorize an integer q into q0*q1 with
 * abs(p0*q0 - p1*q1) -> min */
static void factorize_equal(
    int p0, int p1, int q, 
    int *ptr_q0, int *ptr_q1
    )
{
  int q0, q1;
  int opt_q0  = 1;
  int opt_q1  = q;
  R   min_err = pfft_fabs(p0 * q - p1 * 1.0);

  for(q1 = 1; q1 <= sqrt(q); q1++){
    q0 = q/q1;
    if(q0*q1 == q){
      R err = pfft_fabs(p0*q0 - p1*q1);
      if(err < min_err){
        min_err = err;
        opt_q0 = q0;
        opt_q1 = q1;
      }
    }
  }

  *ptr_q0 = opt_q0;
  *ptr_q1 = opt_q1;
}












