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
static void get_local_n_3d(
    const INT *n, const INT *blks, const MPI_Comm *comms,
    INT *local_n);
static void get_local_start_3d(
    const INT *n, const INT *blks, const MPI_Comm *comms,
    INT *local_start);
static remap_3dto2d_plan remap_3dto2d_mkplan(
    void);




int PX(local_size_remap_3dto2d_transposed)(
    int rnk_n, const INT *n, INT howmany, 
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
  INT iblk[3], mblk[3], oblk[3], pn[3];
  MPI_Comm icomms[3], mcomms[3], ocomms[3];
  MPI_Comm comm_q0, comm_q1;

  /* remap only works for 3d data on 3d procmesh */
  if(rnk_n != 3)
    return 0;

  MPI_Cartdim_get(comm_cart_3d, &rnk_pm);
  if(rnk_pm != 3)
    return 0;

  PX(get_procmesh_dims_2d)(comm_cart_3d, &p0, &p1, &q0, &q1);

  PX(physical_dft_size)(rnk_n, n, trafo_flag,
      pn);

  /* handle r2c and c2r like c2c */
  if(trafo_flag & PFFTI_TRAFO_RDFT)
    trafo_flag = PFFTI_TRAFO_C2C;

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
    get_local_n_3d(pn, iblk, icomms, local_ni);
    get_local_start_3d(pn, iblk, icomms, local_i_start);
    get_local_n_3d(pn, oblk, ocomms, local_no);
    get_local_start_3d(pn, oblk, ocomms, local_o_start);
  } else {
    get_local_n_3d(pn, iblk, icomms, local_no);
    get_local_start_3d(pn, iblk, icomms, local_o_start);
    get_local_n_3d(pn, oblk, ocomms, local_ni);
    get_local_start_3d(pn, oblk, ocomms, local_i_start);
  }

  /* free communicators */
  for(int t=0; t<3; t++){
    MPI_Comm_free(&icomms[t]);
    MPI_Comm_free(&mcomms[t]);
    MPI_Comm_free(&ocomms[t]);
  }

  return mem;
} 


/* ouput is written to 'in', also for outofplace */
remap_3dto2d_plan PX(plan_remap_3dto2d_transposed)(
    int rnk_n, const INT *n, INT howmany, 
    MPI_Comm comm_cart_3d, R *in_user, R *out_user, 
    unsigned transp_flag, unsigned trafo_flag,
    unsigned opt_flag, unsigned io_flag, unsigned fftw_flags 
    )
{
  int rnk_pm;
  INT nb, nt, N0, N1, h0, h1, hm, blk0, blk1;
  INT local_ni[3], local_nm[3], local_no[3];
  INT iblk[3], mblk[3], oblk[3], pn[3];
  MPI_Comm icomms[3], mcomms[3], ocomms[3];
  MPI_Comm comm_q0, comm_q1;
  R *in=in_user, *out=out_user;
  remap_3dto2d_plan ths;

  /* remap only works for 3d data on 3d procmesh */
  if(rnk_n != 3)
    return NULL;

  MPI_Cartdim_get(comm_cart_3d, &rnk_pm);
  if(rnk_pm != 3)
    return NULL;

  ths = remap_3dto2d_mkplan();

  PX(physical_dft_size)(rnk_n, n, trafo_flag,
      pn);

  /* handle r2c and c2r like c2c */
  if(trafo_flag & PFFTI_TRAFO_RDFT)
    trafo_flag = PFFTI_TRAFO_C2C;

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
  N0 = local_ni[1]; h0 = 1;
  N1 = local_nm[2]; h1 = local_ni[0];
  blk0 = mblk[1];
  blk1 = iblk[2];
  hm = howmany * (trafo_flag & PFFTI_TRAFO_C2C) ? 2 : 1;

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
  hm = howmany * (trafo_flag & PFFTI_TRAFO_C2C) ? 2 : 1;

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
  for(int t=0; t<3; t++){
    MPI_Comm_free(&icomms[t]);
    MPI_Comm_free(&mcomms[t]);
    MPI_Comm_free(&ocomms[t]);
  }

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
  get_local_n_3d(n, iblk, icomms,
      local_ni);
  get_local_n_3d(n, mblk, mcomms,
      local_nm);
  get_local_n_3d(n, oblk, ocomms,
      local_no);
}


static void get_local_n_3d(
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

static void get_local_start_3d(
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










void PX(execute_remap_3dto2d)(
    remap_3dto2d_plan ths
    )
{
  if(ths==NULL)
    return;

// #if PFFT_DEBUG_GVARS 
//   int np, myrank;
//   int np0, np1, rnk0, rnk1;
//   MPI_Comm_size(*gdbg_comm_cart, &np);
//   MPI_Comm_rank(*gdbg_comm_cart, &myrank);
//   MPI_Comm_size(gdbg_comms_pm[0], &np0);
//   MPI_Comm_size(gdbg_comms_pm[1], &np1);
//   MPI_Comm_rank(gdbg_comms_pm[0], &rnk0);
//   MPI_Comm_rank(gdbg_comms_pm[1], &rnk1);
// 
//   int dims[3], periods[3], coords[3];
//   MPI_Cart_get(*gdbg_comm_cart, 3,
//       dims, periods, coords);
//   
//   INT local_N[3], local_N_start[3];
// 
//   int p0, p1, q0, q1;
//   p0 = dims[0]; p1 = dims[1];
//   q0 = np0/p0;  q1 = np1/p1;
// 
//   int lerr, m;
// #endif
  
  /* execute all initialized plans */
  PX(execute_sertrafo)(ths->local_transp[0]);

// #if PFFT_DEBUG_GVARS 
//   local_N[0] = 512/p0; local_N_start[0] = 0;
//   local_N[1] = 512/p1; local_N_start[1] = 0;
//   local_N[2] = 512/q0/q1; local_N_start[2] = 0;
//   
//   if(!myrank) fprintf(stderr, "!!! Before 1st remap: check all coefficients !!!\n");
//   if(!myrank) fprintf(stderr, "!!! local_N=[%td, %td, %td], local_N_start = [%td, %td, %td]\n",
//       local_N[0], local_N[1], local_N[2], local_N_start[0], local_N_start[1], local_N_start[2]);
// 
//   lerr=0; 
//   m=0;
//   for(INT k2=local_N_start[2]; k2<local_N_start[2]+local_N[2]; k2++)
//     for(INT k0=local_N_start[0]; k0<local_N_start[0]+local_N[0]; k0++)
//       for(INT k1=local_N_start[1]; k1<local_N_start[1]+local_N[1]; k1++)
//         for(INT h=0; h<2; h++, m++){
//           R ind = h + 2*(k2 + 512*(k1 + 512*k0));
//           R data = ths->global_remap[0]->dbg->in[m];
//           if( (data - ind) > 1e-13){
//             if(!lerr)
//               if(!myrank)
//                 fprintf(stderr, "data[%td] = %e, ind = %e, k0=%td, k1=%td, k2=%td\n", data, m, ind, k0, k1, k2);
//             lerr = 1;
//           }
//         }
// #endif

  PX(execute_gtransp)(ths->global_remap[0]);

// #if PFFT_DEBUG_GVARS 
//   local_N[0] = 512/p0; local_N_start[0] = 0;
//   local_N[1] = 512/p1/q1; local_N_start[1] = 0;
//   local_N[2] = 512/q0; local_N_start[2] = 0;
//   
//   if(!myrank) fprintf(stderr, "!!! Before 2nd remap: check all coefficients !!!\n");
//   if(!myrank) fprintf(stderr, "!!! local_N=[%td, %td, %td], local_N_start = [%td, %td, %td]\n",
//       local_N[0], local_N[1], local_N[2], local_N_start[0], local_N_start[1], local_N_start[2]);
// 
//   lerr=0; 
//   m=0;
//   for(INT k2=local_N_start[2]; k2<local_N_start[2]+local_N[2]; k2++)
//     for(INT k0=local_N_start[0]; k0<local_N_start[0]+local_N[0]; k0++)
//       for(INT k1=local_N_start[1]; k1<local_N_start[1]+local_N[1]; k1++)
//         for(INT h=0; h<2; h++, m++){
//           R ind = h + 2*(k2 + 512*(k1 + 512*k0));
//           R data = ths->global_remap[1]->dbg->in[m];
//           if( (data - ind) > 1e-13){
//             if(!lerr)
//               if(!myrank)
//                 fprintf(stderr, "data[%td] = %e, ind = %e, k0=%td, k1=%td, k2=%td\n", data, m, ind, k0, k1, k2);
//             lerr = 1;
//           }
//         }
// #endif
  
  PX(execute_gtransp)(ths->global_remap[1]);

// #if PFFT_DEBUG_GVARS 
//   local_N[0] = 512/p0/q0; local_N_start[0] = 0;
//   local_N[1] = 512/p1/q1; local_N_start[1] = 0;
//   local_N[2] = 512; local_N_start[2] = 0;
//   
//   if(!myrank) fprintf(stderr, "!!! After 2nd remap: check all coefficients !!!\n");
//   if(!myrank) fprintf(stderr, "!!! local_N=[%td, %td, %td], local_N_start = [%td, %td, %td]\n",
//       local_N[0], local_N[1], local_N[2], local_N_start[0], local_N_start[1], local_N_start[2]);
// 
//   lerr=0; 
//   m=0;
//   for(INT k2=local_N_start[2]; k2<local_N_start[2]+local_N[2]; k2++)
//     for(INT k0=local_N_start[0]; k0<local_N_start[0]+local_N[0]; k0++)
//       for(INT k1=local_N_start[1]; k1<local_N_start[1]+local_N[1]; k1++)
//         for(INT h=0; h<2; h++, m++){
//           R ind = h + 2*(k2 + 512*(k1 + 512*k0));
//           R data = ths->global_remap[1]->dbg->out[m];
//           if( (data - ind) > 1e-13){
//             if(!lerr)
//               if(!myrank)
//                 fprintf(stderr, "data[%td] = %e, ind = %e, k0=%td, k1=%td, k2=%td\n", data, m, ind, k0, k1, k2);
//             lerr = 1;
//           }
//         }
// #endif
  
  PX(execute_sertrafo)(ths->local_transp[1]);
}

static remap_3dto2d_plan remap_3dto2d_mkplan(
    void
    )
{
  remap_3dto2d_plan ths = (remap_3dto2d_plan) malloc(sizeof(remap_3dto2d_plan_s));

  /* initialize to NULL for easy checks */
  for(int t=0; t<2; t++){
    ths->local_transp[t] = NULL;
    ths->global_remap[t] = NULL;
  }
  
  return ths;
}


void PX(remap_3dto2d_rmplan)(
    remap_3dto2d_plan ths
    )
{
  /* plan was already destroyed or never initialized */
  if(ths==NULL)
    return;

  for(int t=0; t<2; t++){
    PX(sertrafo_rmplan)(ths->local_transp[t]);
    PX(gtransp_rmplan)(ths->global_remap[t]);
  }

  /* free memory */
  free(ths);
  /* ths=NULL; would be senseless, since we can not change the pointer itself */
}












