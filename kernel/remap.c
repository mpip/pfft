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

static remap_nd_plan remap_nd_mkplan(void);

int PX(needs_remap_nd)(
    int rnk_n, MPI_Comm comm_cart
    )
{
  int rnk_pm;

  MPI_Cartdim_get(comm_cart, &rnk_pm);

  return (rnk_n == rnk_pm);
}

void PX(local_block_remap_nd_transposed)(
    int rnk_n, const INT *n, 
    MPI_Comm comm_cart, int pid, 
    unsigned transp_flag, unsigned trafo_flag,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start) {

  if(rnk_n == 3) {
    PX(local_block_remap_3dto2d_transposed)(
    rnk_n, n,
    comm_cart, pid, 
    transp_flag, trafo_flag,
    local_ni, local_i_start,
    local_no, local_o_start);
    return;
  }
  if(rnk_n == 2) {
    PX(local_block_remap_2dto1d_transposed)(
    rnk_n, n,
    comm_cart, pid, 
    transp_flag, trafo_flag,
    local_ni, local_i_start,
    local_no, local_o_start);
    return;
  }
  abort();
}

int PX(local_size_remap_nd_transposed)(
    int rnk_n, const INT *n, INT howmany, 
    MPI_Comm comm_cart, 
    unsigned transp_flag, unsigned trafo_flag,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start) {
  if(rnk_n == 3) {
    return
      PX(local_size_remap_3dto2d_transposed)(
    rnk_n, n, howmany, 
    comm_cart, 
    transp_flag, trafo_flag,
    local_ni, local_i_start,
    local_no, local_o_start);
  }
  if(rnk_n == 2) {
    return
      PX(local_size_remap_2dto1d_transposed)(
    rnk_n, n, howmany, 
    comm_cart, 
    transp_flag, trafo_flag,
    local_ni, local_i_start,
    local_no, local_o_start);
  }
  abort();
}

remap_nd_plan PX(plan_remap_nd_transposed)(
    int rnk_n, const INT *n, INT howmany, 
    MPI_Comm comm_cart, R *in, R *out, 
    unsigned transp_flag, unsigned trafo_flag,
    unsigned opt_flag, unsigned io_flag, unsigned fftw_flags)
{
  if(rnk_n == 3) {
    remap_nd_plan ths = remap_nd_mkplan();

    return
      PX(plan_remap_3dto2d_transposed)(
        ths,
        rnk_n, n, howmany, 
        comm_cart, in, out, 
        transp_flag, trafo_flag,
        opt_flag, io_flag, fftw_flags);
  }
  if(rnk_n == 2) {
    remap_nd_plan ths = remap_nd_mkplan();

    return
      PX(plan_remap_2dto1d_transposed)(
        ths,
        rnk_n, n, howmany, 
        comm_cart, in, out, 
        transp_flag, trafo_flag,
        opt_flag, io_flag, fftw_flags);
  }
  abort();
}

void PX(execute_remap_nd)(
    remap_nd_plan ths, R *plannedin, R *plannedout, R *in, R *out)
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
  PX(execute_sertrafo)(ths->local_transp[0], plannedin, plannedout, in, out);

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

  if(ths->global_remap[0])
    PX(execute_gtransp)(ths->global_remap[0], plannedin, plannedout, in, out);

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
  
  if(ths->global_remap[1])
    PX(execute_gtransp)(ths->global_remap[1], plannedin, plannedout, in, out);

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
  
  PX(execute_sertrafo)(ths->local_transp[1], plannedin, plannedout, in, out);

}

static remap_nd_plan remap_nd_mkplan(
    void
    )
{
  remap_nd_plan ths = (remap_nd_plan) malloc(sizeof(remap_nd_plan_s));

  /* initialize to NULL for easy checks */
  for(int t=0; t<2; t++){
    ths->local_transp[t] = NULL;
    ths->global_remap[t] = NULL;
  }
  
  return ths;
}


void PX(remap_nd_rmplan)(
    remap_nd_plan ths
    )
{
  /* plan was already destroyed or never initialized */
  if(ths==NULL)
    return;

  for(int t=0; t<2; t++){
    PX(sertrafo_rmplan)(ths->local_transp[t]);
    if (ths->global_remap[t])
      PX(gtransp_rmplan)(ths->global_remap[t]);
  }

  /* free memory */
  free(ths);
  /* ths=NULL; would be senseless, since we can not change the pointer itself */
}


void PX(remap_nd_split_cart_procmesh)(
    const INT rnk_n,
    MPI_Comm comm_cart,
    MPI_Comm * comms_pm)
{
  if(rnk_n == 3) {
    PX(split_cart_procmesh_3dto2d_p0q0)(comm_cart,
        comms_pm + 0);
    PX(split_cart_procmesh_3dto2d_p1q1)(comm_cart,
        comms_pm + 1);
    return;
  }
  if(rnk_n == 2) {
    PX(split_cart_procmesh_2dto1d_p0q0)(comm_cart,
        comms_pm + 0);
    return;
  }
  abort();
}

void PX(remap_nd_get_coords)(
    const INT rnk_n,
    const INT pid,
    MPI_Comm comm_cart,
    int *np_pm, int *coords_pm)
{
  if(rnk_n == 3) {
    int p0, p1, q0, q1, coords_3d[3];
    PX(get_procmesh_dims_2d)(comm_cart, &p0, &p1, &q0, &q1);
    MPI_Cart_coords(comm_cart, pid, 3, coords_3d);
    PX(coords_3dto2d)(q0, q1, coords_3d, coords_pm);
    (np_pm)[0] = p0*q0; (np_pm)[1] = p1*q1;
    return;
  }
  if(rnk_n == 2) {
    int p0, q0, coords_2d[3];
    PX(get_procmesh_dims_1d)(comm_cart, &p0, &q0);
    MPI_Cart_coords(comm_cart, pid, 2, coords_2d);
    PX(coords_2dto1d)(q0, coords_2d, coords_pm);
    (np_pm)[0] = p0*q0;
    return;
  }
  abort();
}

/* used in API as a short cut to pin the block size; notice that only iblk is needed. */
void PX(remap_nd_calculate_blocks)(
    const INT rnk_n,
    const INT *n, MPI_Comm comm_cart,
    INT *iblk
    )
{
  if(rnk_n == 3) {
    int q0, q1, p0, p1;
    INT mblk[3], oblk[3];

    PX(get_procmesh_dims_2d)(comm_cart, &p0, &p1, &q0, &q1);
    PX(default_block_size_3dto2d)(n, p0, p1, q0, q1,
          iblk, mblk, oblk);
    return;
  }
  if(rnk_n == 3) {
    int q0, p0;
    INT oblk[3];

    PX(get_procmesh_dims_1d)(comm_cart, &p0, &q0);
    PX(default_block_size_2dto1d)(n, p0, q0,
          iblk, oblk);
    return;
  }
  abort();
}

