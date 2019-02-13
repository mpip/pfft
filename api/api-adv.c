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

INT PX(local_size_many_dft)(
    int rnk_n, const INT *n, const INT *ni, const INT *no,
    INT howmany, const INT *iblock, const INT *oblock,
    MPI_Comm comm_cart, unsigned pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  return PX(local_size_partrafo)(
      rnk_n, n, ni, no, howmany, iblock, oblock,
      comm_cart, PFFTI_TRAFO_C2C, pfft_flags,
      local_ni, local_i_start, local_no, local_o_start);
}

INT PX(local_size_many_dft_r2c)(
    int rnk_n, const INT *n, const INT *ni, const INT *no,
    INT howmany, const INT *iblock, const INT *oblock,
    MPI_Comm comm_cart, unsigned pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  unsigned trafo_flag = (pfft_flags & PFFT_PADDED_R2C) ? PFFTI_TRAFO_R2C_PADDED : PFFTI_TRAFO_R2C;

  return PX(local_size_partrafo)(
      rnk_n, n, ni, no, howmany, iblock, oblock,
      comm_cart, trafo_flag, pfft_flags,
      local_ni, local_i_start, local_no, local_o_start);
}

INT PX(local_size_many_dft_c2r)(
    int rnk_n, const INT *n, const INT *ni, const INT *no,
    INT howmany, const INT *iblock, const INT *oblock,
    MPI_Comm comm_cart, unsigned pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  unsigned trafo_flag = (pfft_flags & PFFT_PADDED_C2R) ? PFFTI_TRAFO_C2R_PADDED : PFFTI_TRAFO_C2R;

  return PX(local_size_partrafo)(
      rnk_n, n, ni, no, howmany, iblock, oblock,
      comm_cart, trafo_flag, pfft_flags,
      local_ni, local_i_start, local_no, local_o_start);
}

INT PX(local_size_many_r2r)(
    int rnk_n, const INT *n, const INT *ni, const INT *no,
    INT howmany, const INT *iblock, const INT *oblock,
    MPI_Comm comm_cart, unsigned pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  return PX(local_size_partrafo)(
      rnk_n, n, ni, no, howmany, iblock, oblock,
      comm_cart, PFFTI_TRAFO_R2R, pfft_flags,
      local_ni, local_i_start, local_no, local_o_start);
}


PX(plan) PX(plan_many_dft)(
    int rnk_n, const INT *n, const INT *ni, const INT *no,
    INT howmany, const INT *iblock, const INT *oblock,
    C *in, C *out, MPI_Comm comm_cart,
    int sign, unsigned pfft_flags
    )
{
  X(r2r_kind) *kinds=NULL;
  int *skip_trafos=NULL;

  return PX(plan_partrafo)(
      rnk_n, n, ni, no, howmany, iblock, oblock,
      (R*) in, (R*) out, comm_cart, sign, kinds, skip_trafos,
      PFFTI_TRAFO_C2C, pfft_flags);
}

PX(plan) PX(plan_many_dft_r2c)(
    int rnk_n, const INT *n, const INT *ni, const INT *no,
    INT howmany, const INT *iblock, const INT *oblock,
    R *in, C *out, MPI_Comm comm_cart,
    int sign, unsigned pfft_flags
    )
{
  unsigned trafo_flag = (pfft_flags & PFFT_PADDED_R2C) ? PFFTI_TRAFO_R2C_PADDED : PFFTI_TRAFO_R2C;
  X(r2r_kind) *kinds=NULL;
  int *skip_trafos=NULL;

  return PX(plan_partrafo)(
      rnk_n, n, ni, no, howmany, iblock, oblock,
      in, (R*) out, comm_cart, sign, kinds, skip_trafos,
      trafo_flag, pfft_flags);
}

PX(plan) PX(plan_many_dft_c2r)(
    int rnk_n, const INT *n, const INT *ni, const INT *no,
    INT howmany, const INT *iblock, const INT *oblock,
    C *in, R *out, MPI_Comm comm_cart,
    int sign, unsigned pfft_flags
    )
{
  unsigned trafo_flag = (pfft_flags & PFFT_PADDED_C2R) ? PFFTI_TRAFO_C2R_PADDED : PFFTI_TRAFO_C2R;
  X(r2r_kind) *kinds=NULL;
  int *skip_trafos=NULL;

  return PX(plan_partrafo)(
      rnk_n, n, ni, no, howmany, iblock, oblock,
      (R*) in, out, comm_cart, sign, kinds, skip_trafos,
      trafo_flag, pfft_flags);
}

PX(plan) PX(plan_many_r2r)(
    int rnk_n, const INT *n, const INT *ni, const INT *no,
    INT howmany, const INT *iblock, const INT *oblock,
    R *in, R *out, MPI_Comm comm_cart,
    const PX(r2r_kind) *kinds, unsigned pfft_flags
    )
{
  int sign=0;
  int *skip_trafos=NULL;

  return PX(plan_partrafo)(
      rnk_n, n, ni, no, howmany, iblock, oblock,
      in, out, comm_cart, sign, kinds, skip_trafos,
      PFFTI_TRAFO_R2R, pfft_flags);
}

PX(plan) PX(plan_many_dft_skipped)(
    int rnk_n, const INT *n, const INT *ni, const INT *no,
    INT howmany, const INT *iblock, const INT *oblock,
    const int *skip_trafos, C *in, C *out, MPI_Comm comm_cart,
    int sign, unsigned pfft_flags
    )
{
  X(r2r_kind) *kinds=NULL;

  return PX(plan_partrafo)(
      rnk_n, n, ni, no, howmany, iblock, oblock,
      (R*) in, (R*) out, comm_cart, sign, kinds, skip_trafos,
      PFFTI_TRAFO_C2C, pfft_flags);
}

PX(plan) PX(plan_many_dft_r2c_skipped)(
    int rnk_n, const INT *n, const INT *ni, const INT *no,
    INT howmany, const INT *iblock, const INT *oblock,
    const int *skip_trafos, R *in, C *out, MPI_Comm comm_cart,
    int sign, unsigned pfft_flags
    )
{
  unsigned trafo_flag = (pfft_flags & PFFT_PADDED_R2C) ? PFFTI_TRAFO_R2C_PADDED : PFFTI_TRAFO_R2C;
  X(r2r_kind) *kinds=NULL;

  return PX(plan_partrafo)(
      rnk_n, n, ni, no, howmany, iblock, oblock,
      in, (R*) out, comm_cart, sign, kinds, skip_trafos,
      trafo_flag, pfft_flags);
}

PX(plan) PX(plan_many_dft_c2r_skipped)(
    int rnk_n, const INT *n, const INT *ni, const INT *no,
    INT howmany, const INT *iblock, const INT *oblock,
    const int *skip_trafos, C *in, R *out, MPI_Comm comm_cart,
    int sign, unsigned pfft_flags
    )
{
  unsigned trafo_flag = (pfft_flags & PFFT_PADDED_C2R) ? PFFTI_TRAFO_C2R_PADDED : PFFTI_TRAFO_C2R;
  X(r2r_kind) *kinds=NULL;

  return PX(plan_partrafo)(
      rnk_n, n, ni, no, howmany, iblock, oblock,
      (R*) in, out, comm_cart, sign, kinds, skip_trafos,
      trafo_flag, pfft_flags);
}

PX(plan) PX(plan_many_r2r_skipped)(
    int rnk_n, const INT *n, const INT *ni, const INT *no,
    INT howmany, const INT *iblock, const INT *oblock,
    const int *skip_trafos, R *in, R *out, MPI_Comm comm_cart,
    const PX(r2r_kind) *kinds, unsigned pfft_flags
    )
{
  int sign=0;

  return PX(plan_partrafo)(
      rnk_n, n, ni, no, howmany, iblock, oblock,
      in, out, comm_cart, sign, kinds, skip_trafos,
      PFFTI_TRAFO_R2R, pfft_flags);
}


INT PX(local_size_many_gc)(
    int rnk_n, const INT *local_n, const INT *local_start,
    INT howmany, const INT *gc_below, const INT *gc_above,
    INT *local_ngc, INT *local_gc_start
    )
{
  return PX(local_size_gc_internal)(
      rnk_n, local_n, local_start, howmany, gc_below, gc_above,
      local_ngc, local_gc_start);
}


PX(gcplan) PX(plan_many_rgc)(
    int rnk_n, const INT *n,
    INT howmany, const INT *block,
    const INT *gc_below, const INT *gc_above,
    R *data, MPI_Comm comm,
    unsigned gc_flags
    )
{
  int rnk_pm;
  INT blk_for_remap[3];
  MPI_Comm *comms_pm;
  PX(gcplan) ths;
  MPI_Comm comm_cart = PX(assure_cart_comm)(comm);

  INT *pn = PX(malloc_INT)(rnk_n);
  for(int t=0; t<rnk_n; t++) 
    pn[t] = n[t];

  /* compute physical array size for padded real valued arrays */
  if(gc_flags & PFFT_GC_PADDED)
    pn[rnk_n-1] = 2*(n[rnk_n-1]/2+1);

  MPI_Cartdim_get(comm_cart, &rnk_pm);
  /* split comm cart into 1d comms */
  comms_pm = (MPI_Comm*) malloc(sizeof(MPI_Comm) * (size_t) rnk_pm);
  PX(split_cart_procmesh)(comm_cart, comms_pm);

  if( PX(needs_remap_nd)(rnk_n, comm_cart) ){
    /* 3d to 2d remap results in complicated blocks.
     * We ignore users input and use default block size. */
    PX(remap_nd_calculate_blocks)(rnk_n, pn, comm_cart,
        blk_for_remap);
    ths = PX(plan_rgc_internal)(
        rnk_n, pn, howmany, blk_for_remap, gc_below, gc_above,
        data, rnk_pm, comms_pm, comm_cart, gc_flags);
  } else
    ths = PX(plan_rgc_internal)(
        rnk_n, pn, howmany, block, gc_below, gc_above,
        data, rnk_pm, comms_pm, comm_cart, gc_flags);
  
  /* free one-dimensional comms */
  for(int t=0; t<rnk_pm; t++)
    MPI_Comm_free(comms_pm + t);
  free(comms_pm);  
  MPI_Comm_free(&comm_cart);
  free(pn);

  return ths;
}

PX(gcplan) PX(plan_many_cgc)(
    int rnk_n, const INT *n,
    INT howmany, const INT *block,
    const INT *gc_below, const INT *gc_above,
    C *data, MPI_Comm comm_cart,
    unsigned gc_flags
    )
{
  PX(gcplan) plan;
  INT *pn = PX(malloc_INT)(rnk_n);

  /* ignore padding for complex valued arrays */
  gc_flags &= ~PFFT_GC_PADDED;

  for(int t=0; t<rnk_n; t++)
    pn[t] = n[t];

  if(gc_flags & PFFT_GC_R2C){
    /* r2c: use c2c ghost cells with physical array size */
    pn[rnk_n-1] = n[rnk_n-1]/2+1;
    gc_flags &= ~PFFT_GC_R2C;
  }

  plan = PX(plan_many_rgc)(
      rnk_n, pn, 2*howmany, block, gc_below, gc_above,
      (R*) data, comm_cart, gc_flags);

  free(pn);
  return plan;
}


void PX(local_block_many_dft)(
    int rnk_n, const INT *ni, const INT *no,
    const INT *iblock, const INT *oblock,
    MPI_Comm comm_cart, int pid, unsigned pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  PX(local_block_partrafo)(
      rnk_n, ni, no, iblock, oblock,
      comm_cart, pid, PFFTI_TRAFO_C2C, pfft_flags,
      local_ni, local_i_start, local_no, local_o_start);
}

void PX(local_block_many_dft_r2c)(
    int rnk_n, const INT *ni, const INT *no,
    const INT *iblock, const INT *oblock,
    MPI_Comm comm_cart, int pid, unsigned pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  unsigned trafo_flag = (pfft_flags & PFFT_PADDED_R2C) ? PFFTI_TRAFO_R2C_PADDED : PFFTI_TRAFO_R2C;

  PX(local_block_partrafo)(
      rnk_n, ni, no, iblock, oblock,
      comm_cart, pid, trafo_flag, pfft_flags,
      local_ni, local_i_start, local_no, local_o_start);
}

void PX(local_block_many_dft_c2r)(
    int rnk_n, const INT *ni, const INT *no,
    const INT *iblock, const INT *oblock,
    MPI_Comm comm_cart, int pid, unsigned pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  unsigned trafo_flag = (pfft_flags & PFFT_PADDED_C2R) ? PFFTI_TRAFO_C2R_PADDED : PFFTI_TRAFO_C2R;

  PX(local_block_partrafo)(
      rnk_n, ni, no, iblock, oblock,
      comm_cart, pid, trafo_flag, pfft_flags,
      local_ni, local_i_start, local_no, local_o_start);
}

void PX(local_block_many_r2r)(
    int rnk_n, const INT *ni, const INT *no,
    const INT *iblock, const INT *oblock,
    MPI_Comm comm_cart, int pid, unsigned pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  PX(local_block_partrafo)(
      rnk_n, ni, no, iblock, oblock,
      comm_cart, pid, PFFTI_TRAFO_R2R, pfft_flags,
      local_ni, local_i_start, local_no, local_o_start);
}

