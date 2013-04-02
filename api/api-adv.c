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

static void calculate_3dto2d_blocks(
    const INT *n, MPI_Comm comm_cart_3d,
    INT *iblk);

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
  return PX(local_size_partrafo)(
    rnk_n, n, ni, no, howmany, iblock, oblock,
    comm_cart, PFFTI_TRAFO_R2C, pfft_flags,
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
  return PX(local_size_partrafo)(
    rnk_n, n, ni, no, howmany, iblock, oblock,
    comm_cart, PFFTI_TRAFO_C2R, pfft_flags,
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
  X(r2r_kind) *kinds=NULL;
  int *skip_trafos=NULL;

  return PX(plan_partrafo)(
    rnk_n, n, ni, no, howmany, iblock, oblock,
    in, (R*) out, comm_cart, sign, kinds, skip_trafos,
    PFFTI_TRAFO_R2C, pfft_flags);
}

PX(plan) PX(plan_many_dft_c2r)(
    int rnk_n, const INT *n, const INT *ni, const INT *no,
    INT howmany, const INT *iblock, const INT *oblock,
    C *in, R *out, MPI_Comm comm_cart,
    int sign, unsigned pfft_flags
    )
{
  X(r2r_kind) *kinds=NULL;
  int *skip_trafos=NULL;

  return PX(plan_partrafo)(
    rnk_n, n, ni, no, howmany, iblock, oblock,
    (R*) in, out, comm_cart, sign, kinds, skip_trafos,
    PFFTI_TRAFO_C2R, pfft_flags);
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

INT PX(local_size_many_gc)(
    int rnk_n, const INT *local_n, const INT *local_start, INT alloc_local,
    INT howmany, const INT *gc_below, const INT *gc_above,
    INT *local_ngc, INT *local_gc_start
    )
{
  INT alloc_local_gc;

  alloc_local_gc = PX(local_size_gc_internal)(
      rnk_n, local_n, local_start, howmany, gc_below, gc_above,
      local_ngc, local_gc_start);

  return MAX(alloc_local, alloc_local_gc);
}


PX(gcplan) PX(plan_many_rgc)(
    int rnk_n, const INT *n,
    INT howmany, const INT *block,
    const INT *gc_below, const INT *gc_above,
    R *data, MPI_Comm comm_cart,
    unsigned gc_flags
    )
{
  int rnk_pm;
  INT blk_3dto2d[3];
  MPI_Comm *comms_pm;
  PX(gcplan) ths;

  MPI_Cartdim_get(comm_cart, &rnk_pm);
  /* split comm cart into 1d comms */
  comms_pm = (MPI_Comm*) malloc(sizeof(MPI_Comm) * (size_t) rnk_pm);
  PX(split_cart_procmesh)(comm_cart, comms_pm);

  if( PX(needs_3dto2d_remap)(rnk_n, comm_cart) ){
    /* 3d to 2d remap results in complicated blocks.
     * We ignore users input and use default block size. */
    calculate_3dto2d_blocks(n, comm_cart,
        blk_3dto2d);
    ths = PX(plan_rgc_internal)(
        rnk_n, n, howmany, blk_3dto2d, gc_below, gc_above,
        data, rnk_pm, comms_pm, comm_cart, gc_flags);
  } else
    ths = PX(plan_rgc_internal)(
        rnk_n, n, howmany, block, gc_below, gc_above,
        data, rnk_pm, comms_pm, comm_cart, gc_flags);
  
  /* free one-dimensional comms */
  for(int t=0; t<rnk_pm; t++)
    MPI_Comm_free(comms_pm + t);
  free(comms_pm);  

  return ths;
}

static void calculate_3dto2d_blocks(
    const INT *n, MPI_Comm comm_cart_3d,
    INT *iblk
    )
{
  int q0, q1, p0, p1;
  INT mblk[3], oblk[3];

  PX(get_procmesh_dims_2d)(comm_cart_3d, &p0, &p1, &q0, &q1);
  PX(default_block_size_3dto2d)(n, p0, p1, q0, q1,
      iblk, mblk, oblk);
}



PX(gcplan) PX(plan_many_cgc)(
    int rnk_n, const INT *n,
    INT howmany, const INT *block,
    const INT *gc_below, const INT *gc_above,
    C *data, MPI_Comm comm_cart,
    unsigned gc_flags
    )
{
  return PX(plan_many_rgc)(
      rnk_n, n, 2*howmany, block, gc_below, gc_above,
      (R*) data, comm_cart, gc_flags);
}
