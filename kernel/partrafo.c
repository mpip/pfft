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

static PX(plan) mkplan(
    int rnk_n, int rnk_pm);

static void malloc_and_split_cart_procmesh(
    int rnk_n, MPI_Comm comm_cart, 
    int *rnk_pm, MPI_Comm **comms_pm);

static void init_param_local_size(
    INT *lni, INT *lis, INT *dummy_ln, INT *dummy_ls, INT *lno, INT *los,
    unsigned transp_flag,
    INT **lni_to, INT **lis_to, INT **lno_to, INT **los_to,
    INT **lni_ti, INT **lis_ti, INT **lno_ti, INT **los_ti);
static void save_param_into_plan(
    int rnk_n, const INT *n, const INT *ni, const INT *no,
    INT howmany, const INT *iblock, const INT *mblock, const INT *oblock,
    MPI_Comm comm_cart, int rnk_pm, MPI_Comm *comms_pm,
    R *in, R *out, int sign, const X(r2r_kind) *kinds,
    unsigned transp_flag, unsigned trafo_flag,
    unsigned opt_flag, unsigned fftw_flags,
    unsigned pfft_flags, 
    PX(plan) ths);
static void init_param_size_and_trafo_flags(
    int rnk_n, const INT *n, const INT *ni, const INT *no, int rnk_pm,
    unsigned trafo_flag, unsigned transp_flag, const int *skip_trafos,
    INT *n_to, INT *ni_to, INT *no_to, unsigned *trafo_flags_to,
    INT *n_ti, INT *ni_ti, INT *no_ti, unsigned *trafo_flags_ti);
static unsigned get_skip_flag(
    const int *skip_trafos, int t);
static void set_plans_to_null(
    int rnk_pm, unsigned transp_flag,
    outrafo_plan *trafos, gtransp_plan *remaps);

static void evaluate_blocks(
    int rnk_n, const INT *ni, const INT *no,
    const INT *iblock_user, const INT *oblock_user,
    int rnk_pm, const MPI_Comm *comms_pm,
    unsigned trafo_flag, unsigned transp_flag,
    INT *iblk, INT *mblk, INT *oblk);

static unsigned extract_transp_flag(
    unsigned pfft_flags);
static unsigned extract_opt_flag(
    unsigned pfft_flags);
static unsigned extract_io_flag(
    unsigned pfft_flags);
static unsigned extract_fftw_flags(
    unsigned pfft_flags);

/* TRANPOSED_OUT:
 * - iblock gives block size of nontransp. input
 * - oblock gives block size of transp. output
 * TRANPOSED_IN:
 * - iblock gives block size of transp. input
 * - oblock gives block size of nontransp. output
 * NOT TRANSPOSED:
 * - iblock gives block size of nontransp. input
 * - defblock gives block size of intermediate transp. data layout
 * - oblock gives block size of nontransp. output
 */


INT PX(local_size_partrafo)(
    int rnk_n, const INT *n, const INT *ni, const INT *no,
    INT howmany, const INT *iblock_user, const INT *oblock_user,
    MPI_Comm comm_cart,
    unsigned trafo_flag_user, unsigned pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  unsigned transp_flag = extract_transp_flag(pfft_flags);
  unsigned *trafo_flags_to, *trafo_flags_ti;
  INT mem=1, mem_tmp;
  INT *ni_to, *n_to, *no_to, *ni_ti, *n_ti, *no_ti;
  INT *iblk, *mblk, *oblk;
  INT *dummy_ln, *dummy_ls;
  INT *lni_to, *lis_to, *lno_to, *los_to; 
  INT *lni_ti, *lis_ti, *lno_ti, *los_ti;
  int rnk_pm;
  MPI_Comm *comms_pm;

  malloc_and_split_cart_procmesh(rnk_n, comm_cart,
      &rnk_pm, &comms_pm);

  ni_to = PX(malloc_INT)(rnk_n);
  n_to  = PX(malloc_INT)(rnk_n);
  no_to = PX(malloc_INT)(rnk_n);
  ni_ti = PX(malloc_INT)(rnk_n);
  n_ti  = PX(malloc_INT)(rnk_n);
  no_ti = PX(malloc_INT)(rnk_n);
  
  trafo_flags_to = PX(malloc_unsigned)(rnk_pm + 1);
  trafo_flags_ti = PX(malloc_unsigned)(rnk_pm + 1);

  dummy_ln = PX(malloc_INT)(rnk_n);
  dummy_ls = PX(malloc_INT)(rnk_n);

  iblk = PX(malloc_INT)(rnk_pm);
  mblk = PX(malloc_INT)(rnk_pm);
  oblk = PX(malloc_INT)(rnk_pm);
  
  /* calculate blocksizes according to trafo and transp flags */
  evaluate_blocks(rnk_n, ni, no, iblock_user, oblock_user,
      rnk_pm, comms_pm, trafo_flag_user, transp_flag,
      iblk, mblk, oblk);

  /* split trafo into transposed out and transposed in step */
  init_param_size_and_trafo_flags(
      rnk_n, n, ni, no, rnk_pm, trafo_flag_user, transp_flag, NULL,
      n_to, ni_to, no_to, trafo_flags_to,
      n_ti, ni_ti, no_ti, trafo_flags_ti);
 
  init_param_local_size(
      local_ni, local_i_start, dummy_ln, dummy_ls,
      local_no, local_o_start, transp_flag,
      &lni_to, &lis_to, &lno_to, &los_to,
      &lni_ti, &lis_ti, &lno_ti, &los_ti);

  /* plan with transposed output */
  if( ~transp_flag & PFFT_TRANSPOSED_IN ){
    mem_tmp = PX(local_size_partrafo_transposed)(
        rnk_n, n_to, ni_to, no_to, howmany, iblk, mblk,
        rnk_pm, comms_pm, PFFT_TRANSPOSED_OUT, trafo_flags_to, 
        lni_to, lis_to, lno_to, los_to);
    mem = MAX(mem, mem_tmp);

    mem_tmp = PX(local_size_remap_3dto2d_transposed)(
        rnk_n, ni_to, howmany, comm_cart, PFFT_TRANSPOSED_OUT, trafo_flags_to[rnk_pm],
        lni_to, lis_to, dummy_ln, dummy_ls);
    mem = MAX(mem, mem_tmp);
  }

  /* plan with transposed input */
  if( ~transp_flag & PFFT_TRANSPOSED_OUT ){
    mem_tmp = PX(local_size_partrafo_transposed)(
        rnk_n, n_ti, ni_ti, no_ti, howmany, mblk, oblk,
        rnk_pm, comms_pm, PFFT_TRANSPOSED_IN, trafo_flags_ti, 
        lni_ti, lis_ti, lno_ti, los_ti);
    mem = MAX(mem, mem_tmp);

    mem_tmp = PX(local_size_remap_3dto2d_transposed)(
        rnk_n, no_ti, howmany, comm_cart, PFFT_TRANSPOSED_IN, trafo_flags_ti[rnk_pm],
        dummy_ln, dummy_ls, lno_ti, los_ti);
    mem = MAX(mem, mem_tmp);
  }

  /* free one-dimensional comms */
  for(int t=0; t<rnk_pm; t++)
    MPI_Comm_free(&comms_pm[t]);
  free(comms_pm);

  free(iblk); free(mblk); free(oblk);
  free(dummy_ln); free(dummy_ls);
  free(trafo_flags_to); free(trafo_flags_ti);
  free(ni_to); free(n_to); free(no_to);
  free(ni_ti); free(n_ti); free(no_ti);

  /* local size of r2c input in units of real */
  if(trafo_flag_user & PFFTI_TRAFO_R2C){
    local_ni[rnk_n-1] *= 2;
    local_i_start[rnk_n-1] *= 2;
  }
    
  /* local size of c2r output in units of real */
  if(trafo_flag_user & PFFTI_TRAFO_C2R){
    local_no[rnk_n-1] *= 2;
    local_o_start[rnk_n-1] *= 2;
  }
  
  return mem;
}


PX(plan) PX(plan_partrafo)(
    int rnk_n, const INT *n, const INT *ni, const INT *no,
    INT howmany, const INT *iblock_user, const INT *oblock_user,
    R *in, R *out, MPI_Comm comm_cart,
    int sign, const X(r2r_kind) *kinds, const int *skip_trafos_user,
    unsigned trafo_flag, unsigned pfft_flags
    )
{
  unsigned transp_flag = extract_transp_flag(pfft_flags);
  unsigned opt_flag = extract_opt_flag(pfft_flags);
  unsigned io_flag = extract_io_flag(pfft_flags);
  unsigned fftw_flags = extract_fftw_flags(pfft_flags);
  unsigned *trafo_flags_to, *trafo_flags_ti;
  INT *ni_to, *n_to, *no_to, *ni_ti, *n_ti, *no_ti;
  INT *iblk, *mblk, *oblk;
  PX(plan) ths;
  int rnk_pm;
  MPI_Comm *comms_pm;

  /* transposed input and output in one plan is not supported */
  if(transp_flag & PFFT_TRANSPOSED_OUT)
    if(transp_flag & PFFT_TRANSPOSED_IN)
      return NULL;

  /* transposed in R2C not supported */
  if(trafo_flag & PFFTI_TRAFO_R2C)
    if(transp_flag & PFFT_TRANSPOSED_IN)
      return NULL;

  /* transposed out C2R not supported */
  if(trafo_flag & PFFTI_TRAFO_C2R)
    if(transp_flag & PFFT_TRANSPOSED_OUT)
      return NULL;

  MPI_Cartdim_get(comm_cart, &rnk_pm);

  /* dimension of FFT is not allowed to be smaller than procmesh dimension */
  if(rnk_n < rnk_pm)
    return NULL;

  /* equal dimension of FFT and procmesh only implemented for 3 dimensions */
  if(rnk_n == rnk_pm)
    if(rnk_n != 3)
      return NULL;
  
  malloc_and_split_cart_procmesh(rnk_n, comm_cart,
      &rnk_pm, &comms_pm);

  ths = mkplan(rnk_n, rnk_pm);

  ni_to = PX(malloc_INT)(rnk_n);
  n_to  = PX(malloc_INT)(rnk_n);
  no_to = PX(malloc_INT)(rnk_n);
  ni_ti = PX(malloc_INT)(rnk_n);
  n_ti  = PX(malloc_INT)(rnk_n);
  no_ti = PX(malloc_INT)(rnk_n);
  
  trafo_flags_to = PX(malloc_unsigned)(rnk_pm + 1);
  trafo_flags_ti = PX(malloc_unsigned)(rnk_pm + 1);

  iblk = PX(malloc_INT)(rnk_pm);
  mblk = PX(malloc_INT)(rnk_pm);
  oblk = PX(malloc_INT)(rnk_pm);
  
  /* calculate blocksizes according to trafo and transp flags */
  evaluate_blocks(rnk_n, ni, no, iblock_user, oblock_user,
      rnk_pm, comms_pm, trafo_flag, transp_flag,
      iblk, mblk, oblk);
  
  /* Avoid recalculation of the same parameters all the time. */
  save_param_into_plan(rnk_n, n, ni, no, howmany, iblk, mblk, oblk,
      comm_cart, rnk_pm, comms_pm, in, out, sign, kinds,
      transp_flag, trafo_flag, opt_flag, fftw_flags, pfft_flags,
      ths);

  /* split trafo into transposed out and transposed in step */
  init_param_size_and_trafo_flags(
      rnk_n, n, ni, no, rnk_pm, trafo_flag, transp_flag, skip_trafos_user,
      n_to, ni_to, no_to, trafo_flags_to,
      n_ti, ni_ti, no_ti, trafo_flags_ti);

  /* For C2R trafos the output of the forward (transpose) step ends
   * up in pointer 'in', since we skip the last local transposition.
   * For all other trafos the input of the backward step is given by
   * pointer 'out', since we ommit its first local transposition.
   * So, for forward and backward steps we use 'in' for input and
   * 'out' for output. */

  /* plan with transposed output */
  if(transp_flag & PFFT_TRANSPOSED_IN){
    ths->remap_3dto2d[0] = NULL;
    set_plans_to_null(rnk_pm, PFFT_TRANSPOSED_OUT,
        ths->serial_trafo, ths->global_remap);
  } else {
    ths->remap_3dto2d[0] = PX(plan_remap_3dto2d_transposed)(
        rnk_n, ni_to, howmany, comm_cart, in, out,
        PFFT_TRANSPOSED_OUT, trafo_flags_to[rnk_pm], opt_flag, io_flag, fftw_flags);

    /* If remap_3dto2d exists, go on with in-place transforms in order to preserve input. */
    if( (ths->remap_3dto2d[0] != NULL) && (~io_flag & PFFT_DESTROY_INPUT) )
      in = out;

    PX(plan_partrafo_transposed)(rnk_n, n_to, ni_to, no_to,
        howmany, iblk, mblk,
        rnk_pm, comms_pm, in, out, sign, kinds,
        PFFT_TRANSPOSED_OUT, trafo_flags_to, opt_flag, io_flag, fftw_flags, 
        ths->serial_trafo, ths->global_remap);

    /* Go on with in-place transforms in order to preserve input. */
    if( ~io_flag & PFFT_DESTROY_INPUT)
      in = out;
  }

  /* plan with transposed input */
  if(transp_flag & PFFT_TRANSPOSED_OUT){
    set_plans_to_null(rnk_pm, PFFT_TRANSPOSED_IN,
        ths->serial_trafo, ths->global_remap);
    ths->remap_3dto2d[1] = NULL;
  } else {
    PX(plan_partrafo_transposed)(rnk_n, n_ti, ni_ti, no_ti,
        howmany, mblk, oblk,
        rnk_pm, comms_pm, in, out, sign, kinds,
        PFFT_TRANSPOSED_IN, trafo_flags_ti, opt_flag, io_flag, fftw_flags, 
        ths->serial_trafo, ths->global_remap);

    /* Go on with in-place transforms in order to preserve input. */
    if( ~io_flag & PFFT_DESTROY_INPUT )
      in = out;

    ths->remap_3dto2d[1] = PX(plan_remap_3dto2d_transposed)(
        rnk_n, no_ti, howmany, comm_cart, out, in,
        PFFT_TRANSPOSED_IN, trafo_flags_ti[rnk_pm], opt_flag, io_flag, fftw_flags);
  }

  /* free one-dimensional comms */
  for(int t=0; t<rnk_pm; t++)
    MPI_Comm_free(comms_pm + t);
  free(comms_pm);  

  free(iblk); free(mblk); free(oblk);
  free(trafo_flags_to); free(trafo_flags_ti);
  free(ni_to); free(n_to); free(no_to);
  free(ni_ti); free(n_ti); free(no_ti);
  return ths;
}

/* Split comm cart into 1d comms.
 * Special case: 3d data with 3d decomposition
 * => First remap to 2d decomposition. */
static void malloc_and_split_cart_procmesh(
    int rnk_n, MPI_Comm comm_cart, 
    int *rnk_pm, MPI_Comm **comms_pm
    )
{
  MPI_Cartdim_get(comm_cart, rnk_pm);

  if( PX(needs_3dto2d_remap)(rnk_n, comm_cart) )
    *rnk_pm = 2;

  *comms_pm = (MPI_Comm*) malloc(sizeof(MPI_Comm) * (size_t) *rnk_pm);

  if( PX(needs_3dto2d_remap)(rnk_n, comm_cart) ){
    PX(split_cart_procmesh_3dto2d_p0q0)(comm_cart,
        *comms_pm + 0);
    PX(split_cart_procmesh_3dto2d_p1q1)(comm_cart,
        *comms_pm + 1);
  } else
    PX(split_cart_procmesh)(comm_cart, *comms_pm);
}






static void save_param_into_plan(
    int rnk_n, const INT *n, const INT *ni, const INT *no,
    INT howmany, const INT *iblock, const INT *mblock, const INT *oblock,
    MPI_Comm comm_cart, int rnk_pm, MPI_Comm *comms_pm,
    R *in, R *out, int sign, const X(r2r_kind) *kinds,
    unsigned transp_flag, unsigned trafo_flag,
    unsigned opt_flag, unsigned fftw_flags, 
    unsigned pfft_flags,
    PX(plan) ths
    )
{
  ths->rnk_n = rnk_n;
  for(int t=0; t<rnk_n; t++){
    ths->n[t]  = n[t];
    ths->ni[t] = ni[t];
    ths->no[t] = no[t];
  }
  ths->howmany = howmany;
  for(int t=0; t<rnk_pm; t++){
    ths->iblock[t] = iblock[t];
    ths->mblock[t] = mblock[t];
    ths->oblock[t] = oblock[t];
  }

  MPI_Comm_dup(comm_cart, &ths->comm_cart);
  
  ths->rnk_pm = rnk_pm;
  for(int t=0; t<rnk_pm; t++){
    MPI_Comm_dup(comms_pm[t], &ths->comms_pm[t]);
    MPI_Comm_size(ths->comms_pm[t], &ths->np[t]);
  }

  ths->in = in;
  ths->out = out;
  ths->sign = sign;
  if(kinds != NULL){
    ths->kinds = (X(r2r_kind)*) malloc(sizeof(X(r2r_kind)) * (size_t) rnk_n);
    for(int t=0; t<rnk_n; t++)
      ths->kinds[t] = kinds[t];
  } else
    ths->kinds = NULL;
  ths->fftw_flags = fftw_flags;
  ths->transp_flag = transp_flag;
  ths->trafo_flag = trafo_flag;
  ths->opt_flag = opt_flag;
  ths->pfft_flags = pfft_flags;
}


static void evaluate_blocks(
    int rnk_n, const INT *ni, const INT *no,
    const INT *iblock_user, const INT *oblock_user,
    int rnk_pm, const MPI_Comm *comms_pm,
    unsigned trafo_flag, unsigned transp_flag,
    INT *iblk, INT *mblk, INT *oblk
    )
{
  INT *pni, *pnm, *pno;
  const INT *mblock_user;

  pni = PX(malloc_INT)(rnk_n);
  pno = PX(malloc_INT)(rnk_n);

  PX(physical_dft_size)(rnk_n, ni, trafo_flag,
      pni);
  PX(physical_dft_size)(rnk_n, no, trafo_flag,
      pno);

  if( (trafo_flag & PFFTI_TRAFO_C2R) || (transp_flag & PFFT_TRANSPOSED_IN) )
    pnm = pni;
  else
    pnm = pno;

  /* Default: user gives block sizes of nontransposed input and output.
   * We choose default block size for intermediate transposed layout */
  mblock_user = PFFT_DEFAULT_BLOCKS;

  if( transp_flag & PFFT_TRANSPOSED_IN )
    mblock_user = iblock_user;

  if( transp_flag & PFFT_TRANSPOSED_OUT )
    mblock_user = oblock_user;

  /* initialize all blocks to zero */
  for(int t=0; t<rnk_pm; t++)
    iblk[t] = mblk[t] = oblk[t] = 0;

  /* not needed for transposed input */
  if( ~transp_flag & PFFT_TRANSPOSED_IN )
    PX(evaluate_user_block_size)(rnk_pm, pni, iblock_user, comms_pm, iblk);

  PX(evaluate_user_block_size)(rnk_pm, pnm+1, mblock_user, comms_pm, mblk);

  /* not needed for transposed output */
  if( ~transp_flag & PFFT_TRANSPOSED_OUT )
    PX(evaluate_user_block_size)(rnk_pm, pno, oblock_user, comms_pm, oblk);

  free(pni); free(pno);
}

static void init_param_local_size(
    INT *lni, INT *lis, INT *dummy_ln, INT *dummy_ls, INT *lno, INT *los,
    unsigned transp_flag,
    INT **lni_to, INT **lis_to, INT **lno_to, INT **los_to,
    INT **lni_ti, INT **lis_ti, INT **lno_ti, INT **los_ti
    )
{
  /* default does trafo at forward (transp. out) step */
  /* c2r does trafo at backward (transp. in) step */
  *lni_to = lni; *lis_to = lis;
  *lno_to = dummy_ln; *los_to = dummy_ls;
  *lni_ti = dummy_ln; *lis_ti = dummy_ls;
  *lno_ti = lno; *los_ti = los;

  /* TRANSPOSED_OUT */
  if(transp_flag & PFFT_TRANSPOSED_OUT){
    *lni_to = lni; *lis_to = lis;
    *lno_to = lno; *los_to = los;
    *lni_ti = dummy_ln; *lis_ti = dummy_ls;
    *lno_ti = dummy_ln; *los_ti = dummy_ls;
  }

  /* TRANSPOSED_IN */
  if(transp_flag & PFFT_TRANSPOSED_IN){
    *lni_to = dummy_ln; *lis_to = dummy_ls;
    *lno_to = dummy_ln; *los_to = dummy_ls;
    *lni_ti = lni; *lis_ti = lis;
    *lno_ti = lno; *los_ti = los;
  }
}


static void init_param_size_and_trafo_flags(
    int rnk_n, const INT *n, const INT *ni, const INT *no, int rnk_pm,
    unsigned trafo_flag, unsigned transp_flag, const int *skip_trafos,
    INT *n_to, INT *ni_to, INT *no_to,
    unsigned *trafo_flags_to,
    INT *n_ti, INT *ni_ti, INT *no_ti,
    unsigned *trafo_flags_ti
    )
{
  int t;

  /* Default: Compute trafo at forward (transp. out) step and perform a backward (transp. in) transposition. */
  for(t=0; t<rnk_pm+1; t++)
    trafo_flags_to[t] = trafo_flag | get_skip_flag(skip_trafos, t);
  PX(vcopy_INT)(rnk_n, ni, ni_to);
  PX(vcopy_INT)(rnk_n, n,  n_to);
  PX(vcopy_INT)(rnk_n, no, no_to);

  for(t=0; t<rnk_pm+1; t++)
    trafo_flags_ti[t] = trafo_flag | PFFTI_TRAFO_SKIP;
  PX(vcopy_INT)(rnk_n, no, ni_ti);
  PX(vcopy_INT)(rnk_n, no, n_ti);
  PX(vcopy_INT)(rnk_n, no, no_ti);

  /* TRANSPOSED_IN: Compute trafo on backward step */
  if(transp_flag & PFFT_TRANSPOSED_IN){
    for(t=0; t<rnk_pm+1; t++)
      trafo_flags_ti[t] = trafo_flag | get_skip_flag(skip_trafos, t);
    PX(vcopy_INT)(rnk_n, ni, ni_ti);
    PX(vcopy_INT)(rnk_n, n,  n_ti);
    PX(vcopy_INT)(rnk_n, no, no_ti);
  }
    
  /* R2C: Compute trafo always at forward (transp. out) step.
   * Therefore it doesn't work with TRANSPOSED_IN. */
  if( trafo_flag & PFFTI_TRAFO_R2C ){
    /* Compute trafo on forward step */
    for(t=0; t<rnk_pm+1; t++)
      trafo_flags_to[t] = ((t==rnk_pm) ? PFFTI_TRAFO_R2C : PFFTI_TRAFO_C2C) | get_skip_flag(skip_trafos, t);
    PX(vcopy_INT)(rnk_n, ni, ni_to);
    PX(vcopy_INT)(rnk_n, n,  n_to);
    PX(vcopy_INT)(rnk_n, no, no_to);
    
    /* Do backward transpose on complex array */
    for(t=0; t<rnk_pm+1; t++)
      trafo_flags_ti[t] = PFFTI_TRAFO_C2C | PFFTI_TRAFO_SKIP;
    for(t=0; t<rnk_n-1; t++)
      ni_ti[t] = n_ti[t] = no_ti[t] = no[t];
    ni_ti[t] = n_ti[t] = no_ti[t] = no[t]/2+1;
  }
  
  /* C2R: Compute trafo always at backward (transp. in) step.
   * Therefore it doesn't work with TRANSPOSED_OUT. */
  if( trafo_flag & PFFTI_TRAFO_C2R ){
    /* Do forward transpose on complex array */
    for(t=0; t<rnk_pm+1; t++)
      trafo_flags_to[t] = PFFTI_TRAFO_C2C | PFFTI_TRAFO_SKIP;
    for(t=0; t<rnk_n-1; t++)
      ni_to[t] = n_to[t] = no_to[t] = ni[t];
    ni_to[t] = n_to[t] = no_to[t] = ni[t]/2+1;

    /* Compute trafo on backward step */
    for(t=0; t<rnk_pm+1; t++)
      trafo_flags_ti[t] = ((t==rnk_pm) ? PFFTI_TRAFO_C2R : PFFTI_TRAFO_C2C) | get_skip_flag(skip_trafos, t);
    PX(vcopy_INT)(rnk_n, ni, ni_ti);
    PX(vcopy_INT)(rnk_n, n,  n_ti);
    PX(vcopy_INT)(rnk_n, no, no_ti);
  }

  /* TRANSPOSED_OUT: Ignore the '*_ti' parameters */
  if(transp_flag & PFFT_TRANSPOSED_OUT){
    for(t=0; t<rnk_pm+1; t++)
      trafo_flags_ti[t] = 0;
    for(t=0; t<rnk_n; t++)
      ni_ti[t] = n_ti[t] = no_ti[t] = 0;
  }
  
  /* TRANSPOSED_IN: Ignore the '*_to' parameters */
  if(transp_flag & PFFT_TRANSPOSED_IN){
    for(t=0; t<rnk_pm+1; t++)
      trafo_flags_to = 0;
    for(t=0; t<rnk_n; t++)
      ni_to[t] = n_to[t] = no_to[t] = 0;
  }

  /* If none of the flags TRANSPOSED_IN and TRANSPOSED_OUT is set, 
   * the results of out-of-place plans would be stored in 'in'.
   * However, we must ensure that output of out-of-place plans ends up in 'out'.
   * Therefore, we skip the (needless) first local transform of the TRANSPOSED_IN step. 
   * C2R case:
   * Here, the trafo is computed during the TRANSPOSED_IN step, so we skip
   * the last local transform of the TRANSPOSED_OUT step.
   * */
  if(~transp_flag & (PFFT_TRANSPOSED_IN || PFFT_TRANSPOSED_OUT) ){
    if(trafo_flag & PFFTI_TRAFO_C2R)
      trafo_flags_to[0] = PFFTI_TRAFO_PHANTOM;
    else
      trafo_flags_ti[0] = PFFTI_TRAFO_PHANTOM;
  }
}

static unsigned get_skip_flag(
    const int *skip_trafos, int t
    )
{
  if(skip_trafos == NULL)
    return 0;

  return (skip_trafos[t]) ? PFFTI_TRAFO_SKIP : 0;
}

static void set_plans_to_null(
    int rnk_pm, unsigned transp_flag,
    outrafo_plan *trafos, gtransp_plan *remaps
    )
{
  int s;

  for(int t=0; t<rnk_pm+1; t++){
    s = (transp_flag & PFFT_TRANSPOSED_IN) ? rnk_pm+1+t : t;
    trafos[s] = NULL;
  }
  
  for(int t=0; t<rnk_pm; t++){
    s = (transp_flag & PFFT_TRANSPOSED_IN) ? rnk_pm+t : t;
    remaps[s] = NULL;
  }
}


static PX(plan) mkplan(
    int rnk_n, int rnk_pm
    )
{
  PX(plan) ths = (plan_s*) malloc(sizeof(plan_s));

  ths->n  = PX(malloc_INT)(rnk_n);
  ths->ni = PX(malloc_INT)(rnk_n);
  ths->no = PX(malloc_INT)(rnk_n);
  ths->iblock = PX(malloc_INT)(rnk_pm);
  ths->mblock = PX(malloc_INT)(rnk_pm);
  ths->oblock = PX(malloc_INT)(rnk_pm);

  ths->comms_pm = (MPI_Comm*) malloc(sizeof(MPI_Comm) * (size_t) rnk_pm);
  ths->np = PX(malloc_int)(rnk_pm);
  ths->kinds = NULL; /* allocate later if needed */

  /* allocate array of plans */
  ths->serial_trafo = (outrafo_plan *)
      malloc(sizeof(outrafo_plan) * (size_t) (2*rnk_pm+2));
  ths->global_remap = (gtransp_plan *)
      malloc(sizeof(gtransp_plan) * (size_t) (2*rnk_pm));

  /* allocate timer and set all times to zero */
  ths->timer = PX(mktimer)(rnk_pm);

  return ths;
}


void PX(rmplan)(
    PX(plan) ths
    )
{
  /* plan was already destroyed or never initialized */
  if(ths==NULL)
    return;

  for(int t=0; t<2*ths->rnk_pm+2; t++)
    PX(outrafo_rmplan)(ths->serial_trafo[t]);
  free(ths->serial_trafo);

  for(int t=0; t<2*ths->rnk_pm; t++)
    PX(gtransp_rmplan)(ths->global_remap[t]);
  free(ths->global_remap);

  free(ths->ni); free(ths->n); free(ths->no);
  free(ths->iblock); free(ths->mblock); free(ths->oblock);

  MPI_Comm_free(&ths->comm_cart);
  for(int t=0; t<ths->rnk_pm; t++)
    MPI_Comm_free(&ths->comms_pm[t]);
  free(ths->np);
  free(ths->comms_pm);
  if(ths->kinds != NULL)
    free(ths->kinds);

  for(int t=0; t<2; t++)
    PX(remap_3dto2d_rmplan)(ths->remap_3dto2d[t]);

  PX(destroy_timer)(ths->timer);

  /* free memory */
  free(ths);
  /* ths=NULL; would be senseless, since we can not change the pointer itself */
}

static unsigned extract_transp_flag(
    unsigned pfft_flags
    )
{
  unsigned transp_flag = PFFT_TRANSPOSED_NONE;

  if(pfft_flags & PFFT_TRANSPOSED_IN)
    transp_flag |= PFFT_TRANSPOSED_IN;
  if(pfft_flags & PFFT_TRANSPOSED_OUT)
    transp_flag |= PFFT_TRANSPOSED_OUT;

  return transp_flag;
}

static unsigned extract_opt_flag(
    unsigned pfft_flags
    )
{
  return (pfft_flags & PFFT_NO_TUNE) ? PFFT_NO_TUNE : 0;
}

static unsigned extract_io_flag(
    unsigned pfft_flags
    )
{
  unsigned flag = 0;

  if( (pfft_flags & PFFT_DESTROY_INPUT) && (~pfft_flags & PFFT_PRESERVE_INPUT))
    flag |= PFFT_DESTROY_INPUT;
  else
    flag |= PFFT_PRESERVE_INPUT;

  return flag;
}

/* Assure that only PFFT-compatible FFTW flags are used. */
static unsigned extract_fftw_flags(
    unsigned pfft_flags
    )
{
  unsigned fftw_flags = FFTW_MEASURE;

  if(pfft_flags & PFFT_ESTIMATE)
    fftw_flags = FFTW_ESTIMATE;
  if(pfft_flags & PFFT_PATIENT)
    fftw_flags = FFTW_PATIENT;
  if(pfft_flags & PFFT_EXHAUSTIVE)
    fftw_flags = FFTW_EXHAUSTIVE;

  /* PFFT_PRESERVE_INPUT needs one out-of-place FFT with
   * FFTW_PRESERVE_INPUT followed by several in-place FFTs */
  if( (pfft_flags & PFFT_DESTROY_INPUT) && (~pfft_flags & PFFT_PRESERVE_INPUT))
    fftw_flags |= FFTW_DESTROY_INPUT;
  else
    fftw_flags |= FFTW_PRESERVE_INPUT;

  return fftw_flags;
}

