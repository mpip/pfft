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


static INT initialize_Nb(
    int rnk_pm, const INT *local_ni, const INT *local_no, unsigned transp_flag);
static INT calculate_tuple_size(
    int rnk_n, const INT *ni, const INT *no, int rnk_pm, INT howmany,
    unsigned trafo_flag, unsigned transp_flag);
static void local_size_transposed(
    int rnk_n, const INT *ni, const INT *no,
    const INT *iblock, const INT *oblock,
    int rnk_pm, const MPI_Comm *comms_pm,
    unsigned trafo_flag, unsigned transp_flag,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start);
static void decompose_nontransposed(
    int rnk_n, const INT *n, const INT *blk,
    int rnk_pm, const MPI_Comm *comms_pm,
    unsigned trafo_flag,
    INT *local_n, INT *local_start);
static void decompose_transposed(
    int rnk_n, const INT *n, const INT *blk,
    int rnk_pm, const MPI_Comm *comms_pm,
    unsigned trafo_flag,
    INT *local_n, INT *local_start); 
static void decompose(
    const INT *pn, const INT *block,
    int rnk_pm, const MPI_Comm *comms_pm,
    INT *local_n, INT *local_start);
static void get_coords(
    int rnk_pm, const MPI_Comm *comms_pm,
    int *coords_pm);


/* call these routines with transp_flag
 * PFFT_TRANSPOSED_IN or PFFT_TRANSPOSED_OUT */


INT  PX(local_size_partrafo_transposed)(
    int rnk_n, const INT *n, const INT *ni, const INT *no,
    INT howmany, const INT *iblock, const INT *oblock,
    int rnk_pm, MPI_Comm *comms_pm,
    unsigned transp_flag, const unsigned *trafo_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  unsigned trafo_flag = trafo_flags[rnk_pm];
  INT mem=1, mem_tmp;
  INT Nb, tuple_size, N0, N1, h0, h1, hm, blk0, blk1, N, Ni, No;
  INT *pni, *pno;
  X(r2r_kind) kind, *kinds=NULL;
  
  /* get initial and final data decomposition */
  local_size_transposed(rnk_n, ni, no, iblock, oblock,
      rnk_pm, comms_pm, trafo_flag, transp_flag,
      local_ni, local_i_start, local_no, local_o_start);

  /* parameter of canonicalized trafo */
  Nb = initialize_Nb(rnk_pm, local_ni, local_no, transp_flag);
  tuple_size = calculate_tuple_size(
      rnk_n, ni, no, rnk_pm, howmany, trafo_flag, transp_flag);

  /* plan forward trafo of last dims */
  mem_tmp = PX(local_size_outrafo)(
      Nb, rnk_n - rnk_pm,
      &n[rnk_pm], &ni[rnk_pm], &no[rnk_pm],
      howmany, trafo_flag);
  mem = MAX(mem, mem_tmp);
  
  pni = PX(malloc_INT)(rnk_n);
  pno = PX(malloc_INT)(rnk_n);
  PX(physical_dft_size)(rnk_n, no, trafo_flag,
      pno);
  PX(physical_dft_size)(rnk_n, no, trafo_flag,
      pni);
  
  /* only trafo of last dimensions is r2c or c2r */
  if(trafo_flag & PFFTI_TRAFO_RDFT)
    trafo_flag = PFFTI_TRAFO_C2C;
  
  for(int t=0; t<rnk_pm; t++){
    PX(get_global_transp_param)(
        t, rnk_pm, pni, pno, local_ni, local_no,
        tuple_size, iblock, oblock, trafo_flag, transp_flag, 
	&N0, &N1, &h0, &h1, &hm, &blk0, &blk1);

    /* set hm to 1 since mem will be in units of real/complex */
    hm = 1;
    
    mem_tmp = PX(local_size_global_transp)(
        N0, N1, h0, h1, hm, blk0, blk1, comms_pm[rnk_pm-1-t]);
    mem = MAX(mem, mem_tmp);

    /* Read the parameters of the current step. */
    PX(get_outrafo_param)(
        t+1, rnk_pm, n, ni, no, local_ni, local_no, kinds, transp_flag, trafo_flags,
	&Nb, &N, &Ni, &No, &kind, &trafo_flag);

    mem_tmp = PX(local_size_outrafo)(
        Nb, 1, &N, &Ni, &No, tuple_size, trafo_flag);
    mem = MAX(mem, mem_tmp);
  }

  free(pni); free(pno);
  
  return mem;
}


void PX(plan_partrafo_transposed)(
    int rnk_n, const INT *n, const INT *ni, const INT *no,
    INT howmany, const INT *iblock, const INT *oblock,
    int rnk_pm, MPI_Comm *comms_pm,
    R *in_user, R *out_user, int sign, const X(r2r_kind) *kinds,
    unsigned transp_flag_user, const unsigned *trafo_flags,
    unsigned opt_flag, unsigned io_flag, unsigned fftw_flags, 
    outrafo_plan *trafos, gtransp_plan *remaps
    )
{
  unsigned transp_flag = transp_flag_user;
  unsigned trafo_flag = trafo_flags[rnk_pm];
  int s;
  INT *local_ni, *local_i_start, *local_no, *local_o_start;
  INT Nb, tuple_size, N0, N1, h0, h1, hm, blk0, blk1, N, Ni, No;
  INT *pni, *pno;
  R *in=in_user, *out=out_user;
  X(r2r_kind) kind;
 
  /* perform last trafos of backward plan in-place on 'out' in order to preserve input */
  if(transp_flag_user & PFFT_TRANSPOSED_IN)
    if(~io_flag & PFFT_DESTROY_INPUT)
      in=out;

  local_ni = PX(malloc_INT)(rnk_n);
  local_no = PX(malloc_INT)(rnk_n);
  local_i_start = PX(malloc_INT)(rnk_n);
  local_o_start = PX(malloc_INT)(rnk_n);

  /* get initial and final data decomposition */
  local_size_transposed(rnk_n, ni, no, iblock, oblock,
      rnk_pm, comms_pm, trafo_flag, transp_flag,
      local_ni, local_i_start, local_no, local_o_start);

  /* parameter of canonicalized trafo */
  Nb = initialize_Nb(rnk_pm, local_ni, local_no, transp_flag);
  tuple_size = calculate_tuple_size(
      rnk_n, ni, no, rnk_pm, howmany, trafo_flag, transp_flag);

  /* plan forward trafo of last dims */
  s = (transp_flag_user & PFFT_TRANSPOSED_IN) ? 2*rnk_pm+1 : 0;
  trafos[s] = PX(plan_outrafo)(
      Nb, rnk_n - rnk_pm,
      &n[rnk_pm], &ni[rnk_pm], &no[rnk_pm],
      howmany, in, out, sign, (kinds!=NULL) ? &kinds[rnk_pm] : NULL,
      trafo_flag, transp_flag, opt_flag, fftw_flags);

  /* perform last trafos of forward plan in-place on 'out' in order to preserve input */
  if(transp_flag_user & PFFT_TRANSPOSED_OUT)
    if(~io_flag & PFFT_DESTROY_INPUT)
      in=out;

  pni = PX(malloc_INT)(rnk_n);
  pno = PX(malloc_INT)(rnk_n);
  PX(physical_dft_size)(rnk_n, ni, trafo_flag,
      pni);
  PX(physical_dft_size)(rnk_n, no, trafo_flag,
      pno);
  
  /* only trafo of last dimensions is r2c or c2r */
  if(trafo_flag & PFFTI_TRAFO_RDFT)
    trafo_flag = PFFTI_TRAFO_C2C;
  
  for(int t=0; t<rnk_pm; t++){
    s = (transp_flag_user & PFFT_TRANSPOSED_IN) ? (2*rnk_pm)-t-1 : t;

    PX(get_global_transp_param)(
        t, rnk_pm, pni, pno, local_ni, local_no,
        tuple_size, iblock, oblock, trafo_flag, transp_flag,
	&N0, &N1, &h0, &h1, &hm, &blk0, &blk1);

    remaps[s] = PX(plan_global_transp)(
        N0, N1, h0, h1, hm, blk0, blk1, comms_pm[rnk_pm-1-t],
	out, in, transp_flag, fftw_flags);
    
    s = (transp_flag_user & PFFT_TRANSPOSED_IN) ? (2*rnk_pm+2)-t-2 : t+1;

    /* Read the parameters of the current step. */
    PX(get_outrafo_param)(
        t+1, rnk_pm, n, ni, no, local_ni, local_no, kinds, transp_flag, trafo_flags,
	&Nb, &N, &Ni, &No, &kind, &trafo_flag);

    /* last trafo without local transpose */
    if(t == rnk_pm-1)
      transp_flag = PFFT_TRANSPOSED_NONE;
    
    /* perform first trafo of backward plan out-of-place in order to preserve input */
    if(t == rnk_pm-1)
      if(transp_flag_user & PFFT_TRANSPOSED_IN)
        if(~io_flag & PFFT_DESTROY_INPUT)
          in=in_user;

    trafos[s] = PX(plan_outrafo)(
        Nb, 1, &N, &Ni, &No, tuple_size,
        in, out, sign, &kind, 
        trafo_flag, transp_flag, opt_flag, fftw_flags);
  }

  free(pni); free(pno);
  free(local_ni); free(local_i_start);
  free(local_no); free(local_o_start);
}


static INT initialize_Nb(
    int rnk_pm, const INT *local_ni, const INT *local_no,
    unsigned transp_flag
    )
{
  const INT *ln;

  /* TRANSPOSED_IN is planned in backward direction.
   * Therefore switch input and output parameters. */
  ln = (transp_flag & PFFT_TRANSPOSED_IN) ? local_no : local_ni;
  
  INT Nb = ln[0];
  for(int t=1; t<rnk_pm; t++)
    Nb *= ln[t];
  return Nb;
}


static INT calculate_tuple_size(
    int rnk_n, const INT *ni, const INT *no, int rnk_pm, INT howmany,
    unsigned trafo_flag, unsigned transp_flag
    )
{
  const INT *n;
  INT tpl = howmany;

  /* TRANSPOSED_IN is planned in backward direction.
   * Therefore switch input and output parameters. */
  n = (transp_flag & PFFT_TRANSPOSED_OUT) ? no : ni;
  
  if(rnk_pm+1 < rnk_n)
    tpl *= PX(physical_dft_size_1d)(n[rnk_n-1], trafo_flag);

  /* only trafo of last dimensions is r2c or c2r */
  if(trafo_flag & PFFTI_TRAFO_RDFT){
    trafo_flag &= ~PFFTI_TRAFO_RDFT;
    trafo_flag ^= PFFTI_TRAFO_C2C;
  }

  for(int t=rnk_pm+1; t<rnk_n-1; t++)
      tpl *= n[t];
  
  return tpl;
}



static void local_size_transposed(
    int rnk_n, const INT *ni, const INT *no,
    const INT *iblock, const INT *oblock,
    int rnk_pm, const MPI_Comm *comms_pm,
    unsigned trafo_flag, unsigned transp_flag,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  if( transp_flag & PFFT_TRANSPOSED_OUT ){
    decompose_nontransposed(rnk_n, ni, iblock, rnk_pm, comms_pm, trafo_flag,
        local_ni, local_i_start);
    decompose_transposed(rnk_n, no, oblock, rnk_pm, comms_pm, trafo_flag,
        local_no, local_o_start);
  } else { /* TRANSPOSED IN */
    decompose_transposed(rnk_n, ni, iblock, rnk_pm, comms_pm, trafo_flag,
        local_ni, local_i_start);
    decompose_nontransposed(rnk_n, no, oblock, rnk_pm, comms_pm, trafo_flag,
        local_no, local_o_start);
  }
}

static void decompose_nontransposed(
    int rnk_n, const INT *n, const INT *blk,
    int rnk_pm, const MPI_Comm *comms_pm,
    unsigned trafo_flag,
    INT *local_n, INT *local_start 
    )
{
  INT *pn = PX(malloc_INT)(rnk_n); 
  PX(physical_dft_size)(rnk_n, n, trafo_flag,
    pn);

  /* init all dims to undistributed case */
  for(int t=0; t<rnk_n; t++){
    local_n[t] = pn[t];
    local_start[t] = 0;
  }

  /* all dims from 0 to rnk_pm-1 are distributed */
  decompose(pn, blk, rnk_pm, comms_pm,
      local_n, local_start);

  free(pn);
}

static void decompose_transposed(
    int rnk_n, const INT *n, const INT *blk,
    int rnk_pm, const MPI_Comm *comms_pm,
    unsigned trafo_flag,
    INT *local_n, INT *local_start 
    )
{
  INT *pn = PX(malloc_INT)(rnk_n); 
  PX(physical_dft_size)(rnk_n, n, trafo_flag,
    pn);

  /* init all dims to undistributed case */
  for(int t=0; t<rnk_n; t++){
    local_n[t] = pn[t];
    local_start[t] = 0;
  }

  /* transposed distribution shifts by one dimension */
  decompose(pn+1, blk, rnk_pm, comms_pm,
      local_n+1, local_start+1);

  free(pn);
}


/* calculate block sizes from physical array size
 * for the 'rnk_pm' distributed dimensions */
static void decompose(
    const INT *pn, const INT *block,
    int rnk_pm, const MPI_Comm *comms_pm,
    INT *local_n, INT *local_start
    )
{
  int *coords_pm = PX(malloc_int)(rnk_pm);

  get_coords(rnk_pm, comms_pm,
      coords_pm);

  PX(decompose)(pn, block, rnk_pm, coords_pm,
      local_n, local_start);

  free(coords_pm);
}


static void get_coords(
    int rnk_pm, const MPI_Comm *comms_pm,
    int *coords_pm
    )
{
  int rnk;
  
  for(int t=0; t<rnk_pm; t++){
    MPI_Comm_rank(comms_pm[t], &rnk);
    MPI_Cart_coords(comms_pm[t], rnk, 1, &coords_pm[t]); 
  }
}

