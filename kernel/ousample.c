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
#include <string.h> /* memcpy, memset */

static ousam_plan_dd ousam_dd_mkplan(
    int rnk);
static ousam_plan_1d ousam_1d_mkplan(void);
static void ousam_1d_rmplan(
    ousam_plan_1d ths);
static void execute_ousam_1d(
    ousam_plan_1d ths,
    R *planned_in, R *planned_out,
    R *executed_in, R *executed_out);

static int illegal_params(
    int rnk, const INT *ni, const INT *no,
    unsigned ousam_flag);
static ousam_plan_1d plan_ousam_1d(
    INT n0, INT n1i, INT n1o, INT howmany,
    R *in, R *out, unsigned trafo_flag,
   unsigned si_flag, unsigned ousam_flag);
static int is_complex_data(
    unsigned trafo_flag, unsigned ousam_flag);
static void execute_embed(
    ousam_plan_1d ths,
    R *planned_in, R *planned_out,
    R *executed_in, R *executed_out);
static void execute_trunc(
    ousam_plan_1d ths,
    R *planned_in, R *planned_out,
    R *executed_in, R *executed_out);


INT PX(local_size_ousam_dd)(
    INT nb, int rnk, const INT *ni, const INT *no, INT howmany, 
    unsigned trafo_flag
    )
{
  INT *pni = (INT*) malloc(sizeof(INT) * (size_t) rnk);
  INT *pno = (INT*) malloc(sizeof(INT) * (size_t) rnk);

  PX(physical_dft_size)(rnk, ni, trafo_flag,
      pni);
  PX(physical_dft_size)(rnk, no, trafo_flag,
      pno);

  INT memi = nb * PX(prod_INT)(rnk, pni) * howmany;
  INT memo = nb * PX(prod_INT)(rnk, pno) * howmany;

  free(pni);
  free(pno);

  return MAX(memi, memo);
}


/* generalizes 'plan_ousam_1d' to multiple dimensions */
/* default:
 * nb x ni0 x ni1 x ... nir x hm -> nb x no0 x no1 x ... nor x hm 
 * transposed (PFFTI_OUSAM_TRANSPOSED):
 * ni0 x nb x ni1 x ... nir x hm -> no0 x nb x no1 x ... nor x hm 
 * 
 * */

ousam_plan_dd PX(plan_ousam_dd)(
    INT nb, int rnk, const INT *ni, const INT *no, INT howmany, 
    R *in, R *out, unsigned trafo_flag_user, unsigned si_flag, unsigned ousam_flag
    )
{
  INT n0, tuple_size, n1i, n1o;
  INT *pni, *pno;
  unsigned trafo_flag = trafo_flag_user;
  ousam_plan_dd ths;

  if( illegal_params(rnk, ni, no, ousam_flag) )
    return NULL;

  ths = ousam_dd_mkplan(rnk);
  
  pni = (INT*) malloc(sizeof(INT) * (size_t) rnk);
  pno = (INT*) malloc(sizeof(INT) * (size_t) rnk);

  PX(physical_dft_size)(rnk, ni, trafo_flag,
      pni);
  PX(physical_dft_size)(rnk, no, trafo_flag,
      pno);

  /* Index 'r' determines the order of plan execution.
   * Index 't' determines the dimensional order 1d embed/trunc. 
   * Note that 't' runs backward for embed and forward for trunc. */  
  for(int r=0; r<rnk; r++){
    /* plan 1d trunc 'rnk' times from first to last dimension */
    /* plan 1d embed 'rnk' times from last to first dimension */
    int t = (ousam_flag & PFFTI_OUSAM_EMBED) ? rnk-1-r : r;

    /* First 1d-r2c-embed creates padding. Afterward, we use C2C trunc on the physical array size. */
    if(trafo_flag_user & PFFTI_TRAFO_R2C)
      trafo_flag = (t == rnk-1) ? trafo_flag_user : PFFTI_TRAFO_C2C;

    /* First 'rnk-1' 1d-c2r-truncations work like C2C on the physical array size.
     * Last 1d-c2r-truncation deletes the padding */
    if(trafo_flag_user & PFFTI_TRAFO_C2R)
      trafo_flag = (t == rnk-1) ? trafo_flag_user : PFFTI_TRAFO_C2C;

    tuple_size = howmany;
    for(int s=t+1; s<rnk; s++)
      tuple_size *= (ousam_flag & PFFTI_OUSAM_EMBED) ? pno[s] : pni[s];

    n0 = nb;
    for(int s=0; s<t; s++)
      n0 *= (ousam_flag & PFFTI_OUSAM_EMBED) ? ni[s] : no[s];
    
    n1i = ni[t];
    n1o = no[t];

    /* handle transposed layout for first dimension */
    if( (ousam_flag & PFFTI_OUSAM_TRANSPOSED) && (t==0) ){
      tuple_size *= n0;
      n0 = 1;
    }
 
    ths->ousam_1d[r] = plan_ousam_1d(
        n0, n1i, n1o, tuple_size, in, out, trafo_flag, si_flag, ousam_flag);
  }

  free(pni);
  free(pno);

  return ths;
}


static int illegal_params(
    int rnk, const INT *ni, const INT *no,
    unsigned ousam_flag
    )
{
  if( ousam_flag & PFFTI_OUSAM_EMBED ){
    /* check ni < no for all dims  */
    for(int t=0; t<rnk; t++)
      if( ni[t] > no[t] )
        return 1;
  } else if( ousam_flag & PFFTI_OUSAM_TRUNC ){
    /* check ni > no for all dims  */
    for(int t=0; t<rnk; t++)
      if( ni[t] < no[t] )
        return 1;
  } else {
    /* no ousam flag was set */
    return 1;
  }

  return 0;
}

static int work_on_r2c_input(
    unsigned trafo_flag, unsigned ousam_flag
    )
{
  if(trafo_flag & PFFTI_TRAFO_R2C)
    if(ousam_flag & PFFTI_OUSAM_EMBED)
      return 1;

  return 0;
}

static int work_on_c2r_output(
    unsigned trafo_flag, unsigned ousam_flag
    )
{
  if(trafo_flag & PFFTI_TRAFO_C2R)
    if(ousam_flag & PFFTI_OUSAM_TRUNC)
      return 1;

  return 0;
}

/* 1D embed/truncate:
 * n0 x n1i x h -> n0 x n1o x h */
/* case n1i < n1o:
 * map array of size n0 x n1i x h into bigger array
 * of size n0 x n1o x h and fill the gaps with zeros */
/* case n1i > n1o:
 * truncate array of size n0 x n1i x h to smaller array 
 * n0 x n1o x h and put it contiguous into memory */

static ousam_plan_1d plan_ousam_1d(
    INT n0, INT n1i, INT n1o, INT howmany,
    R *in, R *out, unsigned trafo_flag,
    unsigned si_flag, unsigned ousam_flag
    )
{
  INT pni, pno;
  ousam_plan_1d ths;

  /* Return empty plan if nothing to do and inplace.
   * Attention: ousam generates the padding for r2c and c2r trafos also for n1i == n1o. */
  if( (n1i == n1o) && (in == out) )
    if( !work_on_r2c_input(trafo_flag, ousam_flag) )
      if( !work_on_c2r_output(trafo_flag, ousam_flag) )
        return NULL;
  
  ths = ousam_1d_mkplan();

  pni = PX(physical_dft_size_1d)(n1i, trafo_flag);
  pno = PX(physical_dft_size_1d)(n1o, trafo_flag);

  /* double 'howmany' for complex data type */
  if( is_complex_data(trafo_flag, ousam_flag) )
    howmany *= 2;
 
  ths->N0 = n0;
  ths->N1i = n1i;
  ths->N1o = n1o;

  /* default: no padding (only necessary for r2c and c2r) */
  ths->Pi = 0;
  ths->Po = 0;

  if(ousam_flag & PFFTI_OUSAM_EMBED)
  {
    ths->Zl = 0;
    ths->D  = ths->N1i;
    ths->Zr = ths->N1o - ths->N1i;

    if(si_flag & PFFT_SHIFTED_IN)
      ths->Zl = ths->Zr = ths->Zr/2;
  }

  if(ousam_flag & PFFTI_OUSAM_TRUNC)
  {
    ths->Zl = 0;
    ths->D = ths->N1o;
    ths->Zr = ths->N1i - ths->N1o;

    if(si_flag & PFFT_SHIFTED_OUT)
      ths->Zl = ths->Zr = ths->Zr/2;
  }

  /* special case r2c */
  if( trafo_flag & PFFTI_TRAFO_R2C )
  {
    /* use physical array size for r2c output */
    if( ousam_flag & PFFTI_OUSAM_TRUNC )
    {
      ths->N1i = pni;
      ths->N1o = pno;
      ths->Zl = pni - pno;
      ths->D  = pno;
      ths->Zr = 0;
    }

    /* generate padding for r2c input */
    if( ousam_flag & PFFTI_OUSAM_EMBED ){
      ths->Zl = 0;
      ths->D  = n1i;
      ths->Zr = n1o - n1i;
      ths->Po  = 2*pno-n1o;
      ths->N1o = 2*pno;

      /* skip padding generation if the inputs are already padded */
      if( trafo_flag & PFFTI_TRAFO_PADDED ){
        ths->Pi  = 2*pni-n1i;
        ths->N1i = 2*pni;
      }
    
      if( si_flag & PFFT_SHIFTED_IN )
        ths->Zl = ths->Zr = ths->Zr/2;
    }
  }

  /* special case c2r */
  if( trafo_flag & PFFTI_TRAFO_C2R )
  {
    /* use physical array size for c2r input */
    if( ousam_flag & PFFTI_OUSAM_EMBED )
    {
      ths->N1i = pni;
      ths->N1o = pno;
      ths->Zl = pno - pni;
      ths->D  = pni;
      ths->Zr = 0;
    }

    /* truncate padding for c2r output */
    if( ousam_flag & PFFTI_OUSAM_TRUNC ){
      ths->Zl = 0;
      ths->D  = n1o;
      ths->Zr = n1i - n1o;
      ths->Pi  = 2*pni-n1i;
      ths->N1i = 2*pni;
   
      /* do not truncate padding, if the user explicitly wants it */
      if( trafo_flag & PFFTI_TRAFO_PADDED ){
        ths->Po = 2*pno-n1o;
        ths->N1o = 2*pno;
      }
    
      if( si_flag & PFFT_SHIFTED_OUT )
        ths->Zl = ths->Zr = ths->Zr/2;
    }
  }

  /* include howmany into parameters */
  ths->N1i *= howmany;
  ths->N1o *= howmany;
  ths->Zl  *= howmany;
  ths->D   *= howmany;
  ths->Zr  *= howmany;
  ths->Pi  *= howmany;
  ths->Po  *= howmany;

  ths->in = in;
  ths->out = out;
  ths->ousam_flag = ousam_flag;

  return ths;
}

static int is_complex_data(
    unsigned trafo_flag, unsigned ousam_flag
    )
{
  if( trafo_flag & PFFTI_TRAFO_C2C )
    return 1;

  if( trafo_flag & PFFTI_TRAFO_R2C)
    if( ousam_flag & PFFTI_OUSAM_TRUNC )
      return 1;

  if( trafo_flag & PFFTI_TRAFO_C2R )
    if( ousam_flag & PFFTI_OUSAM_EMBED )
      return 1;
 
  return 0;
}


void PX(execute_ousam_dd)(
    ousam_plan_dd ths,
    R *planned_in, R *planned_out,
    R *executed_in, R *executed_out
    )
{
  if(ths == NULL)
    return;

  for(int t=0; t<ths->rnk; t++)
    execute_ousam_1d(ths->ousam_1d[t], planned_in, planned_out, executed_in, executed_out);
}

static void execute_ousam_1d(
    ousam_plan_1d ths,
    R *planned_in, R *planned_out,
    R *executed_in, R *executed_out
    )
{
  if(ths == NULL)
    return;

  if( ths->ousam_flag & PFFTI_OUSAM_EMBED )
    execute_embed(ths, planned_in, planned_out, executed_in, executed_out);
  else
    execute_trunc(ths, planned_in, planned_out, executed_in, executed_out);
}


/* oversampling in one local dimension:
 * split the array at Cl and add Cm zeros in the gap */
static void execute_embed(
    ousam_plan_1d ths,
    R *planned_in, R *planned_out,
    R *executed_in, R *executed_out
    )
{
  INT mi, mo, k0, k1;
  const INT N0 = ths->N0, N1i = ths->N1i, N1o = ths->N1o;
  const INT Zl = ths->Zl, Zr = ths->Zr; 
  const INT D = ths->D; 
  const INT Pi = ths->Pi, Po = ths->Po;
  R *in  = (ths->in  == planned_in)  ? executed_in  : executed_out;
  R *out = (ths->out == planned_out) ? executed_out : executed_in;

  if( (Zl == 0) && (Zr == 0) && (Pi == 0) && (Po == 0) ){
    if( in != out )
      if(N0*N1i > 0)
        memcpy(in, out, sizeof(R) * (size_t) (N0*N1i));
    return;
  } 

  if( in != out ){
    mi = mo = 0;
    for(k0 = 0; k0 < N0; k0++){
      /* fill in left zero block */
      if(Zl > 0) memset(&out[mo], 0, sizeof(R) * (size_t) Zl);
      mo += Zl;

      /* copy data block */
      for(k1 = 0; k1 < D; k1++)
        out[mo++] = in[mi++];
      
      /* fill in right zero block */
      if(Zr > 0) memset(&out[mo], 0, sizeof(R) * (size_t) Zr);
      mo += Zr;

      /* padding for r2c/c2r */
      mi += Pi; mo += Po;
    }
  } else {
    /* copy from last to first element to allow inplace execute */
    mi = N0 * N1i - 1;
    mo = N0 * N1o - 1;
    
    for(k0 = 0; k0 < N0; k0++){
      /* padding for r2c/c2r */
      mi -= Pi; mo -= Po;
  
      /* fill in right zero block */
      mo -= Zr;
      if(Zr > 0) memset(&out[mo+1], 0, sizeof(R) * (size_t) Zr);
  
      /* special case:  N0 == 1 and Zl == 0
       * left data block stays where it is */
      if( (N0 != 1) || (Zl != 0) ){
        /* copy data block */
        for(k1 = 0; k1 < D; k1++)
          out[mo--] = in[mi--];

        /* fill in left zero block */
        mo -= Zl;
        if(Zl > 0) memset(&out[mo+1], 0, sizeof(R) * (size_t) Zl);
      }
    }
  }
}


/* undersampling in one local dimension:
 * skip Zl coefficients starting at 0 and Zr coefficients starting at D */
static void execute_trunc(
    ousam_plan_1d ths,
    R *planned_in, R *planned_out,
    R *executed_in, R *executed_out
    )
{
  INT mi, mo, k0, k1;
  const INT N0 = ths->N0, N1i = ths->N1i;
  const INT Zl = ths->Zl, Zr = ths->Zr;
  const INT D = ths->D; 
  const INT Pi = ths->Pi, Po = ths->Po;
  R *in  = (ths->in  == planned_in)  ? executed_in  : executed_out;
  R *out = (ths->out == planned_out) ? executed_out : executed_in;

  if( (Zl == 0) && (Zr == 0) && (Pi == 0) && (Po == 0) ){
    if( in != out )
      if(N0*N1i > 0)
        memcpy(in, out, sizeof(R) * (size_t) (N0*N1i));
    return;
  } 

  mi = mo = 0;
  if( (in == out) && (N0 == 1) && (Zl == 0) ){
    /* special case:  (N0 == 1) and (Zl == 0) and inplace
     * data block stays where it is */
    mi += D; mo += D;

    /* skip right zero block */
    mi += Zr;
  } else {
    for(k0 = 0; k0 < N0; k0++){
      /* skip left zero block */
      mi += Zl;

      /* copy data block */
      for(k1 = 0; k1 < D; k1++)
        out[mo++] = in[mi++];
    
      /* skip right zero block */
      mi += Zr;

      /* skip padding for inplace r2c */
      mi += Pi; mo += Po;
    }
  }
}


static ousam_plan_dd ousam_dd_mkplan(
    int rnk
    )
{
  ousam_plan_dd ths = (ousam_plan_dd) malloc(sizeof(ousam_plan_dd_s));
  ths->ousam_1d = (ousam_plan_1d*) malloc(sizeof(ousam_plan_1d) * (size_t)rnk);

  ths->rnk = rnk;  
        
  return ths;
}

static ousam_plan_1d ousam_1d_mkplan(void)
{
  return (ousam_plan_1d) malloc(sizeof(ousam_plan_1d_s));
}

void PX(ousam_dd_rmplan)(
    ousam_plan_dd ths
    )
{
  /* plan was already destroyed or never initialized */
  if(ths==NULL)
    return;
  
  /* destroy all ousam plans and free the array of pointers */
  if(ths->ousam_1d != NULL){
    for(int t=0; t<ths->rnk; t++)
      ousam_1d_rmplan(ths->ousam_1d[t]);
    free(ths->ousam_1d);
  }

  /* free struct */
  free(ths);
}

static void ousam_1d_rmplan(
    ousam_plan_1d ths
    )
{
  /* plan was already destroyed or never initialized */
  if(ths==NULL)
    return;
  
  free(ths);
}




