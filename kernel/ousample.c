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
    ousam_plan_1d ths);

static int illegal_params(
    int rnk, const INT *ni, const INT *no,
    unsigned ousam_flag);
static INT local_size_ousam_1d(
    INT n0, INT n1i, INT n1o, INT howmany,
    unsigned trafo_flag);
static ousam_plan_1d plan_ousam_1d(
    INT n0, INT n1i, INT n1o, INT howmany,
    R *in, R *out, unsigned trafo_flag, unsigned ousam_flag);
static int is_complex_data(
    unsigned trafo_flag, unsigned ousam_flag);
static void execute_embed(
    ousam_plan_1d ths);
static void execute_trunc(
    ousam_plan_1d ths);
// static INT logical_dft_size(
//     INT n, unsigned trafo_flag);
static INT padding_size(
    INT n, unsigned trafo_flag, unsigned ousam_flag);



/* generalizes 'local_size_ousam_1d' to multiple dimensions */
INT PX(local_size_ousam_dd)(
    INT nb, int rnk, const INT *ni, const INT *no, INT howmany, 
    unsigned trafo_flag
    )
{
  INT n0, tuple_size, n1i, n1o, mem=0;
  INT *pno;

  pno = (INT*) malloc(sizeof(INT) * (size_t) rnk);

  PX(physical_dft_size)(rnk, no, trafo_flag,
      pno);

  for(int t=0; t<rnk; t++){
    tuple_size = howmany;
    for(int s=rnk-t; s<rnk; s++)
      tuple_size *= pno[s];

    n0 = nb;
    for(int s=0; s<rnk-t-1; s++)
      n0 *= ni[s];
    
    n1i = ni[rnk-t-1];
    n1o = no[rnk-t-1];

    mem = local_size_ousam_1d(
        n0, n1i, n1o, tuple_size, trafo_flag);

    /* only first onedimensional trafo is r2c or c2r */
    if(trafo_flag & PFFTI_TRAFO_RDFT)
      trafo_flag = PFFTI_TRAFO_C2C;
  }

  free(pno);

  return MAX(mem,1);
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
    R *in, R *out, unsigned trafo_flag, unsigned ousam_flag
    )
{
  INT n0, tuple_size, n1i, n1o;
  INT *pno;
  ousam_plan_dd ths;

  if( illegal_params(rnk, ni, no, ousam_flag) )
    return NULL;

  ths = ousam_dd_mkplan(rnk);
  
  pno = (INT*) malloc(sizeof(INT) * (size_t) rnk);

  PX(physical_dft_size)(rnk, no, trafo_flag,
      pno);

  /* plan 1d ousam 'rnk' times from last to first dimension */
  for(int t=0; t<rnk; t++){
    tuple_size = howmany;
    for(int s=rnk-t; s<rnk; s++)
      tuple_size *= pno[s];

    n0 = nb;
    for(int s=0; s<rnk-t-1; s++)
      n0 *= ni[s];
    
    n1i = ni[rnk-t-1];
    n1o = no[rnk-t-1];

    /* handle transposed layout for first dimension */
    if( (ousam_flag & PFFTI_OUSAM_TRANSPOSED) && (t==rnk-1) ){
      tuple_size *= n0;
      n0 = 1;
    }
  
    ths->ousam_1d[t] = plan_ousam_1d(
        n0, n1i, n1o, tuple_size, in, out, trafo_flag, ousam_flag);

    /* only first onedimensional trafo is r2c or c2r */
    if(trafo_flag & PFFTI_TRAFO_RDFT)
      trafo_flag = PFFTI_TRAFO_C2C;
  }

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


static INT local_size_ousam_1d(
    INT n0, INT n1i, INT n1o, INT howmany,
    unsigned trafo_flag
    )
{
  INT pni, pno, mem=1;

  pni = PX(physical_dft_size_1d)(n1i, trafo_flag);
  pno = PX(physical_dft_size_1d)(n1o, trafo_flag);

  mem *= n0 * MAX(pni, pno) * howmany;

  return mem;
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
    R *in, R *out, unsigned trafo_flag, unsigned ousam_flag
    )
{
  INT pni, pno;
  ousam_plan_1d ths;

  /* nothing to do and inplace */
  if( (n1i == n1o) && (in == out) )
    return NULL;
  
  ths = ousam_1d_mkplan();

//   INT lni = logical_dft_size(n1i, trafo_flag);
//   INT lno = logical_dft_size(n1o, trafo_flag);

  pni = PX(physical_dft_size_1d)(n1i, trafo_flag);
  pno = PX(physical_dft_size_1d)(n1o, trafo_flag);

  /* double phys. size for r2c input and c2r output */
  if( trafo_flag & PFFTI_TRAFO_RDFT )
    if( ! is_complex_data(trafo_flag, ousam_flag) ){
      pni *= 2; pno *= 2;
    }
  
  /* double 'howmany' for complex data type */
  if( is_complex_data(trafo_flag, ousam_flag) )
    howmany *= 2;
  
  ths->N0 = n0;
  ths->N1i = pni;
  ths->N1o = pno;

  /* calculate padding for r2c input and c2r output */
  ths->Pi = padding_size(n1i, trafo_flag, ousam_flag);
  ths->Po = padding_size(n1o, trafo_flag, ousam_flag);

  /* FIXME: Do not need Cr anymore */
  if(ousam_flag & PFFTI_OUSAM_EMBED){
    ths->Cl = pni - ths->Pi;
    ths->Cm = (pno - ths->Po) - (pni - ths->Pi);
    ths->Cr = 0;
  } else {
    ths->Cl = pno - ths->Po;
    ths->Cm = (pni - ths->Pi) - (pno - ths->Po);
    ths->Cr = 0;
  } 

  /* include howmany into parameters */
  ths->N1i *= howmany;
  ths->N1o *= howmany;
  ths->Cl  *= howmany;
  ths->Cm  *= howmany;
  ths->Cr  *= howmany;
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


// static INT logical_dft_size(
//     INT n, unsigned trafo_flag
//     )
// {
//   if( trafo_flag & PFFTI_TRAFO_RE00 ){
//     return 2*(n - 1);
//   } else if( trafo_flag & PFFTI_TRAFO_RO00 ) {
//     return 2*(n + 1);
//   } else if( trafo_flag & PFFTI_TRAFO_R2R ){
//     return 2*n;
//   } else {
//     return n;
//   }
// }

static INT padding_size(
    INT n, unsigned trafo_flag, unsigned ousam_flag
    )
{
  if( trafo_flag & PFFTI_TRAFO_RDFT )
    if( ! is_complex_data(trafo_flag, ousam_flag) )
      return 2*PX(physical_dft_size_1d)(n, trafo_flag) - n;
  
  return 0;
}


void PX(execute_ousam_dd)(
    ousam_plan_dd ths
    )
{
  if(ths == NULL)
    return;

  for(int t=0; t<ths->rnk; t++)
    execute_ousam_1d(ths->ousam_1d[t]);
}

static void execute_ousam_1d(
    ousam_plan_1d ths
    )
{
  if(ths == NULL)
    return;

  if( ths->ousam_flag & PFFTI_OUSAM_EMBED )
    execute_embed(ths);
  else
    execute_trunc(ths);
}


/* oversampling in one local dimension:
 * split the array at Cl and add Cm zeros in the gap */
static void execute_embed(
    ousam_plan_1d ths
    )
{
  INT mi, mo, k0, k1;
  const INT N0 = ths->N0, N1i = ths->N1i, N1o = ths->N1o;
  const INT Cl = ths->Cl, Cm = ths->Cm, Cr = ths->Cr; 
  const INT Pi = ths->Pi, Po = ths->Po;
  R *const in = ths->in, *const out = ths->out;

  if( (Cm == 0) && (Pi == 0) && (Po == 0) ){
    if( in != out )
      memcpy(in, out, sizeof(R) * (size_t) (N0*N1i));
    return;
  } 

  if( in != out ){
    mi = mo = 0;
    for(k0 = 0; k0 < N0; k0++){
      /* copy first half */
      for(k1 = 0; k1 < Cl; k1++)
        out[mo++] = in[mi++];
      
      /* fill gap with zeros */
      memset(&out[mo], 0, sizeof(R) * (size_t) Cm);
      mo += Cm;

      /* copy last half */
      for(k1 = 0; k1 < Cr; k1++)
        out[mo++] = in[mi++];

      /* padding for inplace r2c */
      mi += Pi; mo += Po;
    }
  } else {
    /* copy form last to first element to allow inplace execute */
    mi = N0 * N1i - 1;
    mo = N0 * N1o - 1;
    
    for(k0 = 0; k0 < N0; k0++){
      /* padding for r2c */
      mi -= Pi; mo -= Po;
  
      /* copy last half of array */
      for(k1 = 0; k1 < Cr; k1++)
        out[mo--] = in[mi--];
      
      /* fill gap with zeros */
      mo -= Cm;
      memset(&out[mo+1], 0, sizeof(R) * (size_t) Cm);
  
      /* special case:  N0 == 1
       * first half of array stays where it is */
      if( (N0 != 1) ){
        /* copy first half */
        for(k1 = 0; k1 < Cl; k1++)
          out[mo--] = in[mi--];
      }
    }
  }
}


/* undersampling in one local dimension:
 * skip Cm coefficients starting at Cl */
static void execute_trunc(
    ousam_plan_1d ths
    )
{
  INT mi, mo, k0, k1;
  const INT N0 = ths->N0, N1i = ths->N1i;
  const INT Cl = ths->Cl, Cm = ths->Cm, Cr = ths->Cr; 
  const INT Pi = ths->Pi, Po = ths->Po;
  R *const in = ths->in, *const out = ths->out;

  if( (Cm == 0) && (Pi == 0) && (Po == 0) ){
    if( in != out )
      memcpy(in, out, sizeof(R) * (size_t) (N0*N1i));
    return;
  } 

  mi = mo = 0;
  if( (in == out) && (N0 == 1) ){
    /* special case:  N0 == 1 and inplace
     *  first half of array stays where it is */
    mi += Cl + Cm; mo += Cl;

    for(k1 = 0; k1 < Cr; k1++)
      out[mo++] = in[mi++];
  } else{
    for(k0 = 0; k0 < N0; k0++){
      /* copy first half of array */
      for(k1 = 0; k1 < Cl; k1++)
        out[mo++] = in[mi++];
    
      /* copy last half of array */
      mi += Cm;
      for(k1 = 0; k1 < Cr; k1++)
        out[mo++] = in[mi++];

      /* padding for inplace r2c */
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




