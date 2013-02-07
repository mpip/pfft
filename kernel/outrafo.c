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

static void get_outrafo_param(
    int step, int rnk_pm,
    const INT *n, const INT *ni, const INT *no,
    const INT *local_ni, const INT *local_no,
    const X(r2r_kind) *kind_in, const unsigned *trafo_flags,
    INT *Nb, INT *N, INT *Ni, INT *No,
    X(r2r_kind) *kind, unsigned *trafo_flag);
static outrafo_plan outrafo_mkplan(void);

/* combine multidimensional serial trafo (c2c, r2c, c2r, r2r)
 * with over- and/or under sampling */
/* input array:
 *   nb x ni0 x ni1 x ... x nir x howmany
 * oversample to array of size:
 *   nb x n0 x n1 x ... x nr x howmany
 * plan trafo of last dimensions:
 *   nb x N0 x N1 x ... x Nr x howmany
 * truncate to array of size:
 *   nb x no0 x no1 x ... x nor x howmany
 */ 


void PX(get_outrafo_param)(
    int step, int rnk_pm,
    const INT *n, const INT *ni, const INT *no,
    const INT *local_ni, const INT *local_no,
    const X(r2r_kind) *kinds, unsigned transp_flag, const unsigned *trafo_flags,
    INT *Nb, INT *N, INT *Ni, INT *No,
    X(r2r_kind) *kind, unsigned *trafo_flag
    )
{
  /* TRANSPOSED_IN is planned in backward direction.
   * Therefore switch input and output parameters. */
  if(transp_flag & PFFT_TRANSPOSED_IN)
    get_outrafo_param(
        step, rnk_pm, n, ni, no, local_no, local_ni, kinds, trafo_flags,
	Nb, N, Ni, No, kind, trafo_flag);
  else
    get_outrafo_param(
        step, rnk_pm, n, ni, no, local_ni, local_no, kinds, trafo_flags,
	Nb, N, Ni, No, kind, trafo_flag);
}


static void get_outrafo_param(
    int step, int rnk_pm,
    const INT *n, const INT *ni, const INT *no,
    const INT *local_ni, const INT *local_no,
    const X(r2r_kind) *kinds, const unsigned *trafo_flags,
    INT *Nb, INT *N, INT *Ni, INT *No,
    X(r2r_kind) *kind, unsigned *trafo_flag
    )
{
  *Nb = 1;
  for(int t=0; t<step; t++)
    *Nb *= local_no[rnk_pm - t];
  for(int t=0; t<rnk_pm-step; t++)
    *Nb *= local_ni[t];

  *N  =  n[rnk_pm - step];
  *Ni = ni[rnk_pm - step];
  *No = no[rnk_pm - step];

  /* pick any kind, if not needed */
  *kind = (kinds != NULL) ? kinds[rnk_pm - step] : FFTW_REDFT00;

  *trafo_flag = trafo_flags[rnk_pm - step];
}

INT PX(local_size_outrafo)(
    INT nb, int rnk, const INT *n, const INT *ni, const INT *no, INT howmany,
    unsigned trafo_flag
    )
{
  INT mem=1, mem_tmp; /* Never return zero */

  /* embed */ 
  mem_tmp = PX(local_size_ousam_dd)(nb, rnk, ni, n, howmany, 
      trafo_flag);
  mem = MAX(mem, mem_tmp);

  /* serial trafo */
  mem_tmp = PX(local_size_sertrafo)(nb, rnk, n, howmany, trafo_flag);
  mem = MAX(mem, mem_tmp);

  /* trunc */ 
  mem_tmp = PX(local_size_ousam_dd)(nb, rnk, n, no, howmany, 
      trafo_flag);
  mem = MAX(mem, mem_tmp);

  return mem;
}


outrafo_plan PX(plan_outrafo)(
    INT nb, int rnk, const INT *n, const INT *ni, const INT *no, INT howmany,
    R *in, R *out, int sign, const X(r2r_kind) *kinds,
    unsigned trafo_flag, unsigned transp_flag,
    unsigned opt_flag, unsigned fftw_flags
    )
{
  unsigned io_flag = 0; /* implementation not finished */
  unsigned ousam_flag;
  outrafo_plan ths = outrafo_mkplan();

  /* plan embed */ 
  ousam_flag = (transp_flag & PFFT_TRANSPOSED_IN) ?
    PFFTI_OUSAM_TRANSPOSED : 0;
  ths->embed = PX(plan_ousam_dd)(nb, rnk, ni, n, howmany, 
    in, in, trafo_flag, PFFTI_OUSAM_EMBED | ousam_flag);

  /* plan serial trafo */
  ths->trafo = PX(plan_sertrafo)(
    nb, rnk, n, howmany, in, out, sign, kinds,
    trafo_flag, transp_flag, io_flag, opt_flag, fftw_flags);

  /* plan trunc */ 
  ousam_flag = (transp_flag & PFFT_TRANSPOSED_OUT) ?
    PFFTI_OUSAM_TRANSPOSED : 0;
  ths->trunc = PX(plan_ousam_dd)(nb, rnk, n, no, howmany, 
    out, out, trafo_flag, PFFTI_OUSAM_TRUNC | ousam_flag);

  /* give back NULL if all contained palns are NULL */
  if( (ths->trafo == NULL) && (ths->embed == NULL) && (ths->trunc == NULL) ){
    PX(outrafo_rmplan)(ths); 
    return NULL;
  }

  return ths;
}


void PX(execute_outrafo)(
    outrafo_plan ths
    )
{
  if(ths == NULL)
    return;

  PX(execute_ousam_dd)(ths->embed);
  PX(execute_sertrafo)(ths->trafo);
  PX(execute_ousam_dd)(ths->trunc);
}


static outrafo_plan outrafo_mkplan(void)
{
  outrafo_plan ths = (outrafo_plan) malloc(sizeof(outrafo_plan_s));
  
  ths->embed = NULL;
  ths->trafo = NULL;
  ths->trunc = NULL;
  
  return ths;
}    


void PX(outrafo_rmplan)(
    outrafo_plan ths
    )
{
  /* plan was already destroyed or never initialized */
  if(ths==NULL)
    return;

  /* destroy trafo and ousam plans */
  PX(ousam_dd_rmplan)(ths->embed);
  PX(sertrafo_rmplan)(ths->trafo);
  PX(ousam_dd_rmplan)(ths->trunc);

  /* finally free mem of the plan itself */
  free(ths);
}

