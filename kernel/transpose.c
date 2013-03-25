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

static void get_global_transp_param(
    int step, int rnk_pm, const INT *ni, const INT *no,
    const INT *local_ni, const INT *local_no,
    INT tuple_size, const INT *iblock, const INT *oblock,
    unsigned trafo_flag,
    INT *N0, INT *N1, INT *h0, INT *h1,
    INT *hm, INT *blk0, INT *blk1);
static gtransp_plan gtransp_mkplan(
    void);

#if PFFT_DEBUG_GTRANSP
static gtransp_dbg gtransp_mkdbg(
    INT *N, INT hm, INT *blk, R *in, R *out,
    MPI_Comm comm, unsigned fftw_flags);
static void gtransp_rmdbg(
    gtransp_dbg ths);
static void print_dbg(
    gtransp_dbg ths);
static void print_parameters(
    INT N0, INT N1, INT hm, INT blk0, INT blk1, R *in, R *out,
    INT local_N0, INT local_N0_start,
    INT local_N1, INT local_N1_start,
    MPI_Comm comm, unsigned fftw_flags);
#endif


/* MPI Transpose for real data. Use hm==2 for complex data.  */
/* Input:
 *   (N0 x h0) x (N1/P x h1) x hm  (default)
 *   (N1/P x h1) x (N0 x h0) x hm (TRANSPOSED_IN)
 * Output:
 *   (N1 x h1) x (N0/P x h0) x hm (default)
 *   (N0/P x h0) x (N1 x h1) x hm (TRANSPOSED_OUT)
 *  
 *  This is a collective routine.
 *  N0, h0, N1, h1 must be equal on all calling processes.
 *   */


void PX(get_global_transp_param)(
    int step, int rnk_pm, const INT *ni, const INT *no,
    const INT *local_ni, const INT *local_no,
    INT tuple_size, const INT *iblock, const INT *oblock,
    unsigned trafo_flag, unsigned transp_flag,
    INT *N0, INT *N1, INT *h0, INT *h1,
    INT *hm, INT *blk0, INT *blk1
    )
{
  /* TRANSPOSED_IN is planned in backward direction.
   * Therefore switch input and output parameters. */
  if(transp_flag & PFFT_TRANSPOSED_IN)
    get_global_transp_param(
        step, rnk_pm, no, ni, local_no, local_ni,
        tuple_size, oblock, iblock, trafo_flag,
	N1, N0, h1, h0, hm, blk1, blk0);
  else
    get_global_transp_param(
        step, rnk_pm, ni, no, local_ni, local_no,
        tuple_size, iblock, oblock, trafo_flag,
	N0, N1, h0, h1, hm, blk0, blk1);
}


static void get_global_transp_param(
    int step, int rnk_pm, const INT *ni, const INT *no,
    const INT *local_ni, const INT *local_no,
    INT tuple_size, const INT *iblock, const INT *oblock,
    unsigned trafo_flag,
    INT *N0, INT *N1, INT *h0, INT *h1,
    INT *hm, INT *blk0, INT *blk1
    )
{
  *N0 = no[rnk_pm     - step];
  *N1 = ni[rnk_pm - 1 - step];

  *h0 = 1;
  for(int t=0; t<step; t++)
    *h0 *= local_no[rnk_pm - t];
  for(int t=0; t<rnk_pm-1-step; t++)
    *h0 *= local_ni[t];
  *h1 = tuple_size;
  
  /* double hm for complex data */
  *hm = (~trafo_flag & PFFTI_TRAFO_R2R) ? 2 : 1;

  *blk0 = oblock[rnk_pm - 1 - step];
  *blk1 = iblock[rnk_pm - 1 - step];
}


INT PX(local_size_global_transp)(
    INT N0, INT N1, INT h0, INT h1, INT hm, INT blk0, INT blk1,
    MPI_Comm comm 
    )
{
  INT N[2], blk[2];
  INT dummy0, dummy1, dummy2, dummy3;
  INT mem = 1; /* never return 0 */  

  N[0] = N1 * h1; blk[0] = blk1 * h1;
  N[1] = N0 * h0; blk[1] = blk0 * h0;

  /* For strange distributions (e.g. 4 data points on 3 processes)
   * all processes in a row get no data. */
  if(N[0]*N[1]*h0*h1*hm==0)
    return 1;
  
  mem = XM(local_size_many_transposed)(
      2, N, hm, blk[0], blk[1], comm,
      &dummy0, &dummy1, &dummy2, &dummy3);

  return MAX(mem, 1);
}


gtransp_plan PX(plan_global_transp)(
    INT N0, INT N1, INT h0, INT h1, INT hm, INT blk0, INT blk1,
    MPI_Comm comm, R *in, R *out,
    unsigned transp_flag, unsigned fftw_flags
    )
{
  gtransp_plan ths = NULL;
  
  INT N[2], blk[2];

  N[0] = N1 * h1; blk[0] = blk1 * h1;
  N[1] = N0 * h0; blk[1] = blk0 * h0;

  /* For strange distributions (e.g. 4 data points on 3 processes)
   * all processes in a row get no data. */
  if(N[0]*N[1]*h0*h1*hm==0)
    return NULL;

  ths = gtransp_mkplan();

  /* PFFTs and FFTWs transpose flags are complementary */
  if( ~transp_flag & PFFT_TRANSPOSED_OUT)
    fftw_flags |= FFTW_MPI_TRANSPOSED_OUT;
  if( ~transp_flag & PFFT_TRANSPOSED_IN)
    fftw_flags |= FFTW_MPI_TRANSPOSED_IN;

#if PFFT_BUGFIX_FORGET_PARALLEL_FFTW_WISDOM
  X(forget_wisdom)();
#endif

#if PFFT_DEBUG_GTRANSP
  ths->dbg = gtransp_mkdbg(N, hm, blk, in, out, comm, fftw_flags);
#endif

  ths->plan = XM(plan_many_transpose)(
      N[0], N[1], hm, blk[0], blk[1], in, out,
      comm, fftw_flags);

  return ths;
}


static gtransp_plan gtransp_mkplan(
    void
    )
{
  gtransp_plan ths = (gtransp_plan) malloc(sizeof(gtransp_plan_s));
  
  /* initialize to NULL for easy checks */
  ths->plan=NULL;

  /* initialize debug info */
#if PFFT_DEBUG_GTRANSP
  ths->dbg = NULL;
#endif
  
  return ths;
}

void PX(gtransp_rmplan)(
    gtransp_plan ths
    )
{
  /* plan was already destroyed or never initialized */
  if(ths==NULL)
    return;

  /* take care of unsuccessful FFTW planing */
  if(ths->plan != NULL)
    X(destroy_plan)(ths->plan);

#if PFFT_DEBUG_GTRANSP
  if(ths->dbg != NULL)
    gtransp_rmdbg(ths->dbg);
#endif
  
  /* free memory */
  free(ths);
  /* ths=NULL; would be senseless, since we can not change the pointer itself */
}


void PX(execute_gtransp)(
    gtransp_plan ths
    )
{
#if PFFT_DEBUG_GTRANSP
  static int counter=0;
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  R lsum, gsum;
  INT n_total;

  if(!myrank) fprintf(stderr, "\n");
  if(!myrank){
    if(ths != NULL){
      if(ths->plan != NULL){
        fprintf(stderr, "PFFT_DBG_GTRANSP: counter = %d\n", counter);
        print_dbg(ths->dbg);
      } else
        fprintf(stderr, "PFFT_DBG_GTRANSP: nothing to do for FFTW plan\n");
    }
  }
 
  n_total = (ths != NULL) ? ths->dbg->local_N0 * ths->dbg->N1 * ths->dbg->hm : 0;
  
  /* Checksum inputs */ 
  lsum=0.0;
  if(ths != NULL)
    if(ths->plan != NULL)
      for(INT k=0; k<n_total; k++)
        lsum += fabs(ths->dbg->in[k]);
  MPI_Reduce(&lsum, &gsum, 1, PFFT_MPI_REAL_TYPE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(ths != NULL)
    if(ths->plan != NULL)
      if(!myrank) fprintf(stderr, "PFFT_DBG_GTRANSP: counter = %d, Checksum(in) = %e\n", counter, gsum);
#endif

  /* Global transposition */ 
  if(ths != NULL)
    if(ths->plan != NULL)
      X(execute)(ths->plan);
    
#if PFFT_DEBUG_GTRANSP
  n_total = (ths != NULL) ? ths->dbg->N0 * ths->dbg->local_N1 * ths->dbg->hm : 0;

  /* Checksum outputs */ 
  lsum=0.0;
  if(ths != NULL)
    if(ths->plan != NULL)
      for(INT k=0; k<n_total; k++)
        lsum += fabs(ths->dbg->out[k]);
  MPI_Reduce(&lsum, &gsum, 1, PFFT_MPI_REAL_TYPE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(ths != NULL)
    if(ths->plan != NULL)
      if(!myrank) fprintf(stderr, "PFFT_DBG_GTRANSP: counter = %d, Checksum(out) = %e\n", counter, gsum);

//   if(counter==3){
//     int np, rnk, np_world;
//     MPI_Comm_size(ths->dbg->comm, &np);
//     MPI_Comm_size(MPI_COMM_WORLD, &np_world);
//     MPI_Comm_rank(ths->dbg->comm, &rnk);
//     R *lsums = (R*) malloc(sizeof(R)*np);
//     R *gsums = (R*) malloc(sizeof(R)*np);
//       
//     /* calculate horizontal checksums */ 
//     for(int t=0; t<np; t++)
//       lsums[t] = 0;
//     lsums[ ths->dbg->local_N1_start/ths->dbg->blk1 ] = lsum;
//           
//     MPI_Allreduce(lsums, gsums, np, PFFT_MPI_REAL_TYPE, MPI_SUM, MPI_COMM_WORLD);
// 
//     for(int t=0; t<np; t++){
//       if(!myrank)
//         fprintf(stderr, "gsums[%d] = %e\n", t, gsums[t]);
//     }
//     for(int t=0; t<np/2; t++){
//       if(!myrank)
//         fprintf(stderr, "sum of gsums[%d] = %e\n", t, gsums[2*t]+gsums[2*t+1]);
//     }
//     free(lsums); free(gsums);
//     MPI_Barrier(MPI_COMM_WORLD);
//     
//     /* calculate vertical checksums */
//     MPI_Reduce(&lsum, &gsum, 1, PFFT_MPI_REAL_TYPE, MPI_SUM, 0, ths->dbg->comm);
//     if(!rnk) fprintf(stderr, "local_N1 = %td, local_N1_start = %td, vertical gsum = %e\n",
//        ths->dbg->local_N1, ths->dbg->local_N1_start, gsum);
//   }
  
  counter++;
#endif
}

#if PFFT_DEBUG_GTRANSP
static gtransp_dbg gtransp_mkdbg(
    INT *N, INT hm, INT *blk, R *in, R *out,
    MPI_Comm comm, unsigned fftw_flags
    )
{
  gtransp_dbg ths = (gtransp_dbg) malloc(sizeof(gtransp_dbg_s));

  ths->N0 = N[0];
  ths->N1 = N[1];
  ths->hm = hm;
  ths->blk0 = blk[0];
  ths->blk1 = blk[1];
  ths->in = in;
  ths->out = out;
  MPI_Comm_dup(comm, &ths->comm);
  ths->fftw_flags = fftw_flags;
  
  ths->mem = XM(local_size_many_transposed)(
      2, N, hm, blk[0], blk[1], comm,
      &ths->local_N0, &ths->local_N0_start,
      &ths->local_N1, &ths->local_N1_start);
  
  return ths;
}

static void gtransp_rmdbg(
    gtransp_dbg ths
    )
{
  /* plan was already destroyed or never initialized */
  if(ths==NULL)
    return;
  
  /* free memory */
  free(ths);
  /* ths=NULL; would be senseless, since we can not change the pointer itself */
}

static void print_dbg(
    gtransp_dbg ths
    )
{
  print_parameters(
      ths->N0, ths->N1, ths->hm, ths->blk0, ths->blk1, ths->in, ths->out,
      ths->local_N0, ths->local_N0_start, ths->local_N1, ths->local_N1_start,
      ths->comm, ths->fftw_flags);
}

static void print_parameters(
    INT N0, INT N1, INT hm, INT blk0, INT blk1, R *in, R *out,
    INT local_N0, INT local_N0_start,
    INT local_N1, INT local_N1_start,
    MPI_Comm comm, unsigned fftw_flags
    )
{
  int np;
  MPI_Comm_size(comm, &np);

  fprintf(stderr, "PFFT_DBG_GTRANSP: N[0] = %td; N[1] = %td; hm = %td; blk[0] = %td; blk[1] = %td; np = %d;\n",
      N0, N1, hm, blk0, blk1, np);
  fprintf(stderr, "PFFT_DBG_GTRANSP: local_N0 = %td; local_N0_start = %td; local_N1 = %td; local_N1_start = %td;\n",
     local_N0, local_N0_start, local_N1, local_N1_start);
  fprintf(stderr, "PFFT_DBG_GTRANSP: fftw_flags = 0");
  if(fftw_flags & FFTW_ESTIMATE)
    fprintf(stderr, "| FFTW_ESTIMATE");
  if(fftw_flags & FFTW_MEASURE)
    fprintf(stderr, "| FFTW_MEASURE");
  if(fftw_flags & FFTW_PATIENT)
    fprintf(stderr, "| FFTW_PATIENT");
  if(fftw_flags & FFTW_EXHAUSTIVE)
    fprintf(stderr, "| FFTW_EXHAUSTIVE");
  if(fftw_flags & FFTW_MPI_TRANSPOSED_IN)
    fprintf(stderr, "| FFTW_MPI_TRANSPOSED_IN");
  if(fftw_flags & FFTW_MPI_TRANSPOSED_OUT)
    fprintf(stderr, "| FFTW_MPI_TRANSPOSED_OUT");
  fprintf(stderr, ";\n");

  fprintf(stderr, "PFFT_DBG_GTRANSP: mem = XM(local_size_many_transposed)(\n");
  fprintf(stderr, "PFFT_DBG_GTRANSP:     2, N, hm, blk[0], blk[1], comm,\n");
  fprintf(stderr, "PFFT_DBG_GTRANSP:     &dummy0, &dummy1, &dummy2, &dummy3);\n");
  fprintf(stderr, "PFFT_DBG_GTRANSP: ths = XM(plan_many_transpose)(\n");
  fprintf(stderr, "PFFT_DBG_GTRANSP:     N[0], N[1], hm, blk[0], blk[1], in,");
  if(in==out)
    fprintf(stderr, " in,\n");
  else
    fprintf(stderr, " out,\n");
  fprintf(stderr, "PFFT_DBG_GTRANSP:     comm, fftw_flags);\n");
}
#endif

