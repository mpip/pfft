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

/* Global infos about procmesh are only enabled in debug mode
 * Otherwise we do not use any global variables. */
#if PFFT_DEBUG_GVARS
  extern MPI_Comm *gdbg_comm_cart;
  extern int gdbg_rnk_pm;
  extern MPI_Comm *gdbg_comms_pm;
#endif

#if PFFT_DEBUG_SERTRAFO
  #define PFFT_DEBUG_SERTRAFO_PTR00 , &(ths0->dbg[0]) 
  #define PFFT_DEBUG_SERTRAFO_PTR01 , &(ths0->dbg[1]) 
  #define PFFT_DEBUG_SERTRAFO_PTR10 , &(ths1->dbg[0]) 
  #define PFFT_DEBUG_SERTRAFO_PTR11 , &(ths1->dbg[1]) 
#else
  #define PFFT_DEBUG_SERTRAFO_PTR00 
  #define PFFT_DEBUG_SERTRAFO_PTR01 
  #define PFFT_DEBUG_SERTRAFO_PTR10 
  #define PFFT_DEBUG_SERTRAFO_PTR11 
#endif

static sertrafo_plan plan_sertrafo_p(
    INT nb, int rnk, const INT *n, INT howmany,
    R *in0, R *out0, R *in1, R *out1,
    int sign, const X(r2r_kind) *kinds,
    unsigned trafo_flag, unsigned transp_flag,
    unsigned opt_flag, unsigned fftw_flags, 
    double *time);
static sertrafo_plan plan_sertrafo_pt(
    INT nb, int rnk, const INT *n, INT howmany,
    R *in0, R *out0, R *in1, R *out1, int sign,
    const X(r2r_kind) *kinds, unsigned trafo_flag,
    unsigned transp_flag0, unsigned transp_flag1,
    unsigned opt_flag, unsigned fftw_flags, 
    double *time);
static sertrafo_plan plan_remap_only(
    INT nb, int rnk, const INT *n, INT howmany, R *in, R *out,
    unsigned trafo_flag, unsigned transp_flag, unsigned fftw_flags);
static X(plan) plan_trafo(
    INT nb, int rnk, const INT *n, INT howmany,
    R *in, R *out, int sign, const X(r2r_kind) *kinds,
    unsigned trafo_flag, unsigned transp_flag, unsigned fftw_flags
#if PFFT_DEBUG_SERTRAFO
    , sertrafo_dbg *dbg_ptr
#endif
    );
static X(plan) plan_remap(
    INT nb, int rnk, const INT *n, INT howmany, R *in, R *out,
    unsigned trafo_flag, unsigned transp_flag,
    unsigned fftw_flags, unsigned ioarray_flag
#if PFFT_DEBUG_SERTRAFO
    , sertrafo_dbg *dbg_ptr
#endif
    );
static void malloc_and_fill_dims_trafo(
    INT nb, int rnk, const INT *n, INT howmany,
    unsigned trafo_flag, unsigned transp_flag,
    int *dims_rnk_ptr, X(iodim64) **dims_ptr,
    int *howmany_rnk_ptr, X(iodim64) **howmany_dims_ptr);
static void malloc_and_fill_dims_remap(
    INT nb, int rnk, const INT *n, INT howmany, 
    unsigned trafo_flag, unsigned transp_flag, unsigned ioarray_flag,
    int *howmany_rnk_ptr, X(iodim64) **howmany_dims_ptr);
static void malloc_and_calculate_strides_remap(
    INT nb, int rnk, const INT *pn, INT howmany,
    unsigned transp_flag,
    INT **sz_ptr, INT **is_ptr, INT **os_ptr);
static void malloc_and_calculate_strides_trafo(
    INT nb, int rnk, const INT *n, INT howmany,
    unsigned trafo_flag, unsigned transp_flag,
    INT **sz_ptr, INT **is_ptr, INT **os_ptr);
static void exchange_INT(
    INT *a, INT *b);
static double measure_time(
    X(plan) plan);
static sertrafo_plan sertrafo_mkplan(
    void);
static int needs_transpose(
    unsigned transp_flag);

#if PFFT_DEBUG_SERTRAFO
static sertrafo_dbg sertrafo_mkdbg(
    int is_fft, INT nb, int rnk, const INT *n, INT howmany,
    R *in, R *out,
    unsigned trafo_flag, unsigned transp_flag, unsigned fftw_flags, unsigned ioarray_flag,
    int dims_rnk, const X(iodim64) *dims, int howmany_rnk, const X(iodim64) *howmany_dims);
void sertrafo_rmdbg(
    sertrafo_dbg ths);
static void print_vec(
    int rnk, INT *vec);
static void print_dbg(
    sertrafo_dbg ths);
INT calc_n_total(
    sertrafo_dbg ths);
#endif

/* input array:
 *   nb x n0 x n1 x ... x nr x howmany (default, PFFT_TRANSPOSED_NONE)
 *   n0 x nb x n1 x ... x nr x howmany (PFFT_TRANSPOSED_IN)
 * transformed output array:
 *   nb x N0 x N1 x ... x Nr x howmany (default, PFFT_TRANSPOSED_NONE)
 *   N0 x nb x N1 x ... x Nr x howmany (PFFT_TRANSPOSED_OUT)
 * !!! first dimension 'nb' will not be transformed !!!
 * - trafo can be c2c, r2c, c2r, r2r
 * - special attention for r2c, c2r padding
 * - r2r trafo needs r2r_kind for every dimension
 *   (different boundary conditions possible)
 * - sign only needed for c2c, otherwise sign==0
 * - kind only needed for r2r, otherwise kind==NULL
 *
 * If trafo_flag 'PFFTI_TRAFO_SKIP' is set, we execute only the
 * local transpose on the input array, but ommit the transform.
 * We use the transform type (C2C, R2C, C2R, R2R) to determine
 * the data type of the input array (and the correct strides).
 * 
 * If trafo_flag 'PFFTI_TRAFO_SKIP' AND transp_flag 'PFFT_TRANSPOSE_NONE'
 * are set, we still copy the input for out-of-place plans. 
 *
 * If trafo_flag 'PFFTI_TRAFO_PHANTOM' is set, we return a NULL pointer.
 * Sometimes its easier to call the planner with PFFTI_TRAFO_PHANTOM than
 * to handle the special case were nothing should be done.
 * 
 * You can use io_flags for detailed influence on the algorithms:
 * PFFT_BUFFERED_INPLACE: similar to a out-of-place plan whose outputs go back to 'in'
 * PFFT_DESTROY_INPUT:    out-of-place plan are allowed to overwrite the inputs (disabled on default)
 * PFFT_PRESERVE_INPUT:   cancels PFFT_DESTROY_INPUT
 * */

INT PX(local_size_sertrafo)(
    INT nb, int rnk, const INT *n, INT howmany,
    unsigned trafo_flag
    )
{
  INT *pn, mem;
  
  pn = PX(malloc_INT)(rnk);

  PX(physical_dft_size)(rnk, n, trafo_flag,
      pn);

  mem = nb * PX(prod_INT)(rnk, pn) * howmany;

  free(pn);

  return mem;
}


/* FFTW sometimes produces very slow plans for 1dFFT with different is and os.
 * Therefore we replace the FFTW 1dFFT by two FFTW plans, one 1d-FFT and one local_transposition.
 * We can choose which one of these plans comes first, which one works out-of-place
 * and which one has different is and os. This results in at most 8 different
 * plans which we compare by time. */

/* DEFAULT SETTINGS:
 *
 * R2C FFTs can not combine FFT and transposition.
 * => compute FFT and transposition in two different steps
 *
 * We want stride-1-FFTs to be the default whenever possible.
 * => T_OUT:           perform FFT before local transposition (maybe input has stride 1)
 *    T_NONE :         perform FFT first, no transposition needed (stride 1)
 *    T_IN:            perform local transposition before FFT (maybe the output has stride 1)
 *    T_IN and T_OUT:  perform FFT as second step, no transposition needed (stride 1 not possible)
 *
 * For symmetry reasons we want the default of a T_IN trafo to be symmetric to a T_OUT trafo.
 * => T_OUT:           perform first step out-of-place
 *    T_IN:            perform second step out-of-place
 * If no transposition is needed, we want the FFT to be out-place.
 * => T_NONE:          perform first step out-of-place 
 *    T_IN and T_OUT:  perform second step out-of-place
 *
 * However, if PFFT_TUNE is set all the other options are also checked by the planner.
 * */


/* At first we decide which plans work in-place or out-of-place. */
sertrafo_plan PX(plan_sertrafo)(
    INT nb, int rnk, const INT *n, INT howmany,
    R *in, R *out, int sign, const X(r2r_kind) *kinds,
    unsigned trafo_flag, unsigned transp_flag, unsigned io_flag,
    unsigned opt_flag, unsigned fftw_flags
    )
{

  double time0, time1;
  sertrafo_plan ths0=NULL, ths1=NULL;

  if(trafo_flag & PFFTI_TRAFO_PHANTOM)
    return NULL;

  if(trafo_flag & PFFTI_TRAFO_SKIP)
    return plan_remap_only(
        nb, rnk, n, howmany, in, out,
        trafo_flag, transp_flag, fftw_flags);

  /* try all allowed sequences of in-place and out-of-place plans */
  if(in == out) { /* in-place */
    ths0 = plan_sertrafo_p(
        nb, rnk, n, howmany, in, in, in, in, sign, kinds,
        trafo_flag, transp_flag, opt_flag, fftw_flags,
        &time0);

    /* assure that plan0 wins, no plan1 possible */
    time1 = time0 * 10.0;

  } else { /* in != out */
    if( io_flag & PFFT_BUFFERED_INPLACE ) { /* out-of-place with outputs back in input array */
      ths0 = plan_sertrafo_p(
          nb, rnk, n, howmany, in, out, out, in, sign, kinds,
          trafo_flag, transp_flag, opt_flag, fftw_flags,
          &time0);

      /* assure that plan0 wins, if plan1 is not possible */
      time1 = time0 * 10.0;

      /* skip extra plan for PFFT_NO_TUNE */
      if( opt_flag & PFFT_TUNE )
        ths1 = plan_sertrafo_p(
            nb, rnk, n, howmany, in, in, in, in, sign, kinds,
            trafo_flag, transp_flag, opt_flag, fftw_flags,
            &time1);

    } else { /* out-of-place */
      if(transp_flag & PFFT_TRANSPOSED_IN){
        /* input does not have stride-1, transpose first */
        /* default: perform second step out-of-place */
        ths0 = plan_sertrafo_p(
            nb, rnk, n, howmany, in, in, in, out, sign, kinds,
            trafo_flag, transp_flag, opt_flag, fftw_flags,
            &time0);

        /* assure that plan0 wins, if plan1 is not possible */
        time1 = time0 * 10.0;

        /* skip extra plan for PFFT_NO_TUNE */
        if( (~io_flag & PFFT_PRESERVE_INPUT) && (opt_flag & PFFT_TUNE) )
          /* compare to first step out-of-place */
          ths1 = plan_sertrafo_p(
              nb, rnk, n, howmany, in, out, out, out, sign, kinds,
              trafo_flag, transp_flag, opt_flag, fftw_flags,
              &time1);

      } else { /* input not transposed */
        /* default: perform first step out-of-place */
        ths0 = plan_sertrafo_p(
            nb, rnk, n, howmany, in, out, out, out, sign, kinds,
            trafo_flag, transp_flag, opt_flag, fftw_flags,
            &time0);

        /* assure that plan0 wins, if plan1 is not possible */
        time1 = time0 * 10.0;

        /* compare to second step out-of-place */
        /* skip extra plan for PFFT_NO_TUNE */
        if( (~io_flag & PFFT_PRESERVE_INPUT) && (opt_flag & PFFT_TUNE) )
          ths1 = plan_sertrafo_p(
              nb, rnk, n, howmany, in, in, in, out, sign, kinds,
              trafo_flag, transp_flag, opt_flag, fftw_flags,
              &time1);
      }
    }
  }
 
  /* choose best plan, delete the other one */
  if(time1 < time0){
    PX(sertrafo_rmplan)(ths0);
    return ths1;
  } else{
    PX(sertrafo_rmplan)(ths1);
    return ths0;
  }
}


/* only perform remap */
static sertrafo_plan plan_remap_only(
    INT nb, int rnk, const INT *n, INT howmany, R *in, R *out,
    unsigned trafo_flag, unsigned transp_flag, unsigned fftw_flags
    )
{
  sertrafo_plan ths=NULL;

  ths = sertrafo_mkplan();

#if PFFT_DEBUG_SERTRAFO
  ths->plan[0] = plan_remap(
      nb, rnk, n, howmany, in, out, trafo_flag, transp_flag, fftw_flags, PFFTI_ARRAY_INPUT,
      &ths->dbg[0]);
#else
  ths->plan[0] = plan_remap(
      nb, rnk, n, howmany, in, out, trafo_flag, transp_flag, fftw_flags, PFFTI_ARRAY_INPUT);
#endif
  ths->plan[1] = NULL;

  return ths;
}


/* Decide if FFT and local transposition are performed at once or in two steps.
 * default: calculate trafo and transposition in two different steps,
 *          since R2C trafo can not combine trafo and transposition */ 
static sertrafo_plan plan_sertrafo_p(
    INT nb, int rnk, const INT *n, INT howmany,
    R *in0, R *out0, R *in1, R *out1,
    int sign, const X(r2r_kind) *kinds,
    unsigned trafo_flag, unsigned transp_flag,
    unsigned opt_flag, unsigned fftw_flags,
    double *time
    )
{
  unsigned transp_flag0, transp_flag1;
  double time0, time1;
  sertrafo_plan ths0=NULL, ths1=NULL;

  if(transp_flag & PFFT_TRANSPOSED_IN){
    /* default: transpose during 1st plan */
      transp_flag0 = transp_flag;
      transp_flag1 = (transp_flag & PFFT_TRANSPOSED_OUT) ?
        PFFT_TRANSPOSED_IN | PFFT_TRANSPOSED_OUT : PFFT_TRANSPOSED_NONE;
  } else {
    /* default: transpose during 2nd plan */
    transp_flag0 = (transp_flag & PFFT_TRANSPOSED_IN) ?
      PFFT_TRANSPOSED_IN | PFFT_TRANSPOSED_OUT : PFFT_TRANSPOSED_NONE;
    transp_flag1 = transp_flag;
  }

  ths0 = plan_sertrafo_pt(
      nb, rnk, n, howmany, in0, out0, in1, out1, sign, kinds,
      trafo_flag, transp_flag0, transp_flag1, opt_flag, fftw_flags,
      &time0);
  
  /* assure that plan0 wins, if plan1 is not possible */
  time1 = time0 * 10.0;
  
  /* try this plan only for altered strides */
  /* skip extra plan for PFFT_NO_TUNE */
  if( (opt_flag & PFFT_TUNE) && needs_transpose(transp_flag) ){
    if(transp_flag & PFFT_TRANSPOSED_IN){
      /* compare to: transpose during 2nd plan */
      transp_flag0 = (transp_flag & PFFT_TRANSPOSED_IN) ?
        PFFT_TRANSPOSED_IN | PFFT_TRANSPOSED_OUT : PFFT_TRANSPOSED_NONE;
      transp_flag1 = transp_flag;
    } else {
      /* compare to: transpose during 1st plan */
        transp_flag0 = transp_flag;
        transp_flag1 = (transp_flag & PFFT_TRANSPOSED_OUT) ?
          PFFT_TRANSPOSED_IN | PFFT_TRANSPOSED_OUT : PFFT_TRANSPOSED_NONE;
    }

    ths1 = plan_sertrafo_pt(
        nb, rnk, n, howmany, in0, out0, in1, out1, sign, kinds,
        trafo_flag, transp_flag0, transp_flag1, opt_flag, fftw_flags,
        &time1);
  }
  
  /* choose best plan, delete the other one */
  if(time1 < time0){
    *time = time1;
    PX(sertrafo_rmplan)(ths0);
    return ths1;
  } else {
    *time = time0;
    PX(sertrafo_rmplan)(ths1);
    return ths0;
  }
}

static int needs_transpose(
    unsigned transp_flag
    )
{
  /* no transposition flag set */
  if(transp_flag == PFFT_TRANSPOSED_NONE)
    return 0;
  
  /* no transposition needed, if both flags are set */
  if(transp_flag & PFFT_TRANSPOSED_IN)
    if(transp_flag & PFFT_TRANSPOSED_OUT)
      return 0;
  
  /* need to transpose, if only one flag is set */
  return 1;
}


/* Decide if FFT or remap is executed first:
 * default: compute FFT first and remap second. */
static sertrafo_plan plan_sertrafo_pt(
    INT nb, int rnk, const INT *n, INT howmany,
    R *in0, R *out0, R *in1, R *out1, int sign,
    const X(r2r_kind) *kinds, unsigned trafo_flag,
    unsigned transp_flag0, unsigned transp_flag1,
    unsigned opt_flag, unsigned fftw_flags, 
    double *time
    ) 
{
  double time0=0.0, time1=0.0;
  sertrafo_plan ths0=NULL, ths1=NULL;
   
  ths0 = sertrafo_mkplan();

  if(transp_flag0 & PFFT_TRANSPOSED_IN){
    /* default: perform remap first, fft second */
    ths0->plan[0] = plan_remap(
        nb, rnk, n, howmany, in0, out0, trafo_flag, transp_flag0, fftw_flags, PFFTI_ARRAY_INPUT
        PFFT_DEBUG_SERTRAFO_PTR00);
    ths0->plan[1] = plan_trafo(
        nb, rnk, n, howmany, in1, out1, sign, kinds, trafo_flag, transp_flag1, fftw_flags
        PFFT_DEBUG_SERTRAFO_PTR01);

  } else {
    /* default: perform FFT first, remap second */
    ths0->plan[0] = plan_trafo(
        nb, rnk, n, howmany, in0, out0, sign, kinds, trafo_flag, transp_flag0, fftw_flags
        PFFT_DEBUG_SERTRAFO_PTR00);
    ths0->plan[1] = plan_remap(
        nb, rnk, n, howmany, in1, out1, trafo_flag, transp_flag1, fftw_flags, PFFTI_ARRAY_OUTPUT
        PFFT_DEBUG_SERTRAFO_PTR01);
  }

  /* only measure times if we need them for comparison to plan1 */
  if(opt_flag & PFFT_TUNE)
    time0 = measure_time(ths0->plan[0]) + measure_time(ths0->plan[1]);
  
  /* this plan differs for altered strides or out-of-place */
  ths1 = sertrafo_mkplan();
  if( opt_flag & PFFT_TUNE ){ /* skip extra plan for PFFT_NO_TUNE */
    if( (in0 != out1)
        || needs_transpose(transp_flag0)
        || needs_transpose(transp_flag1)
        )
    {
      if(transp_flag0 & PFFT_TRANSPOSED_IN){
        /* default: perform remap first, fft second */
        ths1->plan[0] = plan_trafo(
            nb, rnk, n, howmany, in0, out0, sign, kinds, trafo_flag, transp_flag0, fftw_flags
            PFFT_DEBUG_SERTRAFO_PTR10);
        ths1->plan[1] = plan_remap(
            nb, rnk, n, howmany, in1, out1, trafo_flag, transp_flag1, fftw_flags, PFFTI_ARRAY_OUTPUT
            PFFT_DEBUG_SERTRAFO_PTR11);
      } else {
        /* compare to: perform remap first, FFT second */
        ths1->plan[0] = plan_remap(
            nb, rnk, n, howmany, in0, out0, trafo_flag, transp_flag0, fftw_flags, PFFTI_ARRAY_INPUT
            PFFT_DEBUG_SERTRAFO_PTR10);
        ths1->plan[1] = plan_trafo(
            nb, rnk, n, howmany, in1, out1, sign, kinds, trafo_flag, transp_flag1, fftw_flags
            PFFT_DEBUG_SERTRAFO_PTR11);
      }
      time1 = measure_time(ths1->plan[0]) + measure_time(ths1->plan[1]);
    }
  }

  if( (opt_flag & PFFT_TUNE) && (time1 < time0) ){
    *time = time1;
    PX(sertrafo_rmplan)(ths0);
    return ths1;
  } else {
    *time = time0;
    PX(sertrafo_rmplan)(ths1);
    return ths0;
  }
}


static X(plan) plan_trafo(
    INT nb, int rnk, const INT *n, INT howmany,
    R *in, R *out, int sign, const X(r2r_kind) *kinds,
    unsigned trafo_flag, unsigned transp_flag, unsigned fftw_flags
#if PFFT_DEBUG_SERTRAFO
    , sertrafo_dbg *dbg_ptr
#endif
    )
{
  int dims_rnk, howmany_rnk;
  X(iodim64) *dims, *howmany_dims;
  X(plan) ths;
 
  /* R2C can not combine trafo and transposition */
  if( (trafo_flag & PFFTI_TRAFO_RDFT) && needs_transpose(transp_flag) )
    return NULL;

  /* R2R can not combine trafo and transposition */
  if( (trafo_flag & PFFTI_TRAFO_R2R) && needs_transpose(transp_flag) )
    return NULL;
  
  malloc_and_fill_dims_trafo(
    nb, rnk, n, howmany, trafo_flag, transp_flag,
    &dims_rnk, &dims, &howmany_rnk, &howmany_dims);

  /* choose appropriate fftw planner for trafo */
  if(trafo_flag & PFFTI_TRAFO_R2C){
    ths = X(plan_guru64_dft_r2c)(
        dims_rnk, dims, howmany_rnk, howmany_dims,
        in, (C*) out, fftw_flags);
  } else if(trafo_flag & PFFTI_TRAFO_C2R){
    ths = X(plan_guru64_dft_c2r)(
        dims_rnk, dims, howmany_rnk, howmany_dims,
        (C*) in, out, fftw_flags);
  } else if(trafo_flag & PFFTI_TRAFO_R2R){
    ths = X(plan_guru64_r2r)(
        dims_rnk, dims, howmany_rnk, howmany_dims,
	in, out, kinds, fftw_flags);
  } else {
    ths = X(plan_guru64_dft)(
        dims_rnk, dims, howmany_rnk, howmany_dims,
        (C*) in, (C*) out, sign, fftw_flags);
  }
  
#if PFFT_DEBUG_SERTRAFO
  *dbg_ptr = sertrafo_mkdbg(1, nb, rnk, n, howmany, in, out, 
      trafo_flag, transp_flag, fftw_flags, PFFTI_ARRAY_UNDEFINED,
      dims_rnk, dims, howmany_rnk, howmany_dims);
#endif

  free(dims); free(howmany_dims);

  return ths;
}

static void malloc_and_fill_dims_trafo(
    INT nb, int rnk, const INT *n, INT howmany,
    unsigned trafo_flag, unsigned transp_flag,
    int *dims_rnk_ptr, X(iodim64) **dims_ptr,
    int *howmany_rnk_ptr, X(iodim64) **howmany_dims_ptr
    )
{
  int dims_rnk, howmany_rnk;
  INT *sz, *is, *os;
  X(iodim64) *dims, *howmany_dims;
  
  malloc_and_calculate_strides_trafo(
    nb, rnk, n, howmany, trafo_flag, transp_flag,
    &sz, &is, &os);
    
  dims_rnk = rnk;
  dims = (X(iodim64)*) malloc(sizeof(X(iodim64)) * (size_t) dims_rnk);

  for(int t=0; t<rnk; t++){
    dims[t].n  = sz[t+1];
    dims[t].is = is[t+1];
    dims[t].os = os[t+1];
  }

  howmany_rnk = 2;
  howmany_dims = (X(iodim64)*) malloc(sizeof(X(iodim64)) * (size_t) howmany_rnk);

  /* one loop for nb */
  howmany_dims[0].n  = sz[0];
  howmany_dims[0].is = is[0];
  howmany_dims[0].os = os[0];

  /* one loop for howmany */
  howmany_dims[1].n  = sz[rnk+1];
  howmany_dims[1].is = is[rnk+1];
  howmany_dims[1].os = os[rnk+1];

  *dims_ptr = dims;
  *howmany_dims_ptr = howmany_dims;
  *dims_rnk_ptr = dims_rnk;
  *howmany_rnk_ptr = howmany_rnk;

  free(sz); free(is); free(os);
}

/* We use a r2r plan for r2c transposition. For the input
 * we double the physical array size and for the output we
 * double 'howmany'. Analogously done for c2r. */
static X(plan) plan_remap(
    INT nb, int rnk, const INT *n, INT howmany, R *in, R *out,
    unsigned trafo_flag, unsigned transp_flag,
    unsigned fftw_flags, unsigned ioarray_flag
#if PFFT_DEBUG_SERTRAFO
    , sertrafo_dbg *dbg_ptr
#endif
    )
{
  /* wrap FFTWs guru64 FFT planing, only reorder without FFT */
  int sign=1, howmany_rnk, dims_rnk=0;
  X(r2r_kind) *kinds=NULL;
  X(iodim64) *howmany_dims, *dims=NULL;
  X(plan) ths;
    
  malloc_and_fill_dims_remap(
      nb, rnk, n, howmany, trafo_flag, transp_flag, ioarray_flag,
      &howmany_rnk, &howmany_dims);

  /* choose appropriate fftw planner for remap */
  if(trafo_flag & PFFTI_TRAFO_C2C)
    ths = X(plan_guru64_dft)(
        dims_rnk, dims, howmany_rnk, howmany_dims,
        (C*) in, (C*) out, sign, fftw_flags);
  else
    ths = X(plan_guru64_r2r)(
        dims_rnk, dims, howmany_rnk, howmany_dims,
	in, out, kinds, fftw_flags);

#if PFFT_DEBUG_SERTRAFO
  *dbg_ptr = sertrafo_mkdbg(0, nb, rnk, n, howmany, in, out, 
      trafo_flag, transp_flag, fftw_flags, ioarray_flag,
      dims_rnk, dims, howmany_rnk, howmany_dims);
#endif
    
  free(howmany_dims);

  return ths;
}


static void malloc_and_fill_dims_remap(
    INT nb, int rnk, const INT *n, INT howmany, 
    unsigned trafo_flag, unsigned transp_flag, unsigned ioarray_flag,
    int *howmany_rnk_ptr, X(iodim64) **howmany_dims_ptr
    )
{
  int howmany_rnk;
  INT *sz, *is, *os, *pn;
  X(iodim64) *howmany_dims;

  pn = PX(malloc_INT)(rnk);
  PX(physical_dft_size)(rnk, n, trafo_flag, pn);

  if( trafo_flag & PFFTI_TRAFO_R2C ){
    if( ioarray_flag & PFFTI_ARRAY_OUTPUT )
      howmany *= 2;
    else
      pn[rnk-1] *= 2;
  }
  
  if( trafo_flag & PFFTI_TRAFO_C2R ){
    if( ioarray_flag & PFFTI_ARRAY_INPUT )
      howmany *= 2;
    else
      pn[rnk-1] *= 2;
  }

  malloc_and_calculate_strides_remap(
    nb, rnk, pn, howmany, transp_flag,
    &sz, &is, &os);
 
  howmany_rnk = rnk+2;
  howmany_dims = (X(iodim64)*) malloc(sizeof(X(iodim64)) * (size_t) howmany_rnk);

  for(int t=0; t<rnk+2; t++){
    howmany_dims[t].n  = sz[t];
    howmany_dims[t].is = is[t];
    howmany_dims[t].os = os[t];
  }

  *howmany_dims_ptr = howmany_dims;
  *howmany_rnk_ptr = howmany_rnk;

  free(sz); free(is); free(os); free(pn);
}

/* set array size 'sz' for array of size
 * nb x n[0] x n[1] x ... x n[rnk-1] x howmany,
 * set input strides 'is' corresponing to array of size
 * nb x n[0] x n[1] x ... x n[rnk-1] x howmany (default) or
 * n[0] x nb x n[1] x ... x n[rnk-1] x howmany (TRANSPOSED_IN),
 * set output strides 'os' corresponing to array of size
 * nb x n[0] x n[1] x ... x n[rnk-1] x howmany (default) or
 * n[0] x nb x n[1] x ... x n[rnk-1] x howmany (TRANSPOSED_OUT),
 */
static void malloc_and_calculate_strides_remap(
    INT nb, int rnk, const INT *pn, INT howmany,
    unsigned transp_flag,
    INT **sz_ptr, INT **is_ptr, INT **os_ptr
    )
{ 
  INT *sz, *is, *os, *pni, *pno;

  /* allocate return arrays */
  sz = PX(malloc_INT)(rnk+2);
  is = PX(malloc_INT)(rnk+2);
  os = PX(malloc_INT)(rnk+2);

  /* we need the physical array sizes to compute strides */
  pni = PX(malloc_INT)(rnk+2);
  pno = PX(malloc_INT)(rnk+2);
  
  /* set FFTW array size */
  sz[0] = nb;
  for(int t=0; t<rnk; t++)
    sz[t+1] = pn[t];
  sz[rnk+1] = howmany;

  /* set physical array size */
  for(int t=0; t<rnk+2; t++)
    pni[t] = pno[t] = sz[t];

  /* permute first two dims of array size */
  if(transp_flag & PFFT_TRANSPOSED_IN )
    exchange_INT(&pni[0], &pni[1]);
  if(transp_flag & PFFT_TRANSPOSED_OUT )
    exchange_INT(&pno[0], &pno[1]);

  /* calculate strides */
  is[rnk+1] = os[rnk+1] = 1;
  for(int t=rnk; t>=0; t--){
    is[t] = is[t+1] * pni[t+1];
    os[t] = os[t+1] * pno[t+1];
  }
  
  /* permute first two dims of strides */
  if(transp_flag & PFFT_TRANSPOSED_IN )
    exchange_INT(&is[0], &is[1]);
  if(transp_flag & PFFT_TRANSPOSED_OUT )
    exchange_INT(&os[0], &os[1]);

  free(pni); free(pno);

  *sz_ptr = sz;
  *is_ptr = is;
  *os_ptr = os;
}


static void malloc_and_calculate_strides_trafo(
    INT nb, int rnk, const INT *n, INT howmany,
    unsigned trafo_flag, unsigned transp_flag,
    INT **sz_ptr, INT **is_ptr, INT **os_ptr
    )
{ 
  INT *sz, *is, *os, *pni, *pno;

  /* allocate return arrays */
  sz = PX(malloc_INT)(rnk+2);
  is = PX(malloc_INT)(rnk+2);
  os = PX(malloc_INT)(rnk+2);

  /* we need the physical array sizes to compute strides */
  pni = PX(malloc_INT)(rnk+2);
  pno = PX(malloc_INT)(rnk+2);
  
  /* set FFTW array size */
  sz[0] = nb;
  for(int t=0; t<rnk; t++)
    sz[t+1] = n[t];
  sz[rnk+1] = howmany;

  /* set physical array size */
  for(int t=0; t<rnk+2; t++)
    pni[t] = pno[t] = sz[t];

  /* trafo of last dim maybe r2c or c2r */
  /* compute physical size of R2C for real input and complex output */
  if(trafo_flag & PFFTI_TRAFO_R2C){
    pni[rnk] = 2 * PX(physical_dft_size_1d)(sz[rnk], trafo_flag);
    pno[rnk] =  PX(physical_dft_size_1d)(sz[rnk], trafo_flag);
  }

  /* compute physical size of C2R for complex input and real output */
  if(trafo_flag & PFFTI_TRAFO_C2R){
    pni[rnk] = PX(physical_dft_size_1d)(sz[rnk], trafo_flag);
    pno[rnk] = 2 * PX(physical_dft_size_1d)(sz[rnk], trafo_flag);
  }

  /* permute first two dims of array size */
  if(transp_flag & PFFT_TRANSPOSED_IN )
    exchange_INT(&pni[0], &pni[1]);
  if(transp_flag & PFFT_TRANSPOSED_OUT )
    exchange_INT(&pno[0], &pno[1]);

  /* calculate strides */
  is[rnk+1] = os[rnk+1] = 1;
  for(int t=rnk; t>=0; t--){
    is[t] = is[t+1] * pni[t+1];
    os[t] = os[t+1] * pno[t+1];
  }
  
  /* permute first two dims of strides */
  if(transp_flag & PFFT_TRANSPOSED_IN )
    exchange_INT(&is[0], &is[1]);
  if(transp_flag & PFFT_TRANSPOSED_OUT )
    exchange_INT(&os[0], &os[1]);

  free(pni); free(pno);

  *sz_ptr = sz;
  *is_ptr = is;
  *os_ptr = os;
}

static void exchange_INT(
    INT *a, INT *b
    )
{
  INT buf = *a;
  *a = *b;
  *b = buf;
}

static double measure_time(
    X(plan) plan
    )
{
  double time=0;
 
  if(plan != NULL){
    /* 1st execution is often very slow, discard it */
    X(execute)(plan);
    /* measure 2nd execution */
    time = -MPI_Wtime();
    X(execute)(plan);
    time += MPI_Wtime();
  }

  return time;
}


static sertrafo_plan sertrafo_mkplan(
    void
    )
{
  sertrafo_plan ths = (sertrafo_plan) malloc(sizeof(sertrafo_plan_s));
  
  /* initialize to NULL for easy checks */
  for(int t=0; t<2; t++)
    ths->plan[t]=NULL;

  /* initialize debug info */
#if PFFT_DEBUG_SERTRAFO
  for(int t=0; t<2; t++)
    ths->dbg[t] = NULL;
#endif
  
  return ths;
}

void PX(sertrafo_rmplan)(
    sertrafo_plan ths
    )
{
  /* plan was already destroyed or never initialized */
  if(ths==NULL)
    return;

  /* take care of unsuccessful FFTW planing */
  for(int t=0; t<2; t++)
    if(ths->plan[t] != NULL)
      X(destroy_plan)(ths->plan[t]);

#if PFFT_DEBUG_SERTRAFO
  for(int t=0; t<2; t++)
    if(ths->dbg[t] != NULL)
      sertrafo_rmdbg(ths->dbg[t]);
#endif
  
  /* free memory */
  free(ths);
  /* ths=NULL; would be senseless, since we can not change the pointer itself */
}

#if PFFT_DEBUG_SERTRAFO
static sertrafo_dbg sertrafo_mkdbg(
    int is_fft, INT nb, int rnk, const INT *n, INT howmany,
    R *in, R *out,
    unsigned trafo_flag, unsigned transp_flag, unsigned fftw_flags, unsigned ioarray_flag,
    int dims_rnk, const X(iodim64) *dims, int howmany_rnk, const X(iodim64) *howmany_dims
    )
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  
  sertrafo_dbg ths = (sertrafo_dbg) malloc(sizeof(sertrafo_dbg_s));

  ths->is_fft = is_fft;
  ths->nb = nb;
  ths->rnk = rnk;
  if(rnk > 0){
    ths->n = (INT*) malloc(sizeof(INT)*rnk);
    for(int t=0; t<rnk; t++)
      ths->n[t] = n[t];
  } else
    ths->n = NULL;
  ths->howmany = howmany;
  ths->in = in;
  ths->out = out;
  
  ths->trafo_flag = trafo_flag;
  ths->transp_flag = transp_flag;
  ths->fftw_flags = fftw_flags;
  ths->ioarray_flag = ioarray_flag;
  
  ths->dims_rnk = dims_rnk;
  if(dims_rnk > 0){
    ths->dims = (X(iodim64)*) malloc(sizeof(X(iodim64)) * (size_t) dims_rnk);
    for(int t=0; t<dims_rnk; t++){
      ths->dims[t].n = dims[t].n;
      ths->dims[t].is = dims[t].is;
      ths->dims[t].os = dims[t].os;
    }
  } else
    ths->dims = NULL;

  ths->howmany_rnk = howmany_rnk;
  if(howmany_rnk > 0){
    ths->howmany_dims = (X(iodim64)*) malloc(sizeof(X(iodim64)) * (size_t) howmany_rnk);
    for(int t=0; t<howmany_rnk; t++){
      ths->howmany_dims[t].n = howmany_dims[t].n;
      ths->howmany_dims[t].is = howmany_dims[t].is;
      ths->howmany_dims[t].os = howmany_dims[t].os;
    }
  } else 
    ths->howmany_dims = NULL;

  return ths;
}

void sertrafo_rmdbg(
    sertrafo_dbg ths
    )
{
  /* plan was already destroyed or never initialized */
  if(ths==NULL)
    return;

  if(ths->n != NULL)
    free(ths->n);
  if(ths->dims != NULL)
    free(ths->dims);
  if(ths->howmany_dims != NULL)
    free(ths->howmany_dims);
  
  /* free memory */
  free(ths);
  /* ths=NULL; would be senseless, since we can not change the pointer itself */
}


static void print_vec(
    int rnk, INT *vec
    )
{
  fprintf(stderr, "[");
  for(int t=0; t<rnk-1; t++)
    fprintf(stderr, "%td, ", vec[t]);
  fprintf(stderr, "%td]", vec[rnk-1]);
}

static void print_dbg(
    sertrafo_dbg ths
    )
{
  fprintf(stderr, "PFFT_DBG_SERTRAFO: is_fft = %d, nb = %td, rnk = %d,  n = ", ths->is_fft, ths->nb, ths->rnk);
  print_vec(ths->rnk, ths->n);
  fprintf(stderr, ", howmany = %td, in = %p, out = %p\n", ths->howmany, ths->in, ths->out);
  
  fprintf(stderr, "PFFT_DBG_SERTRAFO: transp_flag = ");
  if(ths->transp_flag & PFFT_TRANSPOSED_IN){
    if(ths->transp_flag & PFFT_TRANSPOSED_OUT)
      fprintf(stderr, "PFFT_TRANSPOSED_IN | PFFT_TRANSPOSED_OUT");
    else
      fprintf(stderr, "PFFT_TRANSPOSED_IN");
  } else {
    if(ths->transp_flag & PFFT_TRANSPOSED_OUT)
      fprintf(stderr, "PFFT_TRANSPOSED_OUT");
    else
      fprintf(stderr, "PFFT_TRANSPOSED_NONE");
  } 
  fprintf(stderr, "\n");

  fprintf(stderr, "PFFT_DBG_SERTRAFO: trafo_flag = ");
  if(ths->trafo_flag & PFFTI_TRAFO_R2C){
    fprintf(stderr, "PFFTI_TRAFO_R2C");
  } else if(ths->trafo_flag & PFFTI_TRAFO_C2R){
    fprintf(stderr, "PFFTI_TRAFO_C2R");
  } else if(ths->trafo_flag & PFFTI_TRAFO_R2R){
    fprintf(stderr, "PFFTI_TRAFO_R2R");
  } else if(ths->trafo_flag & PFFTI_TRAFO_C2C){
    fprintf(stderr, "PFFTI_TRAFO_C2C");
  } else {
    fprintf(stderr, "!!! UNKNOWN TRAFO FLAG !!!");
  }
  fprintf(stderr, "\n");
 
  fprintf(stderr, "PFFT_DBG_SERTRAFO: ioarray_flag = ");
  if(ths->ioarray_flag & PFFTI_ARRAY_INPUT){
    fprintf(stderr, "PFFTI_ARRAY_INPUT");
  } else if(ths->ioarray_flag & PFFTI_ARRAY_OUTPUT){
    fprintf(stderr, "PFFTI_ARRAY_OUTPUT");
  } else if(ths->ioarray_flag & PFFTI_ARRAY_UNDEFINED){
    fprintf(stderr, "PFFTI_ARRAY_UNDEFINED");
  } else {
    fprintf(stderr, "!!! UNKNOWN ARRAY FLAG !!!");
  }
  fprintf(stderr, "\n");
  
  fprintf(stderr, "PFFT_DBG_SERTRAFO: dims_rnk = %d, ", ths->dims_rnk);
  for(int t=0; t<ths->dims_rnk; t++)
    fprintf(stderr, "dims[%d] = [n=%td, is=%td, os=%td], ", t, ths->dims[t].n, ths->dims[t].is, ths->dims[t].os);
  fprintf(stderr, "\n");
  
  fprintf(stderr, "PFFT_DBG_SERTRAFO: howmany_rnk = %d, ", ths->howmany_rnk);
  for(int t=0; t<ths->howmany_rnk; t++)
    fprintf(stderr, "howmany_dims[%d] = [n=%td, is=%td, os=%td], ", t, ths->howmany_dims[t].n, ths->howmany_dims[t].is, ths->howmany_dims[t].os);
  fprintf(stderr, "\n");
}

INT calc_n_total(
    sertrafo_dbg ths
    )
{
  INT tuple, n_total;
  
  if(ths==NULL)
    return 0;
  
  tuple = (ths->trafo_flag & PFFTI_TRAFO_C2C) ? 2 : 1;

  n_total = ths->howmany * ths->nb * tuple;
  for(int t=0; t<ths->rnk; t++)
    n_total *= ths->n[t];

  return n_total;
}
#endif



void PX(execute_sertrafo)(
    sertrafo_plan ths
    )
{
  if(ths==NULL)
    return;

#if PFFT_DEBUG_SERTRAFO
  static int counter=0;
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  R lsum, gsum;
  INT n_total;

  if(!counter)
    if(!myrank)
      fprintf(stderr, "PFFT_DBG_SERTRAFO: !!! Attention: checksums for out0 and in1 may differ between different runs, since the order of FFT and remap may change on all processes. Use FFTW_ESTIMATE to be sure that FFT is performed first on all processes. !!!\n");

  if(!myrank) fprintf(stderr, "\n");
  if(!myrank){
    if(ths->plan[0] != NULL){
      fprintf(stderr, "PFFT_DBG_SERTRAFO: counter = %d, plan0\n", counter);
      print_dbg(ths->dbg[0]);
    } else
      fprintf(stderr, "PFFT_DBG_SERTRAFO: nothing to do for FFTW plan0\n");
  }

  n_total = calc_n_total(ths->dbg[0]); 
 
  
  /* Checksum inputs */ 
  lsum=0.0;
  if(ths->plan[0] != NULL)
    for(INT k=0; k<n_total; k++)
      lsum += fabs(ths->dbg[0]->in[k]);
  MPI_Reduce(&lsum, &gsum, 1, PFFT_MPI_REAL_TYPE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(ths->plan[0] != NULL)
    if(!myrank) fprintf(stderr, "PFFT_DBG_SERTRAFO: counter = %d, Checksum(in0) = %e\n", counter, gsum);

// #if PFFT_DEBUG_GVARS 
//   if(counter==1){
//     int np;
//     int np0, np1, rnk0, rnk1;
//     MPI_Comm_size(MPI_COMM_WORLD, &np);
//     MPI_Comm_size(gdbg_comms_pm[0], &np0);
//     MPI_Comm_size(gdbg_comms_pm[1], &np1);
//     MPI_Comm_rank(gdbg_comms_pm[0], &rnk0);
//     MPI_Comm_rank(gdbg_comms_pm[1], &rnk1);
//     INT local_N[3], local_N_start[3];
// 
//     local_N[0] = 512/np0; local_N_start[0] = 512/np0 * rnk0;
//     local_N[1] = 512/np1; local_N_start[1] = 512/np1 * rnk1;
//     local_N[2] = 512; local_N_start[2] = 0;
//     
//     if(!myrank) fprintf(stderr, "!!! check all coefficients !!!\n");
//     if(!myrank) fprintf(stderr, "!!! local_N=[%td, %td, %td], local_N_start = [%td, %td, %td]\n",
//         local_N[0], local_N[1], local_N[2], local_N_start[0], local_N_start[1], local_N_start[2]);
// 
//     int lerr=0; 
//     INT m=0;
//     for(INT k2=local_N_start[2]; k2<local_N_start[2]+local_N[2]; k2++)
//       for(INT k0=local_N_start[0]; k0<local_N_start[0]+local_N[0]; k0++)
//         for(INT k1=local_N_start[1]; k1<local_N_start[1]+local_N[1]; k1++)
//           for(INT h=0; h<2; h++, m++){
//             R ind = h + 2*(k2 + 512*(k1 + 512*k0));
//             R data = ths->dbg[0]->in[m];
//             if( (data - ind) > 1e-13){
//               if(!lerr)
//                 if(!myrank)
//                   fprintf(stderr, "data[%td] = %e, ind = %e, k2=%td, k0=%td, k1=%td\n", data, m, ind, k2, k0, k1);
//               lerr = 1;
//             }
//           }
//     if(!myrank) fprintf(stderr, "coords = [%d, %d], lerr = %d\n", rnk0, rnk1, lerr);
//   }
// #endif
  
  /* Serial trafo */ 
  if(ths->plan[0] != NULL)
    X(execute)(ths->plan[0]);
  
  /* Checksum outputs */ 
  lsum=0.0;
  if(ths->plan[0] != NULL)
    for(INT k=0; k<n_total; k++)
      lsum += fabs(ths->dbg[0]->out[k]);
  MPI_Reduce(&lsum, &gsum, 1, PFFT_MPI_REAL_TYPE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(ths->plan[0] != NULL)
    if(!myrank) fprintf(stderr, "PFFT_DBG_SERTRAFO: counter = %d, Checksum(out0) = %e - Value may change\n", counter, gsum);

// #if PFFT_DEBUG_GVARS 
//   if(counter==1){
//     int np;
//     int np0, np1, rnk0, rnk1;
//     MPI_Comm_size(MPI_COMM_WORLD, &np);
//     MPI_Comm_size(gdbg_comms_pm[0], &np0);
//     MPI_Comm_size(gdbg_comms_pm[1], &np1);
//     MPI_Comm_rank(gdbg_comms_pm[0], &rnk0);
//     MPI_Comm_rank(gdbg_comms_pm[1], &rnk1);
//     INT local_N[3], local_N_start[3];
// 
//     local_N[0] = 512/np0; local_N_start[0] = 512/np0 * rnk0;
//     local_N[1] = 512/np1; local_N_start[1] = 512/np1 * rnk1;
//     local_N[2] = 512; local_N_start[2] = 0;
//     
//     if(!myrank) fprintf(stderr, "!!! check all coefficients !!!\n");
//     if(!myrank) fprintf(stderr, "!!! local_N=[%td, %td, %td], local_N_start = [%td, %td, %td]\n",
//         local_N[0], local_N[1], local_N[2], local_N_start[0], local_N_start[1], local_N_start[2]);
// 
//     int lerr=0; 
//     INT m=0;
//     for(INT k0=local_N_start[0]; k0<local_N_start[0]+local_N[0]; k0++)
//       for(INT k1=local_N_start[1]; k1<local_N_start[1]+local_N[1]; k1++)
//         for(INT k2=local_N_start[2]; k2<local_N_start[2]+local_N[2]; k2++)
//           for(INT h=0; h<2; h++, m++){
//             R ind = h + 2*(k2 + 512*(k1 + 512*k0));
//             R data = ths->dbg[0]->out[m];
//             if( (data - ind) > 1e-13){
//               if(!lerr)
//                 if(!myrank)
//                   fprintf(stderr, "data[%td] = %e, ind = %e, k0=%td, k1=%td, k2=%td\n", data, m, ind, k0, k1, k2);
//               lerr = 1;
//             }
//           }
//     fprintf(stderr, "coords = [%d, %d], lerr = %d\n", rnk0, rnk1, lerr);
//     
//     MPI_Barrier(MPI_COMM_WORLD);
//   }
// #endif
  
  if(!myrank) fprintf(stderr, "\n");
  if(!myrank){
    if(ths->plan[1] != NULL){
      fprintf(stderr, "PFFT_DBG_SERTRAFO: counter = %d, plan1\n", counter);
      print_dbg(ths->dbg[1]);
    } else
      fprintf(stderr, "PFFT_DBG_SERTRAFO: nothing to do for FFTW plan1\n");
  }
  
  n_total = calc_n_total(ths->dbg[1]); 
  
  /* Checksum inputs */ 
  lsum=0.0;
  if(ths->plan[1] != NULL)
    for(INT k=0; k<n_total; k++)
      lsum += fabs(ths->dbg[1]->in[k]);
  MPI_Reduce(&lsum, &gsum, 1, PFFT_MPI_REAL_TYPE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(ths->plan[1] != NULL)
    if(!myrank) fprintf(stderr, "PFFT_DBG_SERTRAFO: counter = %d, Checksum(in1) = %e - Value may change.\n", counter, gsum);
    
  /* Serial trafo */ 
  if(ths->plan[1] != NULL)
    X(execute)(ths->plan[1]);

  /* Checksum outputs */ 
  lsum=0.0;
  if(ths->plan[1] != NULL)
    for(INT k=0; k<n_total; k++)
      lsum += fabs(ths->dbg[1]->out[k]);
  MPI_Reduce(&lsum, &gsum, 1, PFFT_MPI_REAL_TYPE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(ths->plan[1] != NULL)
    if(!myrank) fprintf(stderr, "PFFT_DBG_SERTRAFO: counter = %d, Checksum(out1) = %e\n", counter, gsum);
  
  
  counter++;
#else
  /* execute all initialized serfft plans */
  for(int t=0; t<2; t++)
    if(ths->plan[t] != NULL)
      X(execute)(ths->plan[t]);
#endif
}









