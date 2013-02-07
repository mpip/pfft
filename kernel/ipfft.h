/*
 * Copyright (c) 2003, 2007-8 Matteo Frigo
 * Copyright (c) 2003, 2007-8 Massachusetts Institute of Technology
 * Copyright (c) 2002, 2009 Jens Keiner, Stefan Kunis, Daniel Potts
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

/* This file is based on ifftw3.h, ifftw3-mpi.h and infft.h */


/* PFFT internal header file */
#ifndef IPFFT_H
#define IPFFT_H 1

#include "config.h"

#define PFFT_DEBUG_SERTRAFO 0

/* debug execution of ghostcell send */
#define PFFT_DEBUG_GHOSTCELLS 0

/* debug execution of global transpositions */
#define PFFT_DEBUG_GTRANSP 0

/* use global variables to identify processes */
#define PFFT_DEBUG_GVARS 0

/******************************************************
 * Sometimes parallel FFTW 3.3.2 runs into a deadlock,
 * which can be avoided by calling fftw_forget_wisdom
 * before every parallel FFTW planer.
 * This issue was fixed in FFTW 3.3.3.
 ******************************************************/
#define PFFT_BUGFIX_FORGET_PARALLEL_FFTW_WISDOM 0

/* Begin: This part is based on ifftw3.h */
#include <stdlib.h>		/* size_t */
#include <stdarg.h>		/* va_list */
#include <stddef.h>		/* ptrdiff_t */
#include <stdio.h>		/* fprintf */

#include <math.h>
#include <mpi.h>
#include <fftw3-mpi.h>

#define IPFFT_EXTERN extern

typedef ptrdiff_t INT;
/*
  integral type large enough to contain a stride (what ``int'' should
  have been in the first place.
*/

#define IF(x,a,b) ((x)?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define ABS(x) (((x)>0)?(x):(-(x)))
#define SIGN(a) (((a)>=0)?1:-1)
#define MACRO_PLAIN_INDEX_3D_VECTOR(k, n)     ( k[2] + n[2]*(k[1] + n[1]*k[0]) )
#define MACRO_PLAIN_INDEX_3D(k0,k1,k2,n)      ( (k2) + n[2]*( (k1) + n[1]*(k0) ) )

/* provide FFTW_MANGLE_PREFIX macro */
#ifndef FFTW_MANGLE_PREFIX
# define FFTW_MANGLE_PREFIX(name)  name
#endif

/* all pfft identifiers start with pfft (or pfftf etc.) */
#define CONCAT(prefix, name) prefix ## name
#if defined(PFFT_SINGLE)
  typedef float R;
  typedef FFTW_MANGLE_PREFIX(fftwf_complex) C;
#  define PFFT_MPI_REAL_TYPE MPI_FLOAT
#  define PX(name) CONCAT(pfftf_, name)
#  define X(name)  FFTW_MANGLE_FLOAT(name)
#  define PFFT_MATH(name) CONCAT(name, f)
#elif defined(PFFT_LDOUBLE)
  typedef long double R;
  typedef FFTW_MANGLE_PREFIX(fftwl_complex) C;
#  define PFFT_MPI_REAL_TYPE MPI_LONG_DOUBLE
#  define PX(name) CONCAT(pfftl_, name)
#  define X(name)  FFTW_MANGLE_LONG_DOUBLE(name)
#  define PFFT_MATH(name) CONCAT(name, l)
#else
  typedef double R;
  typedef FFTW_MANGLE_PREFIX(fftw_complex) C;
#  define PFFT_MPI_REAL_TYPE MPI_DOUBLE
#  define PX(name) CONCAT(pfft_, name)
#  define X(name)  FFTW_MANGLE_DOUBLE(name)
#  define PFFT_MATH(name) name
#endif


/* macros for mathematical functions corresponding to the PFFT float data type */
#define pfft_pow(_x_, _y_)   PFFT_MATH(pow)(_x_, _y_)
#define pfft_creal(_x_) PFFT_MATH(creal)(_x_)
#define pfft_cimag(_x_) PFFT_MATH(cimag)(_x_)
#define pfft_cabs(_x_)  PFFT_MATH(cabs)(_x_)
#define pfft_log2(_x_)  PFFT_MATH(log2)(_x_)
#define pfft_sqrt(_x_)  PFFT_MATH(sqrt)(_x_)

#define XM(name)  X(CONCAT(mpi_, name))

#define PFFTI_GC_BORDERS          (0U)
#define PFFTI_GC_CORNERS          (1U<< 0)
#define PFFTI_GC_TRAFO            (0U)
#define PFFTI_GC_ADJOINT          (1U<< 1)

#define PFFTI_PRINT_TIMER_BASIC   (1U<<0)
#define PFFTI_PRINT_TIMER_ADV     (1U<<1)

#define FFT_DATA_INIT(i,j,k,t,n0,n1,n2) \
 ((R)1000 / ((R) ( (i)*(n1)*(n2)*2 + (j)*(n2)*2 + (k)*2 + (t) + 1 )))
  
//#define FFT_DATA_INIT(i,j,k,t,n0,n1,n2) 
// ((R)1000 / (((R)i) + ((R)j)*((R)n0) + ((R)k)*((R)n0)*((R)n1) + ((R)t)*((R)n0)*((R)n1)*((R)n2) + 1) )
// #define FFT_DATA_INIT  (i + j*ni[0] + k*ni[0]*ni[1] + 1)
// #define FFT_DATA_INIT (10*(((R)rand())/RAND_MAX - 0.5))


/* internal flags */
#define PFFTI_TRAFO_C2C            (1U<< 0)
#define PFFTI_TRAFO_R2C            (1U<< 1)
#define PFFTI_TRAFO_C2R            (1U<< 2)
#define PFFTI_TRAFO_RE00           (1U<< 3)
#define PFFTI_TRAFO_RE01           (1U<< 4)
#define PFFTI_TRAFO_RE10           (1U<< 5)
#define PFFTI_TRAFO_RE11           (1U<< 6)
#define PFFTI_TRAFO_RO00           (1U<< 7)
#define PFFTI_TRAFO_RO01           (1U<< 8)
#define PFFTI_TRAFO_RO10           (1U<< 9)
#define PFFTI_TRAFO_RO11           (1U<<10)
#define PFFTI_TRAFO_RC2C           (1U<<11)
#define PFFTI_TRAFO_SKIP           (1U<<12)
#define PFFTI_TRAFO_PHANTOM        (1U<<13)

#define PFFTI_TRAFO_R2R   \
  (PFFTI_TRAFO_REDFT | PFFTI_TRAFO_RODFT)
#define PFFTI_TRAFO_REDFT \
  (PFFTI_TRAFO_RE00 | PFFTI_TRAFO_RE01 | PFFTI_TRAFO_RE10 | PFFTI_TRAFO_RE11)
#define PFFTI_TRAFO_RODFT \
  (PFFTI_TRAFO_RO00 | PFFTI_TRAFO_RO01 | PFFTI_TRAFO_RO10 | PFFTI_TRAFO_RO11)
#define PFFTI_TRAFO_RDFT  \
  (PFFTI_TRAFO_R2C | PFFTI_TRAFO_C2R)

#define PFFTI_OUSAM_EMBED          (1U<< 0)
#define PFFTI_OUSAM_TRUNC          (1U<< 1)
#define PFFTI_OUSAM_TRANSPOSED     (1U<< 2)

#define PFFTI_BUFFERED_INPLACE     (1U<< 0) /* use second array of same size, similar to out-of-place but results end up in input array */


/* For r2c, c2r we sometimes need to know wether we
 * operate on the input or output array */
#define PFFTI_ARRAY_UNDEFINED      (1U<< 0)
#define PFFTI_ARRAY_INPUT          (1U<< 1)
#define PFFTI_ARRAY_OUTPUT         (1U<< 2)





#ifndef PFFT_H

typedef struct PX(plan_s) *PX(plan);
typedef struct PX(gcplan_s) *PX(gcplan);

typedef struct {
  int rnk_pm;
  int rnk_trafo;
  int rnk_remap;
  int iter;
  double whole;
  double *trafo;
  double *remap;
  double remap_3dto2d[2];
} PX(timer_s);
typedef PX(timer_s) *PX(timer);

typedef struct {
  int iter;
  double whole;
  double pad_zeros;
  double exchange;
} PX(gctimer_s);
typedef PX(gctimer_s) *PX(gctimer);

#endif /* !PFFT_H */

/* plan for debug infos of serial trafo (c2c, r2c, c2r, r2r) */
#if PFFT_DEBUG_SERTRAFO
typedef struct{
  int is_fft;
  INT nb;
  int rnk;
  INT *n;
  INT howmany;
  R *in;
  R *out;
  
  unsigned trafo_flag;
  unsigned transp_flag;
  unsigned fftw_flags;
  unsigned ioarray_flag;
  
  int dims_rnk;
  X(iodim64) *dims;
  int howmany_rnk;
  X(iodim64) *howmany_dims;
} sertrafo_dbg_s;
typedef sertrafo_dbg_s *sertrafo_dbg;
#endif

/* plan for serial trafo (c2c, r2c, c2r, r2r) */
typedef struct{
  X(plan) plan[2];
#if PFFT_DEBUG_SERTRAFO
  sertrafo_dbg dbg[2];
#endif
} sertrafo_plan_s;
typedef sertrafo_plan_s *sertrafo_plan;


/* plan for debug infos of global transposition */
#if PFFT_DEBUG_GTRANSP
typedef struct{
  INT N0;
  INT N1;
  INT hm;
  INT blk0;
  INT blk1;
  R *in;
  R *out;
  MPI_Comm comm;
  unsigned fftw_flags;

  INT mem;
  INT local_N0;
  INT local_N0_start;
  INT local_N1;
  INT local_N1_start;
} gtransp_dbg_s;
typedef gtransp_dbg_s *gtransp_dbg;
#endif

/* We need the following wrapper to store debug infos */
/* plan for global transposition */
typedef struct{
  X(plan) plan;
#if PFFT_DEBUG_GTRANSP
  gtransp_dbg dbg;
#endif
} gtransp_plan_s;
typedef gtransp_plan_s *gtransp_plan;



/* plan for 1d over-/undersampling */
typedef struct{
  INT N0;   /* number of loops (first array dimension) */
  INT N1i;  /* 2nd array size before oversampling */
  INT N1o;  /* 2nd array size after oversampling */
  R *in;    /* pointer to input */
  R *out;   /* pointer to output */
  INT Cl;   /* 1st block size of coefficients */
  INT Cm;   /* gap size: number of zeros to add / coefficients to skip */
  INT Cr;   /* 2nd block size of coefficients */
  INT Pi;   /* padding for r2c, c2r before oversampling */
  INT Po;   /* padding for r2c, c2r after oversampling */ 
  unsigned ousam_flag; /* set to PFFTI_EMBED or PFFTI_TRUNC */
} ousam_plan_1d_s;
typedef ousam_plan_1d_s *ousam_plan_1d;

/* plan for d-dim. over-/undersampling */
typedef struct{
  int rnk;
  ousam_plan_1d *ousam_1d;
} ousam_plan_dd_s;
typedef ousam_plan_dd_s *ousam_plan_dd;

/* plan for over-/undersampled serial trafo */
typedef struct{
  ousam_plan_dd embed; /* multidimensional embed */ 
  ousam_plan_dd trunc; /* multidimensional trunc */
  sertrafo_plan trafo; /* multidimensional trafo */
} outrafo_plan_s;
typedef outrafo_plan_s *outrafo_plan;

/* plan for internal remap of 3d to 2d data
 * decomposition for 3d trafos */
typedef struct{
  gtransp_plan global_remap[2];
  sertrafo_plan local_transp[2];
} remap_3dto2d_plan_s;
typedef remap_3dto2d_plan_s *remap_3dto2d_plan;


/* plan for parallel trafo */
struct PX(plan_s){
  int rnk_n;
  INT *n;
  INT *ni;
  INT *no;
  INT howmany;
  INT *iblock;
  INT *mblock;
  INT *oblock;
  MPI_Comm comm_cart;
  int rnk_pm;
  MPI_Comm *comms_pm;
  int *np;
  R *in;
  R *out;
  int sign;
  X(r2r_kind) *kinds;
  unsigned fftw_flags;
  unsigned transp_flag;
  unsigned trafo_flag;
  unsigned opt_flag;

  /* save the init flags for later reference */
  unsigned pfft_flags;

  gtransp_plan *global_remap;
  outrafo_plan *serial_trafo;
  remap_3dto2d_plan remap_3dto2d[2];

  PX(timer) timer;
};
typedef struct PX(plan_s) plan_s;

struct PX(gcplan_s){
  int rnk_n;
  int *np;
  INT *n;
  INT *loc_n;
  INT *gc_below;
  INT *gc_above;
  INT *ngc;
  INT *ngc_prec;
  INT *ngc_succ;
  int *rnk_prec;
  int *rnk_succ;
  INT *blk;
  INT tuple;
  R *data;
  MPI_Comm *comms_pm;
  MPI_Group grp;
  MPI_Win win;
  MPI_Comm comm_cart;
  unsigned alg_flag;
  PX(gctimer) timer_exg;
  PX(gctimer) timer_red;
};
typedef struct PX(gcplan_s) gcplan_s;


/* block.c */

INT PX(global_block_size)(
    INT arraySize, INT userBlockSize, int numProcs);
void PX(local_block_size_and_offset)(
    INT global_array_size, INT global_block_size, int which_block,
    INT *local_block_size, INT *local_block_start);
INT PX(local_block_size_shifted)(
    INT globalArraySize, INT userBlockSize, int shift, MPI_Comm commCart1d);
INT PX(local_block_offset)(
    INT arrayLength, INT globalBlockSize, int whichBlock);
INT PX(local_block_size)(
    INT arraySize, INT globalBlockSize, int whichBlock);
INT PX(num_blocks)(INT global_array_size, INT global_block_size);

/* malloc.c */

void PX(die)(
    const char *s, MPI_Comm comm);

/* partrafo.c */

INT PX(local_size_partrafo)(
    int rnk_n, const INT *n, const INT *ni, const INT *no,
    INT howmany, const INT *iblock, const INT *oblock,
    MPI_Comm comm_cart,
    unsigned trafo_flag_user, unsigned pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start);
PX(plan) PX(plan_partrafo)(
    int rnk_n, const INT *n, const INT *ni, const INT *no,
    INT howmany, const INT *iblock_user, const INT *oblock_user,
    R *in, R *out, MPI_Comm comm_cart,
    int sign, const X(r2r_kind) *kinds, const int *skip_trafos_user,
    unsigned trafo_flag, unsigned pfft_flags);
void PX(rmplan)(
    PX(plan) ths);

/* partrafo-transposed.c */

INT  PX(local_size_partrafo_transposed)(
    int rnk_n, const INT *n, const INT *ni, const INT *no,
    INT howmany, const INT *iblock, const INT *oblock,
    int rnk_pm, MPI_Comm *comms_pm,
    unsigned transp_flag, const unsigned *trafo_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start);
void PX(plan_partrafo_transposed)(
    int rnk_n, const INT *n, const INT *ni, const INT *no,
    INT howmany, const INT *iblock, const INT *oblock,
    int rnk_pm, MPI_Comm *comms_pm,
    R *in, R *out, int sign, const X(r2r_kind) *kinds,
    unsigned transp_flag, const unsigned *trafo_flags,
    unsigned opt_flag, unsigned io_flag, unsigned fftw_flags, 
    outrafo_plan *trafos, gtransp_plan *remaps);

/* transpose.c */

void PX(get_global_transp_param)(
    int step, int rnk_pm, const INT *ni, const INT *no,
    const INT *local_ni, const INT *local_no,
    INT tuple_size, const INT *iblock, const INT *oblock,
    unsigned trafo_flag, unsigned transp_flag,
    INT *N0, INT *N1, INT *h0, INT *h1,
    INT *hm, INT *blk0, INT *blk1);
INT PX(local_size_global_transp)(
    INT N0, INT N1, INT h0, INT h1, INT hm, INT blk0, INT blk1,
    MPI_Comm comm );
gtransp_plan PX(plan_global_transp)(
    INT N0, INT N1, INT h0, INT h1, INT hm, INT blk0, INT blk1,
    MPI_Comm comm, R *in, R *out,
    unsigned transp_flag, unsigned fftw_flags);
void PX(gtransp_rmplan)(
    gtransp_plan ths);
void PX(execute_gtransp)(
    gtransp_plan ths);

/* outrafo.c */


void PX(get_outrafo_param)(
    int step, int rnk_pm,
    const INT *n, const INT *ni, const INT *no,
    const INT *local_ni, const INT *local_no,
    const X(r2r_kind) *kinds, unsigned transp_flag, const unsigned *trafo_flags,
    INT *Nb, INT *N, INT *Ni, INT *No,
    X(r2r_kind) *kind, unsigned *trafo_flag);
INT PX(local_size_outrafo)(
    INT nb, int rnk, const INT *n, const INT *ni, const INT *no, INT howmany,
    unsigned trafo_flag);
outrafo_plan PX(plan_outrafo)(
    INT nb, int rnk, const INT *n, const INT *ni, const INT *no, INT howmany,
    R *in, R *out, int sign, const X(r2r_kind) *kinds,
    unsigned trafo_flag, unsigned transp_flag,
    unsigned opt_flag, unsigned fftw_flags);
void PX(execute_outrafo)(
    outrafo_plan ths);
void PX(outrafo_rmplan)(
    outrafo_plan ths);

/* ousample.c */

INT PX(local_size_ousam_dd)(
    INT nb, int rnk, const INT *ni, const INT *no, INT howmany, 
    unsigned trafo_flag);
ousam_plan_dd PX(plan_ousam_dd)(
    INT nb, int rnk, const INT *ni, const INT *no, INT howmany, 
    R *in, R *out, unsigned trafo_flag, unsigned ousam_flag);
void PX(execute_ousam_dd)(
    ousam_plan_dd ths);
void PX(ousam_dd_rmplan)(
    ousam_plan_dd ths);

/* sertrafo.c */

INT PX(local_size_sertrafo)(
    INT nb, int rnk, const INT *n, INT howmany,
    unsigned trafo_flag);
sertrafo_plan PX(plan_sertrafo)(
    INT nb, int rnk, const INT *n, INT howmany,
    R *in, R *out, int sign, const X(r2r_kind) *kinds,
    unsigned trafo_flag, unsigned transp_flag, unsigned io_flag,
    unsigned opt_flag, unsigned fftw_flags);
void PX(execute_sertrafo)(
    sertrafo_plan ths);
void PX(sertrafo_rmplan)(
    sertrafo_plan ths);

/* procmesh.c */


int PX(is_cart_procmesh)(
    MPI_Comm comm_cart);
int PX(is_cart_procmesh_2d)(
    MPI_Comm comm_cart_2d);
void PX(split_cart_procmesh)(
    MPI_Comm comm_cart,
    MPI_Comm *comms_1d);
void PX(split_cart_procmesh_3dto2d_p0q0)(
    MPI_Comm comm_cart_3d,
    MPI_Comm *comm_1d);
void PX(split_cart_procmesh_3dto2d_p1q1)(
    MPI_Comm comm_cart_3d,
    MPI_Comm *comm_1d);
void PX(get_procmesh_dims_2d)(
    MPI_Comm comm_cart_3d,
    int *p0, int *p1, int *q0, int *q1);
void PX(split_cart_procmesh_for_3dto2d_remap_q0)(
    MPI_Comm comm_cart_3d,
    MPI_Comm *comm_1d);
void PX(split_cart_procmesh_for_3dto2d_remap_q1)(
    MPI_Comm comm_cart_3d,
    MPI_Comm *comms_1d);
int PX(get_mpi_cart_coord_1d)(MPI_Comm comm_cart_1d, int *coord);
int PX(get_mpi_cart_coords)(MPI_Comm comm_cart, int maxdims, int *coords);
int PX(get_mpi_cart_dims)(MPI_Comm comm_cart, int maxdims, int *dims);


/* remap_3dto2d.c */

int PX(local_size_remap_3dto2d_transposed)(
    int rnk_n, const INT *n, INT howmany, 
    MPI_Comm comm_cart_3d, 
    unsigned transp_flag, unsigned trafo_flag,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start);
remap_3dto2d_plan PX(plan_remap_3dto2d_transposed)(
    int rnk_n, const INT *n, INT howmany, 
    MPI_Comm comm_cart_3d, R *in, R *out, 
    unsigned transp_flag, unsigned trafo_flag,
    unsigned opt_flag, unsigned io_flag, unsigned fftw_flags);
void PX(execute_remap_3dto2d)(
    remap_3dto2d_plan ths);
void PX(remap_3dto2d_rmplan)(
    remap_3dto2d_plan ths);
void PX(default_block_size_3dto2d)(
    const INT *n, int p0, int p1, int q0, int q1,
    INT *iblk, INT *mblk, INT *oblk);










/* timer.c */

PX(timer) PX(mktimer)(
    int rnk_pm);
/* use public api to destroy timer with PX(destroy_timer) */





/* gcells_plan.c */
 
INT PX(local_size_gc_internal)(
    int rnk_n, const INT *local_n, const INT *local_start,
    INT tuple, const INT *gc_below_user, const INT *gc_above_user,
    INT *local_ngc, INT *local_gc_start);

PX(gcplan) PX(plan_rgc_internal)(
    int rnk_n, const INT *n, INT tuple_size, const INT *block_user,
    const INT *gc_below_user, const INT *gc_above_user,
    R *data,int rnk_pm, MPI_Comm *comms_pm, MPI_Comm comm_cart,
    unsigned gc_flag);
void PX(exchange_gc)(
    PX(gcplan) ths);
void PX(reduce_gc)(
    PX(gcplan) ths);
void PX(rmplan_gc)(
    PX(gcplan) ths);

/* gcells_sendrecv.c */

void PX(exchange_gc_sendrecv)(
    PX(gcplan) ths);
void PX(reduce_gc_sendrecv)(
    PX(gcplan) ths);

/* gcells_RMA.c */

void PX(exchange_gc_RMA)(
    PX(gcplan) ths);
void PX(reduce_gc_RMA)(
    PX(gcplan) ths);

/* gctimer.c */

PX(gctimer) PX(gc_mktimer)(void);
 

    
    
// /* check.c */
// 
// int PX(check_cart_2d_internal)(
//     MPI_Comm comm_cart_2d, int verbose);
// 
 


#endif /* !IPFFT_H */
