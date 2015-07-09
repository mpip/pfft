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

#include <complex.h>
#include "pfft.h"
#include "ipfft.h"
#include "util.h"


/* Global infos about procmesh are only enabled in debug mode
 * Otherwise we do not use any global variables. */
#if PFFT_DEBUG_GVARS
  MPI_Comm *gdbg_comm_cart;
  int gdbg_rnk_pm=-1;
  MPI_Comm *gdbg_comms_pm;
#endif

static void twiddle_input(
    const PX(plan) ths,
    R *planned_in, R *planned_out,
    R *executed_in, R *executed_out);
static void twiddle_output(
    const PX(plan) ths,
    R *planned_in, R *planned_out,
    R *executed_in, R *executed_out);

static void execute_transposed(
    int rnk_pm, outrafo_plan *trafos, gtransp_plan *remaps,
    double *timer_trafo, double *timer_remap, R * plannedin, R * plannedout, 
    R * in, R * out, MPI_Comm comm_cart);


static INT plain_index(
    int rnk, const INT *n, const INT *kvec);
static void vector_index(
    int rnk, const INT *n, INT k,
    INT *kvec);

static void complex_conjugate(
    const R* in, R* out, const int rnk_n, const INT *local_n);


static C init_scalar(
    int rnk_n, const INT *n, const INT *i
    )
{
  R m = (R) plain_index(rnk_n, n, i);

  /* avoid division by zero */
  if(m==0) return 1500.0 + 1250.0 * _Complex_I;

  return 1000.0/(2*m) + 1000.0/(2*m+1) * _Complex_I;
}

static C init_scalar_periodic(
    int rnk_n, const INT *n, const INT *ind
    )
{
  INT *periodic_ind = PX(malloc_INT)(rnk_n);

  /* assure periodicity in all directions */
  for(int t=0; t<rnk_n; t++)
    periodic_ind[t] = ind[t] % n[t];

  C result = init_scalar(rnk_n, n, periodic_ind);

  free(periodic_ind);

  return result;
}

static void init_array(
    int rnk_n, const INT *n, const INT *local_n, const INT *local_start,
    unsigned arraytype,
    void *data
    )
{
  INT ln_tot;
  INT *kvec_loc, *kvec_glob, *kvec_glob_mirrored;

  kvec_loc  = PX(malloc_INT)(rnk_n);
  kvec_glob = PX(malloc_INT)(rnk_n);
  kvec_glob_mirrored = PX(malloc_INT)(rnk_n);
 
  ln_tot = PX(prod_INT)(rnk_n, local_n);
  for(INT k=0; k<ln_tot; k++){
    vector_index(rnk_n, local_n, k, kvec_loc);
    PX(vadd_INT)(rnk_n, kvec_loc, local_start, kvec_glob);
    PX(vsub_INT)(rnk_n, n, kvec_glob, kvec_glob_mirrored);
    C d1 = init_scalar_periodic(rnk_n, n, kvec_glob);
    C d2 = init_scalar_periodic(rnk_n, n, kvec_glob_mirrored);
    switch (arraytype){
      case PFFTI_ARRAYTYPE_REAL: 
        /* set padding element to zero */
        ((R*)data)[k] = (kvec_glob[rnk_n-1] < n[rnk_n-1]) ? (R) d1 : 0.0; break;
      case PFFTI_ARRAYTYPE_COMPLEX:
        ((C*)data)[k] = d1; break;
      case PFFTI_ARRAYTYPE_HERMITIAN_COMPLEX:
        ((C*)data)[k] = 0.5 * (d1 + conj(d2)); break;
    }
  }

  free(kvec_loc); free(kvec_glob); free(kvec_glob_mirrored);
}

static void clear_array(
    int rnk_n, const INT *n, const INT *local_n, const INT *local_start,
    unsigned arraytype,
    void *data
    )
{
  INT ln_tot;
 
  ln_tot = PX(prod_INT)(rnk_n, local_n);
  for(INT k=0; k<ln_tot; k++){
    switch (arraytype){
      case PFFTI_ARRAYTYPE_REAL: 
        /* set padding element to zero */
        ((R*)data)[k] = 0;
        break;
      case PFFTI_ARRAYTYPE_COMPLEX:
        ((C*)data)[k] = 0;
        break;
      case PFFTI_ARRAYTYPE_HERMITIAN_COMPLEX:
        ((C*)data)[k] = 0;
        break;
    }
  }
}

static R check_array(
    int rnk_n, const INT *n, const INT *local_n, const INT *local_start,
    unsigned arraytype,
    const void *data, MPI_Comm comm
    )
{
  INT ln_tot;
  INT *kvec_loc, *kvec_glob, *kvec_glob_mirrored;
  R err, maxerr, globmaxerr;

  err = maxerr = 0;

  kvec_loc  = PX(malloc_INT)(rnk_n);
  kvec_glob = PX(malloc_INT)(rnk_n);
  kvec_glob_mirrored = PX(malloc_INT)(rnk_n);
 
  ln_tot = PX(prod_INT)(rnk_n, local_n);
  for(INT k=0; k<ln_tot; k++){
    vector_index(rnk_n, local_n, k, kvec_loc);
    PX(vadd_INT)(rnk_n, kvec_loc, local_start, kvec_glob);
    PX(vsub_INT)(rnk_n, n, kvec_glob, kvec_glob_mirrored);
    C d1 = init_scalar_periodic(rnk_n, n, kvec_glob);
    C d2 = init_scalar_periodic(rnk_n, n, kvec_glob_mirrored);
    switch (arraytype){
      case PFFTI_ARRAYTYPE_REAL: 
        /* ignore padding elements */
        err = (kvec_glob[rnk_n-1] < n[rnk_n-1]) ? cabs( ((R*)data)[k] - (R) d1 ) : 0.0; break;
      case PFFTI_ARRAYTYPE_COMPLEX:
        err = cabs( ((C*)data)[k] - d1); break;
      case PFFTI_ARRAYTYPE_HERMITIAN_COMPLEX:
        err = cabs( ((C*)data)[k] - 0.5 * (d1 + conj(d2))); break;
    }

    if( err > maxerr )
      maxerr = err;
  }

  free(kvec_loc); free(kvec_glob); free(kvec_glob_mirrored);

  MPI_Allreduce(&maxerr, &globmaxerr, 1, PFFT_MPI_REAL_TYPE, MPI_MAX, comm);
  return globmaxerr;
}



#ifdef _OPENMP
#include <omp.h>
static int _nthreads;
#endif

void PX(plan_with_nthreads) (int nthreads){
#ifdef _OPENMP
  _nthreads = nthreads;
  omp_set_num_threads(nthreads);
  X(plan_with_nthreads)(nthreads);
#endif
}

int PX(get_nthreads)() {
#ifdef _OPENMP
  return _nthreads;
#else
  return 1;
#endif
}

/* wrappers for fftw init and cleanup */
void PX(init) (void){
#ifdef _OPENMP
  X(init_threads)();
  PX(plan_with_nthreads)(omp_get_max_threads());
#endif
  XM(init)();
}

void PX(cleanup) (void){
#ifdef _OPENMP
  X(cleanup_threads());
#endif
  XM(cleanup)();
}


void PX(destroy_gcplan) (
    PX(gcplan) ths
    )
{
  PX(rmplan_gc)(ths);
}


void PX(destroy_plan)(
    PX(plan) ths
    )
{
  if(ths==NULL){
    PX(fprintf)(MPI_COMM_WORLD, stderr, "!!! Error: Can not finalize PFFT Plan == NULL !!!\n");
    return;
  }

  PX(rmplan)(ths);
}




static INT plain_index(
    int rnk, const INT *n, const INT *kvec
    )
{
  INT k=0;

  for(INT t=0; t<rnk; t++)
    k += k*n[t] + kvec[t];

  return k;
}

  
static void vector_index(
    int rnk, const INT *n, INT k,
    INT *kvec
    )
{
  for(INT t=rnk-1; t>=0; t--){
    kvec[t] = k%n[t];
    k = (k - kvec[t])/n[t];
  }
}


static void complex_conjugate(
    const R* in, R* out, const int rnk_n, const INT *local_n
    )
{
  if (in == NULL)
    return;

  INT local_n_total = PX(prod_INT)(rnk_n, local_n);

  for (INT k=0; k<local_n_total; k++) {
    out[2*k] = in[2*k];
    out[2*k+1] = -in[2*k+1];
  }
}


void PX(init_input_complex_3d)(
    const INT *n, const INT *local_n, const INT *local_start,
    C *data
    )
{
  int rnk_n=3;

  PX(init_input_complex)(rnk_n, n, local_n, local_start,
      data);
}

void PX(init_input_complex)(
    int rnk_n, const INT *n, const INT *local_n, const INT *local_start,
    C *data
    )
{ 
  init_array(rnk_n, n, local_n, local_start, PFFTI_ARRAYTYPE_COMPLEX,
      data);
}

void PX(init_input_complex_hermitian_3d)(
    const INT *n, const INT *local_n, const INT *local_start,
    C *data
    )
{ 
  int rnk_n=3;

  PX(init_input_complex_hermitian)(rnk_n, n, local_n, local_start,
      data);
}

void PX(init_input_complex_hermitian)(
    int rnk_n, const INT *n, const INT *local_n, const INT *local_start,
    C *data
    )
{ 
  init_array(rnk_n, n, local_n, local_start, PFFTI_ARRAYTYPE_HERMITIAN_COMPLEX,
      data);
}


void PX(init_input_real_3d)(
    const INT *n, const INT *local_n, const INT *local_start,
    R *data
    )
{
  int rnk_n=3;

  PX(init_input_real)(rnk_n, n, local_n, local_start,
      data);
}


void PX(init_input_real)(
    int rnk_n, const INT *n, const INT *local_n, const INT *local_start,
    R *data
    )
{ 
  init_array(rnk_n, n, local_n, local_start, PFFTI_ARRAYTYPE_REAL,
      data);
}

/* clear input array; setting it to zero */
void PX(clear_input_complex_3d)(
    const INT *n, const INT *local_n, const INT *local_start,
    C *data
    )
{
  int rnk_n=3;

  PX(clear_input_complex)(rnk_n, n, local_n, local_start,
      data);
}

void PX(clear_input_complex)(
    int rnk_n, const INT *n, const INT *local_n, const INT *local_start,
    C *data
    )
{ 
  clear_array(rnk_n, n, local_n, local_start, PFFTI_ARRAYTYPE_COMPLEX,
      data);
}

void PX(clear_input_complex_hermitian_3d)(
    const INT *n, const INT *local_n, const INT *local_start,
    C *data
    )
{ 
  int rnk_n=3;

  PX(clear_input_complex_hermitian)(rnk_n, n, local_n, local_start,
      data);
}

void PX(clear_input_complex_hermitian)(
    int rnk_n, const INT *n, const INT *local_n, const INT *local_start,
    C *data
    )
{ 
  clear_array(rnk_n, n, local_n, local_start, PFFTI_ARRAYTYPE_HERMITIAN_COMPLEX,
      data);
}


void PX(clear_input_real_3d)(
    const INT *n, const INT *local_n, const INT *local_start,
    R *data
    )
{
  int rnk_n=3;

  PX(clear_input_real)(rnk_n, n, local_n, local_start,
      data);
}


void PX(clear_input_real)(
    int rnk_n, const INT *n, const INT *local_n, const INT *local_start,
    R *data
    )
{ 
  clear_array(rnk_n, n, local_n, local_start, PFFTI_ARRAYTYPE_REAL,
      data);
}



/* Check results after one forward and backward FFT.
 * Only works for intial data layout */
R PX(check_output_complex_3d)(
    const INT *n, const INT *local_n, const INT *local_start,
    const C *data, MPI_Comm comm
    )
{
  int rnk_n = 3; 

  return PX(check_output_complex)(rnk_n, n, local_n, local_start,
      data, comm);
}

R PX(check_output_complex)(
    int rnk_n, const INT *n, const INT *local_n, const INT *local_start,
    const C *data, MPI_Comm comm
    )
{ 
  return check_array(rnk_n, n, local_n, local_start, PFFTI_ARRAYTYPE_COMPLEX, data, comm);
}

R PX(check_output_complex_hermitian_3d)(
    const INT *n, const INT *local_n, const INT *local_start,
    const C *data, MPI_Comm comm
    )
{
  int rnk_n = 3; 

  return PX(check_output_complex_hermitian)(rnk_n, n, local_n, local_start,
      data, comm);
}

R PX(check_output_complex_hermitian)(
    int rnk_n, const INT *n, const INT *local_n, const INT *local_start,
    const C *data, MPI_Comm comm
    )
{ 
  return check_array(rnk_n, n, local_n, local_start, PFFTI_ARRAYTYPE_HERMITIAN_COMPLEX, data, comm);
}

R PX(check_output_real_3d)(
    const INT *n, const INT *local_n, const INT *local_start,
    const R *data, MPI_Comm comm
    )
{
  int rnk_n = 3; 

  return PX(check_output_real)(rnk_n, n, local_n, local_start,
      data, comm);
}

R PX(check_output_real)(
    int rnk_n, const INT *n, const INT *local_n, const INT *local_start,
    const R *data, MPI_Comm comm
    )
{ 
  return check_array(rnk_n, n, local_n, local_start, PFFTI_ARRAYTYPE_REAL, data, comm);
}


INT PX(local_size_gc_3d)(
    const INT *local_n, const INT *local_start,
    const INT *gc_below, const INT *gc_above,
    INT *local_ngc, INT *local_gc_start
    )
{
  int rnk_n = 3;
  
  return PX(local_size_gc)(
      rnk_n, local_n, local_start, gc_below, gc_above,
      local_ngc, local_gc_start);
}

INT PX(local_size_gc)(
    int rnk_n, const INT *local_n, const INT *local_start,
    const INT *gc_below, const INT *gc_above,
    INT *local_ngc, INT *local_gc_start
    )
{
  INT howmany = 1;
  
  return PX(local_size_many_gc)(
      rnk_n, local_n, local_start, howmany, gc_below, gc_above,
      local_ngc, local_gc_start);
}



PX(gcplan) PX(plan_rgc_3d)(
    const INT *n, const INT *gc_below, const INT *gc_above,
    R *data, MPI_Comm comm_cart, unsigned gc_flags
    )
{
  int rnk_n = 3;
  
  return PX(plan_rgc)(rnk_n, n, gc_below, gc_above,
      data, comm_cart, gc_flags);
}

PX(gcplan) PX(plan_cgc_3d)(
    const INT *n, const INT *gc_below, const INT *gc_above,
    C *data, MPI_Comm comm_cart, unsigned gc_flags
    )
{
  int rnk_n = 3;
  
  return PX(plan_cgc)(rnk_n, n, gc_below, gc_above,
      data, comm_cart, gc_flags);
}

PX(gcplan) PX(plan_rgc)(
    int rnk_n, const INT *n, const INT *gc_below, const INT *gc_above,
    R *data, MPI_Comm comm_cart, unsigned gc_flags
    )
{
  INT howmany = 1, *block = PFFT_DEFAULT_BLOCKS;
  
  return PX(plan_many_rgc)(
      rnk_n, n, howmany, block, gc_below, gc_above,
      data, comm_cart, gc_flags);
}

PX(gcplan) PX(plan_cgc)(
    int rnk_n, const INT *n, const INT *gc_below, const INT *gc_above,
    C *data, MPI_Comm comm_cart, unsigned gc_flags
    )
{
  INT howmany = 1, *block = PFFT_DEFAULT_BLOCKS;

  return PX(plan_many_cgc)(
      rnk_n, n, howmany, block, gc_below, gc_above,
      data, comm_cart, gc_flags);
}

void PX(exchange)(
    PX(gcplan) ths
    )
{
  if(ths==NULL)
    return;

  PFFT_START_TIMING(ths->comm_cart, ths->timer_exg->whole);
  PX(exchange_gc)(ths);
  PFFT_FINISH_TIMING(ths->timer_exg->whole);
  ths->timer_exg->iter++;
}


void PX(reduce)(
    PX(gcplan) ths
    )
{
  if(ths==NULL)
    return;

  ths->timer_red->iter++;
  PFFT_START_TIMING(ths->comm_cart, ths->timer_red->whole);
  PX(reduce_gc)(ths);
  PFFT_FINISH_TIMING(ths->timer_red->whole);
}



static void print_array(
     const R *data, const INT *n, const INT *start, const char *name, const int is_complex
     )
{
  INT k0, k1, k2, l=0;
  
  if( PX(prod_INT)(3, n) < 1)
    return;

  for(k0 = 0; k0 < n[0]; k0++){
    printf("%s(%td,%td:%td,%td:%td):\n", name, start[0]+k0, start[1],
        start[1]+n[1]-1, start[2], start[2]+n[2]-1);
    for(k1 = 0; k1 < n[1]; k1++){
      printf("  ");
      for(k2 = 0; k2 < n[2]; k2++, l++){
        if( is_complex )
          printf("  %.2e + %.2ei,", (double) data[2*l], (double) data[2*l+1]);
        else
          printf("  %.2e,", (double) data[l]);
      }
      printf("\n");
    }
  }
  printf("\n");
}


static void print_complex_array(
     const R *data, const INT *n, const INT *start, const char *name
     )
{
  print_array(data, n, start, name, 1);
}


static void print_real_array(
     const R *data, const INT *n, const INT *start, const char *name
     )
{
  print_array(data, n, start, name, 0);
}


void PX(apr_complex_3d)(
     const C *data,
     const INT *local_n, const INT *local_start,
     const char *name, MPI_Comm comm
     )
{
  int num_procs, proc_rank;
  
  MPI_Comm_size(comm, &num_procs);
  MPI_Comm_rank(comm, &proc_rank);

  fflush(stdout);
  for(int t=0; t<num_procs; t++){
    if(proc_rank == t){
      printf("Rank %d:\n", proc_rank);
      print_complex_array((R*) data, local_n, local_start, name);
      fflush(stdout);
    }
    MPI_Barrier(comm);
  }
}


void PX(apr_real_3d)(
     const R *data,
     const INT *local_n, const INT *local_start,
     const char *name, MPI_Comm comm
     )
{
  int num_procs, proc_rank;
  
  MPI_Comm_size(comm, &num_procs);
  MPI_Comm_rank(comm, &proc_rank);

  fflush(stdout);
  for(int t=0; t<num_procs; t++){
    if(proc_rank == t){
      printf("Rank %d:\n", proc_rank);
      print_real_array( data, local_n, local_start, name);
      fflush(stdout);
    }
    MPI_Barrier(comm);
  }
}


void PX(apr_complex_permuted_3d)(
     const C *data,
     const INT *local_n, const INT *local_start,
     int perm0, int perm1, int perm2,
     const char *name, MPI_Comm comm
     )
{
  int perm[3];
  INT local_start_perm[3], local_n_perm[3];
  
  perm[0] = perm0; perm[1] = perm1; perm[2] = perm2;
  
  for(int t=0; t<3; t++){
    local_n_perm[t] = local_n[perm[t]];
    local_start_perm[t] = local_start[perm[t]];
  }
  
  PX(apr_complex_3d)(data, local_n_perm, local_start_perm, name, comm);
}


void PX(apr_real_permuted_3d)(
     const R *data,
     const INT *local_n, const INT *local_start,
     int perm0, int perm1, int perm2,
     const char *name, MPI_Comm comm
     )
{
  int perm[3];
  INT local_start_perm[3], local_n_perm[3];
  
  perm[0] = perm0; perm[1] = perm1; perm[2] = perm2;
  
  for(int t=0; t<3; t++){
    local_n_perm[t] = local_n[perm[t]];
    local_start_perm[t] = local_start[perm[t]];
  }
  
  PX(apr_real_3d)(data, local_n_perm, local_start_perm, name, comm);
}











/****************
 * 3d interface *
 ***************/








INT PX(local_size_dft_3d)(
    const INT *n,
    MPI_Comm comm_cart, unsigned pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  int rnk_n = 3;

  return PX(local_size_dft)(rnk_n, n, comm_cart, pfft_flags,
      local_ni, local_i_start, local_no, local_o_start);
}

INT PX(local_size_dft_r2c_3d)(
    const INT *n,
    MPI_Comm comm_cart, unsigned pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  int rnk_n = 3;

  return PX(local_size_dft_r2c)(rnk_n, n, comm_cart, pfft_flags,
      local_ni, local_i_start, local_no, local_o_start);
}

INT PX(local_size_dft_c2r_3d)(
    const INT *n,
    MPI_Comm comm_cart, unsigned pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  int rnk_n = 3;

  return PX(local_size_dft_c2r)(rnk_n, n, comm_cart, pfft_flags,
      local_ni, local_i_start, local_no, local_o_start);
}

INT PX(local_size_r2r_3d)(
    const INT *n,
    MPI_Comm comm_cart, unsigned pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  int rnk_n = 3;

  return PX(local_size_r2r)(rnk_n, n, comm_cart, pfft_flags,
      local_ni, local_i_start, local_no, local_o_start);
}

PX(plan) PX(plan_dft_3d)(
    const INT *n, C *in, C *out, MPI_Comm comm_cart,
    int sign, unsigned pfft_flags
    )
{
  int rnk_n = 3;

  return PX(plan_dft)(rnk_n, n, in, out,
      comm_cart, sign, pfft_flags);
}

PX(plan) PX(plan_dft_r2c_3d)(
    const INT *n, R *in, C *out, MPI_Comm comm_cart,
    int sign, unsigned pfft_flags
    )
{
  int rnk_n = 3;
  
  return PX(plan_dft_r2c)(rnk_n, n, in, out,
      comm_cart, sign, pfft_flags);
}

PX(plan) PX(plan_dft_c2r_3d)(
    const INT *n, C *in, R *out, MPI_Comm comm_cart,
    int sign, unsigned pfft_flags
    )
{
  int rnk_n = 3;
  
  return PX(plan_dft_c2r)(rnk_n, n, in, out,
      comm_cart, sign, pfft_flags);
}

PX(plan) PX(plan_r2r_3d)(
    const INT *n, R *in, R *out, MPI_Comm comm_cart,
    const PX(r2r_kind) *kinds, unsigned pfft_flags
    )
{
  int rnk_n = 3;

  return PX(plan_r2r)(rnk_n, n, in, out,
      comm_cart, kinds, pfft_flags);
}

/* compute block size and offset for arbitrary process rank */
void PX(local_block_dft_3d)(
    const INT *n,
    MPI_Comm comm_cart, int pid, unsigned pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  int rnk_n = 3;

  PX(local_block_dft)(rnk_n, n, comm_cart, pid, pfft_flags,
      local_ni, local_i_start, local_no, local_o_start);
}

void PX(local_block_dft_r2c_3d)(
    const INT *n,
    MPI_Comm comm_cart, int pid, unsigned pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  int rnk_n = 3;

  PX(local_block_dft_r2c)(rnk_n, n, comm_cart, pid, pfft_flags,
      local_ni, local_i_start, local_no, local_o_start);
}

void PX(local_block_dft_c2r_3d)(
    const INT *n,
    MPI_Comm comm_cart, int pid, unsigned pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  int rnk_n = 3;

  PX(local_block_dft_c2r)(rnk_n, n, comm_cart, pid, pfft_flags,
      local_ni, local_i_start, local_no, local_o_start);
}

void PX(local_block_r2r_3d)(
    const INT *n,
    MPI_Comm comm_cart, int pid, unsigned pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  int rnk_n = 3;

  PX(local_block_r2r)(rnk_n, n, comm_cart, pid, pfft_flags,
      local_ni, local_i_start, local_no, local_o_start);
}


/* basic interface */
INT PX(local_size_dft)(
    int rnk_n, const INT *n,
    MPI_Comm comm_cart, unsigned pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  INT howmany = 1;

  return PX(local_size_many_dft)(rnk_n, n, n, n, howmany,
      PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, comm_cart, pfft_flags,
      local_ni, local_i_start, local_no, local_o_start);
}

INT PX(local_size_dft_r2c)(
    int rnk_n, const INT *n,
    MPI_Comm comm_cart, unsigned pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  INT howmany = 1;

  return PX(local_size_many_dft_r2c)(rnk_n, n, n, n, howmany,
      PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, comm_cart, pfft_flags,
      local_ni, local_i_start, local_no, local_o_start);
}

INT PX(local_size_dft_c2r)(
    int rnk_n, const INT *n,
    MPI_Comm comm_cart, unsigned pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  INT howmany = 1;

  return PX(local_size_many_dft_c2r)(rnk_n, n, n, n, howmany,
      PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, comm_cart, pfft_flags,
      local_ni, local_i_start, local_no, local_o_start);
}

INT PX(local_size_r2r)(
    int rnk_n, const INT *n,
    MPI_Comm comm_cart, unsigned pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  INT howmany = 1;

  return PX(local_size_many_r2r)(rnk_n, n, n, n, howmany,
      PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, comm_cart, pfft_flags,
      local_ni, local_i_start, local_no, local_o_start);
}



PX(plan) PX(plan_dft)(
    int rnk_n, const INT *n, C *in, C *out, MPI_Comm comm_cart,
    int sign, unsigned pfft_flags
    )
{
  INT howmany = 1;
  
  return PX(plan_many_dft)(
    rnk_n, n, n, n, howmany, PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS,
    in, out, comm_cart, sign, pfft_flags);
}

PX(plan) PX(plan_dft_r2c)(
    int rnk_n, const INT *n, R *in, C *out, MPI_Comm comm_cart,
    int sign, unsigned pfft_flags
    )
{
  INT howmany = 1;
  
  return PX(plan_many_dft_r2c)(
    rnk_n, n, n, n, howmany, PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS,
    in, out, comm_cart, sign, pfft_flags);
}

PX(plan) PX(plan_dft_c2r)(
    int rnk_n, const INT *n, C *in, R *out, MPI_Comm comm_cart,
    int sign, unsigned pfft_flags
    )
{
  INT howmany = 1;
  
  return PX(plan_many_dft_c2r)(
    rnk_n, n, n, n, howmany, PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS,
    in, out, comm_cart, sign, pfft_flags);
}

PX(plan) PX(plan_r2r)(
    int rnk_n, const INT *n, R *in, R *out, MPI_Comm comm_cart,
    const PX(r2r_kind) *kinds, unsigned pfft_flags
    )
{
  INT howmany = 1;
  
  return PX(plan_many_r2r)(
    rnk_n, n, n, n, howmany, PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS,
    in, out, comm_cart, kinds, pfft_flags);
}

void PX(local_block_dft)(
    int rnk_n, const INT *n,
    MPI_Comm comm_cart, int pid, unsigned pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  PX(local_block_many_dft)(rnk_n, n, n,
      PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, comm_cart, pid, pfft_flags,
      local_ni, local_i_start, local_no, local_o_start);
}

void PX(local_block_dft_r2c)(
    int rnk_n, const INT *n,
    MPI_Comm comm_cart, int pid, unsigned pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  PX(local_block_many_dft_r2c)(rnk_n, n, n,
      PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, comm_cart, pid, pfft_flags,
      local_ni, local_i_start, local_no, local_o_start);
}

void PX(local_block_dft_c2r)(
    int rnk_n, const INT *n,
    MPI_Comm comm_cart, int pid, unsigned pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  PX(local_block_many_dft_c2r)(rnk_n, n, n,
      PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, comm_cart, pid, pfft_flags,
      local_ni, local_i_start, local_no, local_o_start);
}

void PX(local_block_r2r)(
    int rnk_n, const INT *n,
    MPI_Comm comm_cart, int pid, unsigned pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  PX(local_block_many_r2r)(rnk_n, n, n,
      PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, comm_cart, pid, pfft_flags,
      local_ni, local_i_start, local_no, local_o_start);
}


/* functions to execute and destroy PX(plan) */

static void PX(execute_full)(
    const PX(plan) ths, R *in, R *out
    )
{
  int r;
  
  if(ths==NULL){
    PX(fprintf)(MPI_COMM_WORLD, stderr,
        "!!! Error: Can not execute PFFT Plan == NULL !!!\n");
    return;
  }

  /* set global variables for easy debbuging  */
#if PFFT_DEBUG_GVARS
  gdbg_comm_cart = &ths->comm_cart;
  gdbg_rnk_pm = ths->rnk_pm;
  gdbg_comms_pm = ths->comms_pm;
#endif

  r = ths->rnk_pm;

  if (ths->trafo_flag & PFFTI_TRAFO_C2R && ths->conjugate_in && ths->conjugate_out){
    R* conj_in  = (ths->conjugate_in  == ths->in)  ? in : out;
    R* conj_out = (ths->conjugate_out == ths->out) ? out : in;
    complex_conjugate(conj_in, conj_out, ths->rnk_n, ths->local_ni);
  }

  PFFT_START_TIMING(ths->comm_cart, ths->timer->whole);

  /* twiddle inputs in order to get outputs shifted by n/2 */
  PFFT_START_TIMING(ths->comm_cart, ths->timer->itwiddle);
  if(ths->pfft_flags & PFFT_SHIFTED_OUT)
    twiddle_input(ths, ths->in, ths->out, in, out);
  PFFT_FINISH_TIMING(ths->timer->itwiddle);

  PFFT_START_TIMING(ths->comm_cart, ths->timer->remap_3dto2d[0]);
  PX(execute_remap_3dto2d)(ths->remap_3dto2d[0], ths->in, ths->out, in, out);
  PFFT_FINISH_TIMING(ths->timer->remap_3dto2d[0]);

  execute_transposed(r, ths->serial_trafo, ths->global_remap,
      ths->timer->trafo, ths->timer->remap, ths->in, ths->out, in, out, ths->comm_cart);

  execute_transposed(r, &ths->serial_trafo[r+1], &ths->global_remap[r],
      &ths->timer->trafo[r+1], &ths->timer->remap[r], ths->in, ths->out, in, out, ths->comm_cart);
  
  PFFT_START_TIMING(ths->comm_cart, ths->timer->remap_3dto2d[1]);
  PX(execute_remap_3dto2d)(ths->remap_3dto2d[1], ths->in, ths->out, in, out);
  PFFT_FINISH_TIMING(ths->timer->remap_3dto2d[1]);

  /* twiddle outputs in order to get inputs shifted by n/2 */
  PFFT_START_TIMING(ths->comm_cart, ths->timer->otwiddle);
  if(ths->pfft_flags & PFFT_SHIFTED_IN)
    twiddle_output(ths, ths->in, ths->out, in, out);
  PFFT_FINISH_TIMING(ths->timer->otwiddle);

  if (ths->trafo_flag & PFFTI_TRAFO_R2C && ths->conjugate_in && ths->conjugate_out){
    R* conj_in  = (ths->conjugate_in  == ths->in)  ? in : out;
    R* conj_out = (ths->conjugate_out == ths->out) ? out : in;
    complex_conjugate(conj_in, conj_out, ths->rnk_n, ths->local_no);
  }

  PFFT_FINISH_TIMING(ths->timer->whole);
  ths->timer->iter++;
}

void PX(execute)(
    const PX(plan) ths
    )
{
    PX(execute_full)(ths, ths->in, ths->out);
}

void PX(execute_dft)(
    const PX(plan) ths,
    C * in, C * out
    )
{
    PX(execute_full)(ths, (R*) in, (R*) out);
}


void PX(execute_dft_r2c)(
    const PX(plan) ths,
    R * in, C * out
    )
{
    PX(execute_full)(ths, (R*)in, (R*)out);
}

void PX(execute_dft_c2r)(
    const PX(plan) ths,
    C * in, R * out
    )
{
    PX(execute_full)(ths, (R*)in, (R*)out);
}

void PX(execute_r2r)(
    const PX(plan) ths,
    R * in, R * out
    )
{
    PX(execute_full)(ths, in, out);
}

static int* malloc_and_transpose_int(
    int rnk_n, int rnk_pm, int transp,
    int *n
    )
{
  int *nt = PX(malloc_int)(rnk_n);

  for(int t=0; t<rnk_pm; t++)
    nt[t] = transp ? n[t+1] : n[t];
  
  nt[rnk_pm] = transp ? n[0] : n[rnk_pm];

  for(int t=rnk_pm+1; t<rnk_n; t++)
    nt[t] = n[t];

  return nt;
} 

static INT* malloc_and_transpose_INT(
    int rnk_n, int rnk_pm, int transp,
    INT *n
    )
{
  INT *nt = PX(malloc_INT)(rnk_n);

  for(int t=0; t<rnk_pm; t++)
    nt[t] = transp ? n[t+1] : n[t];
  
  nt[rnk_pm] = transp ? n[0] : n[rnk_pm];

  for(int t=rnk_pm+1; t<rnk_n; t++)
    nt[t] = n[t];

  return nt;
} 

/* twiddle inputs in order to get outputs shifted by n/2 */
static void twiddle_input(
    const PX(plan) ths,
    R *planned_in, R *planned_out,
    R *executed_in, R *executed_out
    )
{
  if(ths->itwiddle_in == NULL)
    return;

  INT howmany = ths->howmany, l;
  R factor;
  INT *n = malloc_and_transpose_INT(ths->rnk_n, ths->rnk_pm, ths->transp_flag & PFFT_TRANSPOSED_IN, ths->n);
  INT *ni = malloc_and_transpose_INT(ths->rnk_n, ths->rnk_pm, ths->transp_flag & PFFT_TRANSPOSED_IN, ths->ni);
  INT *local_ni = malloc_and_transpose_INT(ths->rnk_n, ths->rnk_pm, ths->transp_flag & PFFT_TRANSPOSED_IN, ths->local_ni);
  INT *local_ni_start = malloc_and_transpose_INT(ths->rnk_n, ths->rnk_pm, ths->transp_flag & PFFT_TRANSPOSED_IN, ths->local_ni_start);
  int *skip_trafos = malloc_and_transpose_int(ths->rnk_n, ths->rnk_pm, ths->transp_flag & PFFT_TRANSPOSED_IN, ths->skip_trafos);
  R *in  = (ths->itwiddle_in  == planned_in)  ? executed_in  : executed_out;
  R *out = (ths->itwiddle_out == planned_out) ? executed_out : executed_in;

  if( (ths->trafo_flag & PFFTI_TRAFO_C2C) || (ths->trafo_flag & PFFTI_TRAFO_C2R ) )
    howmany *= 2;
  
  INT local_n_total = PX(prod_INT)(ths->rnk_n, local_ni);
  int tempbool=ths->pfft_flags & PFFT_SHIFTED_IN;
/*#pragma omp parallel for schedule(static,16)*/
  
#pragma omp parallel for schedule(static,8) private(factor,l)
  for(INT k=0; k<local_n_total; k++){
    l = k;
    factor = 1.0;
    for(int t=ths->rnk_n-1; t>=0; t--){
      INT kt = l%local_ni[t];
      /* check for r2c/c2r padding elements and skipped trafos */
      INT temp=kt+local_ni_start[t];
      if(temp < ni[t]/2){
        if(!skip_trafos[t]){
          if(temp%2) factor*=-1.0;          

          /* take care of the extra twiddle (-1)^(n/2) if both (input AND output) are shifted */
          if(tempbool && (n[t]/2)%2)
            factor *= -1.0;
        }
      }
      l /= local_ni[t];
    }

    for(INT h=0; h<howmany; h++)
      out[howmany*k+h] = in[howmany*k+h] * factor;
//     fprintf(stderr, "pfft: api-basic: in[%2td] = %.2e + I* %.2e\n", k, ths->in[howmany*k], ths->in[howmany*k+1]);
  }

  free(n); free(ni); free(local_ni); free(local_ni_start); free(skip_trafos);
}


/* twiddle outputs in order to get inputs shifted by n/2 */
static void twiddle_output(
    const PX(plan) ths,
    R *planned_in, R *planned_out,
    R *executed_in, R *executed_out
    )
{
  if(ths->otwiddle_in == NULL)
    return;

  INT howmany = ths->howmany, l;
  R factor;
  INT *no = malloc_and_transpose_INT(ths->rnk_n, ths->rnk_pm, ths->transp_flag & PFFT_TRANSPOSED_OUT, ths->no);
  INT *local_no = malloc_and_transpose_INT(ths->rnk_n, ths->rnk_pm, ths->transp_flag & PFFT_TRANSPOSED_OUT, ths->local_no);
  INT *local_no_start = malloc_and_transpose_INT(ths->rnk_n, ths->rnk_pm, ths->transp_flag & PFFT_TRANSPOSED_OUT, ths->local_no_start);
  int *skip_trafos = malloc_and_transpose_int(ths->rnk_n, ths->rnk_pm, ths->transp_flag & PFFT_TRANSPOSED_OUT, ths->skip_trafos);
  R *in  = (ths->otwiddle_in  == planned_in)  ? executed_in  : executed_out;
  R *out = (ths->otwiddle_out == planned_out) ? executed_out : executed_in;
  
  if( (ths->trafo_flag & PFFTI_TRAFO_C2C) || (ths->trafo_flag & PFFTI_TRAFO_R2C ) )
    howmany *= 2;

  INT local_n_total = PX(prod_INT)(ths->rnk_n, local_no);
/*#pragma omp parallel for */
#pragma omp parallel for schedule(static,8) private(factor,l)
  for(INT k=0; k<local_n_total; k++){
    l = k;
    factor = 1.0;
    for(int t=ths->rnk_n-1; t>=0; t--){
      INT kt = l%local_no[t];
      INT temp=kt+local_no_start[t];
      /* check for r2c/c2r padding elements and skipped trafos */
      if(temp < no[t]/2 && !skip_trafos[t])
        if(temp%2) 
          factor *= -1.0;
      l /= local_no[t];
    }

    for(INT h=0; h<howmany; h++)
      out[howmany*k+h] = in[howmany*k+h] * factor;
//     fprintf(stderr, "pfft: api-basic: out[%2td] = %.2e + I* %.2e\n", k, ths->out[howmany*k], ths->out[howmany*k+1]);
  }

  free(no); free(local_no); free(local_no_start); free(skip_trafos);
}


static void execute_transposed(
    int rnk_pm, outrafo_plan *trafos, gtransp_plan *remaps,
    double *timer_trafo, double *timer_remap,
    R * plannedin, R * plannedout,
    R * in, R * out,
    MPI_Comm comm_cart
    )
{
  int t;
  
  for(t=0; t<rnk_pm; t++){
    PFFT_START_TIMING(comm_cart, timer_trafo[t]);
    PX(execute_outrafo)(trafos[t], plannedin, plannedout, in, out);
    PFFT_FINISH_TIMING(timer_trafo[t]);

    PFFT_START_TIMING(comm_cart, timer_remap[t]);
    PX(execute_gtransp)(remaps[t], plannedin, plannedout, in, out);
    PFFT_FINISH_TIMING(timer_remap[t]);
  }
  
  PFFT_START_TIMING(comm_cart, timer_trafo[t]);
  PX(execute_outrafo)(trafos[t], plannedin, plannedout, in, out);
  PFFT_FINISH_TIMING(timer_trafo[t]);
}

