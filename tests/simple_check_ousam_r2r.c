#include <complex.h>
#include <pfft.h>

/* This program tests the over-/under-sampled PFFT.
 * We first compute an oversampled forward PFFT and second
 * revert it with an undersampled backward PFFT.
 * This should give back the initial input values
 * (up to the typical scaling factor) */

int main(int argc, char **argv){
  int np[2];
  ptrdiff_t n[3], ni[3], no[3], N[3];
  ptrdiff_t alloc_local_forw, alloc_local_back, alloc_local, howmany;
  ptrdiff_t local_ni[3], local_i_start[3];
  ptrdiff_t local_n[3], local_start[3];
  ptrdiff_t local_no[3], local_o_start[3];
  double err, *in, *out;
  pfft_plan plan_forw=NULL, plan_back=NULL;
  MPI_Comm comm_cart_2d;
  fftw_r2r_kind kinds_forw[3], kinds_back[3];
  
  /* Set size of FFT and process mesh */
  ni[0] = ni[1] = ni[2] = 16;
  n[0] = 29; n[1] = 27; n[2] = 31;
  for(int t=0; t<3; t++)
    no[t] = ni[t];
  np[0] = 2; np[1] = 2;
  howmany = 1;
  
  /* Set PFFT kinds of 1d R2R trafos */
  kinds_forw[0] = PFFT_REDFT00; kinds_back[0] = PFFT_REDFT00;
  kinds_forw[1] = PFFT_REDFT01; kinds_back[1] = PFFT_REDFT10;
  kinds_forw[2] = PFFT_RODFT00; kinds_back[2] = PFFT_RODFT00;

  /* Set logical DFT sizes corresponding to FFTW manual:
   * for REDFT00 N=2*(n-1), for RODFT00 N=2*(n+1), otherwise N=2*n */
  N[0] = 2*(n[0]-1);
  N[1] = 2*n[1];
  N[2] = 2*(n[2]+1); 

  /* Initialize MPI and PFFT */
  MPI_Init(&argc, &argv);
  pfft_init();

  /* Create two-dimensional process grid of size np[0] x np[1], if possible */
  if( pfft_create_procmesh_2d(MPI_COMM_WORLD, np[0], np[1], &comm_cart_2d) ){
    pfft_fprintf(MPI_COMM_WORLD, stderr, "Error: This test file only works with %d processes.\n", np[0]*np[1]);
    MPI_Finalize();
    return 1;
  }

  /* Get parameters of data distribution */
  alloc_local_forw = pfft_local_size_many_r2r(3, n, ni, n, howmany,
      PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS,
      comm_cart_2d, PFFT_TRANSPOSED_NONE,
      local_ni, local_i_start, local_n, local_start);

  alloc_local_back = pfft_local_size_many_r2r(3, n, n, no, howmany,
      PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS,
      comm_cart_2d, PFFT_TRANSPOSED_NONE,
      local_n, local_start, local_no, local_o_start);

  /* Allocate enough memory for both trafos */
  alloc_local = (alloc_local_forw > alloc_local_back) ?
    alloc_local_forw : alloc_local_back;
  in  = fftw_alloc_real(alloc_local);
  out = fftw_alloc_real(alloc_local);

  /* Plan parallel forward FFT */
  plan_forw = pfft_plan_many_r2r(
      3, n, ni, n, howmany, PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS,
      in, out, comm_cart_2d, kinds_forw, PFFT_TRANSPOSED_NONE| PFFT_MEASURE| PFFT_DESTROY_INPUT);

  /* Plan parallel backward FFT */
  plan_back = pfft_plan_many_r2r(
      3, n, n, no, howmany, PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS,
      out, in, comm_cart_2d, kinds_back, PFFT_TRANSPOSED_NONE| PFFT_MEASURE| PFFT_DESTROY_INPUT);

  /* Initialize input with random numbers */
  pfft_init_input_real_3d(ni, local_ni, local_i_start,
      in);

  /* Execute parallel forward FFT */
  pfft_execute(plan_forw);

  /* clear the old input */
  pfft_clear_input_real_3d(ni, local_ni, local_i_start,
      in);
 
  /* execute parallel backward FFT */
  pfft_execute(plan_back);
  
  /* Scale data */
  for(ptrdiff_t l=0; l < local_ni[0] * local_ni[1] * local_ni[2]; l++)
    in[l] /= (N[0]*N[1]*N[2]);

  /* Print error of back transformed data */
  MPI_Barrier(MPI_COMM_WORLD);
  err = pfft_check_output_real_3d(ni, local_ni, local_i_start, in, comm_cart_2d);
  pfft_printf(comm_cart_2d, "Error after one forward and backward trafo of size n=(%td, %td, %td):\n", n[0], n[1], n[2]); 
  pfft_printf(comm_cart_2d, "maxerror = %6.2e;\n", err);
  
  /* free mem and finalize MPI */
  pfft_destroy_plan(plan_forw);
  pfft_destroy_plan(plan_back);
  MPI_Comm_free(&comm_cart_2d);
  fftw_free(in); fftw_free(out);
  MPI_Finalize();
  return 0;
}
