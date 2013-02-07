#include <complex.h>
#include <pfft.h>

/* This program tests the over-/under-sampled PFFT.
 * We first compute an oversampled forward PFFT and second
 * revert it with an undersampled backward PFFT.
 * This should give back the initial input values
 * (up to the typical scaling factor) */

int main(int argc, char **argv){
  int np[2];
  ptrdiff_t n[3], ni[3], no[3];
  ptrdiff_t alloc_local_forw, alloc_local_back, alloc_local, howmany;
  ptrdiff_t local_ni[3], local_i_start[3];
  ptrdiff_t local_n[3], local_start[3];
  ptrdiff_t local_no[3], local_o_start[3];
  double err, *in;
  pfft_complex *out;
  pfft_plan plan_forw=NULL, plan_back=NULL;
  MPI_Comm comm_cart_2d;
  
  /* Set size of FFT and process mesh */
  ni[0] = ni[1] = ni[2] = 16;
  n[0] = 29; n[1] = 27; n[2] = 31;
  for(int t=0; t<3; t++)
    no[t] = ni[t];
  np[0] = 2; np[1] = 2;
  howmany = 1;
  
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
  alloc_local_forw = pfft_local_size_many_dft_r2c(3, n, ni, n, howmany,
      PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS,
      comm_cart_2d, PFFT_TRANSPOSED_OUT,
      local_ni, local_i_start, local_n, local_start);

  alloc_local_back = pfft_local_size_many_dft_c2r(3, n, n, no, howmany,
      PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS,
      comm_cart_2d, PFFT_TRANSPOSED_IN,
      local_n, local_start, local_no, local_o_start);

  /* Allocate enough memory for both trafos */
  alloc_local = (alloc_local_forw > alloc_local_back) ?
    alloc_local_forw : alloc_local_back;
  in  = pfft_alloc_real(2 * alloc_local);
  out = pfft_alloc_complex(alloc_local);

  /* Plan parallel forward FFT */
  plan_forw = pfft_plan_many_dft_r2c(
      3, n, ni, n, howmany, PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS,
      in, out, comm_cart_2d, PFFT_FORWARD, PFFT_TRANSPOSED_OUT| PFFT_MEASURE| PFFT_DESTROY_INPUT);

  /* Plan parallel backward FFT */
  plan_back = pfft_plan_many_dft_c2r(
      3, n, n, no, howmany, PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS,
      out, in, comm_cart_2d, PFFT_BACKWARD, PFFT_TRANSPOSED_IN| PFFT_MEASURE| PFFT_DESTROY_INPUT);

  /* Initialize input with random numbers */
  pfft_init_input_r2c(3, ni, local_ni, local_i_start,
      in);

  /* Execute parallel forward FFT */
  pfft_execute(plan_forw);
 
  /* execute parallel backward FFT */
  pfft_execute(plan_back);
  
  /* Scale data */
  for(ptrdiff_t l=0; l < local_ni[0] * local_ni[1] * local_ni[2]; l++)
    in[l] /= (n[0]*n[1]*n[2]);

  /* Print error of back transformed data */
  MPI_Barrier(MPI_COMM_WORLD);
  err = pfft_check_output_c2r(3, ni, local_ni, local_i_start, in, comm_cart_2d);
  pfft_printf(comm_cart_2d, "Error after one forward and backward trafo of size n=(%td, %td, %td):\n", n[0], n[1], n[2]); 
  pfft_printf(comm_cart_2d, "maxerror = %6.2e;\n", err);
  
  /* free mem and finalize MPI */
  pfft_destroy_plan(plan_forw);
  pfft_destroy_plan(plan_back);
  MPI_Comm_free(&comm_cart_2d);
  pfft_free(in); pfft_free(out);
  MPI_Finalize();
  return 0;
}
