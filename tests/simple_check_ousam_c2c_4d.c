#include <complex.h>
#include <pfft.h>

/* This program tests the over-/under-sampled PFFT.
 * We first compute an oversampled forward PFFT and second
 * revert it with an undersampled backward PFFT.
 * This should give back the initial input values
 * (up to the typical scaling factor) */

int main(int argc, char **argv)
{
  int np[2];
  ptrdiff_t n[4], ni[4], no[4];
  ptrdiff_t alloc_local_forw, alloc_local_back, alloc_local, howmany;
  ptrdiff_t local_ni[4], local_i_start[4];
  ptrdiff_t local_n[4], local_start[4];
  ptrdiff_t local_no[4], local_o_start[4];
  double err;
  pfft_complex *in, *out;
  pfft_plan plan_forw=NULL, plan_back=NULL;
  MPI_Comm comm_cart_2d;
  
  /* Set size of FFT and process mesh */
  ni[0] = ni[1] = ni[2] = ni[3] = 8;
  n[0] = 13; n[1] = 14; n[2] = 19; n[3] = 17;
  for(int t=0; t<4; t++)
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
  alloc_local_forw = pfft_local_size_many_dft(4, n, ni, n, howmany,
      PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS,
      comm_cart_2d, PFFT_TRANSPOSED_NONE,
      local_ni, local_i_start, local_n, local_start);

  alloc_local_back = pfft_local_size_many_dft(4, n, n, no, howmany,
      PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS,
      comm_cart_2d, PFFT_TRANSPOSED_NONE,
      local_n, local_start, local_no, local_o_start);

  /* Allocate enough memory for both trafos */
  alloc_local = (alloc_local_forw > alloc_local_back) ?
    alloc_local_forw : alloc_local_back;
  in  = pfft_alloc_complex(alloc_local);
  out = pfft_alloc_complex(alloc_local);

  /* Plan parallel forward FFT */
  plan_forw = pfft_plan_many_dft(
      4, n, ni, n, howmany, PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS,
      in, out, comm_cart_2d, PFFT_FORWARD, PFFT_TRANSPOSED_NONE| PFFT_MEASURE| PFFT_DESTROY_INPUT);

  /* Plan parallel backward FFT */
  plan_back = pfft_plan_many_dft(
      4, n, n, no, howmany, PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS,
      out, in, comm_cart_2d, PFFT_BACKWARD, PFFT_TRANSPOSED_NONE| PFFT_MEASURE| PFFT_DESTROY_INPUT);

  /* Initialize input with random numbers */
  pfft_init_input_c2c(4, ni, local_ni, local_i_start,
      in);

  /* execute parallel forward FFT */
  pfft_execute(plan_forw);
  
  /* execute parallel backward FFT */
  pfft_execute(plan_back);
  
  /* Scale data */
  for(ptrdiff_t l=0; l < local_ni[0] * local_ni[1] * local_ni[2] * local_ni[3]; l++)
    in[l] /= (n[0]*n[1]*n[2]*n[3]);
  
  /* Print error of back transformed data */
  err = pfft_check_output_c2c(4, ni, local_ni, local_i_start, in, comm_cart_2d);
  pfft_printf(comm_cart_2d, "Error after one forward and backward trafo of size n=(%td, %td, %td, %td):\n", n[0], n[1], n[2], n[3]); 
  pfft_printf(comm_cart_2d, "maxerror = %6.2e;\n", err);

  /* free mem and finalize */
  pfft_destroy_plan(plan_forw);
  pfft_destroy_plan(plan_back);
  MPI_Comm_free(&comm_cart_2d);
  pfft_free(in); pfft_free(out);
  MPI_Finalize();
  return 0;
}
