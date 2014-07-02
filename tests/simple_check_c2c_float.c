#include <complex.h>
#include <pfft.h>

int main(int argc, char **argv)
{
  int np[2];
  ptrdiff_t n[3];
  ptrdiff_t alloc_local;
  ptrdiff_t local_ni[3], local_i_start[3];
  ptrdiff_t local_no[3], local_o_start[3];
  float err;
  pfftf_complex *in, *out;
  pfftf_plan plan_forw=NULL, plan_back=NULL;
  MPI_Comm comm_cart_2d;
  
  /* Set size of FFT and process mesh */
  n[0] = 29; n[1] = 27; n[2] = 31;
  np[0] = 2; np[1] = 2;
  
  /* Initialize MPI and PFFT */
  MPI_Init(&argc, &argv);
  pfftf_init();

  /* Create two-dimensional process grid of size np[0] x np[1], if possible */
  if( pfftf_create_procmesh_2d(MPI_COMM_WORLD, np[0], np[1], &comm_cart_2d) ){
    pfftf_fprintf(MPI_COMM_WORLD, stderr, "Error: This test file only works with %d processes.\n", np[0]*np[1]);
    MPI_Finalize();
    return 1;
  }
  
  /* Get parameters of data distribution */
  alloc_local = pfftf_local_size_dft_3d(n, comm_cart_2d, PFFT_TRANSPOSED_NONE,
      local_ni, local_i_start, local_no, local_o_start);

  /* Allocate memory */
  in  = pfftf_alloc_complex(alloc_local);
  out = pfftf_alloc_complex(alloc_local);

  /* Plan parallel forward FFT */
  plan_forw = pfftf_plan_dft_3d(
      n, in, out, comm_cart_2d, PFFT_FORWARD, PFFT_TRANSPOSED_NONE| PFFT_MEASURE| PFFT_DESTROY_INPUT);
  
  /* Plan parallel backward FFT */
  plan_back = pfftf_plan_dft_3d(
      n, out, in, comm_cart_2d, PFFT_BACKWARD, PFFT_TRANSPOSED_NONE| PFFT_MEASURE| PFFT_DESTROY_INPUT);

  /* Initialize input with random numbers */
  pfftf_init_input_complex_3d(n, local_ni, local_i_start,
      in);
  
  /* execute parallel forward FFT */
  pfftf_execute(plan_forw);
  
  /* execute parallel backward FFT */
  pfftf_execute(plan_back);
  
  /* Scale data */
  ptrdiff_t l;
  for(l=0; l < local_ni[0] * local_ni[1] * local_ni[2]; l++)
    in[l] /= (n[0]*n[1]*n[2]);

  /* Print error of back transformed data */
  err = pfftf_check_output_complex_3d(n, local_ni, local_i_start, in, comm_cart_2d);
  pfftf_printf(comm_cart_2d, "Error after one forward and backward trafo of size n=(%td, %td, %td):\n", n[0], n[1], n[2]); 
  pfftf_printf(comm_cart_2d, "maxerror = %6.2e;\n", err);
  
  /* free mem and finalize */
  pfftf_destroy_plan(plan_forw);
  pfftf_destroy_plan(plan_back);
  MPI_Comm_free(&comm_cart_2d);
  pfftf_free(in); pfftf_free(out);
  MPI_Finalize();
  return 0;
}
