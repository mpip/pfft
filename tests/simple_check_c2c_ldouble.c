#include <complex.h>
#include <pfft.h>

int main(int argc, char **argv)
{
  int np[2];
  ptrdiff_t n[3];
  ptrdiff_t alloc_local;
  ptrdiff_t local_ni[3], local_i_start[3];
  ptrdiff_t local_no[3], local_o_start[3];
  long double err;
  pfftl_complex *in, *out;
  pfftl_plan plan_forw=NULL, plan_back=NULL;
  MPI_Comm comm_cart_2d;
  
  /* Set size of FFT and process mesh */
  n[0] = 29; n[1] = 27; n[2] = 31;
  np[0] = 2; np[1] = 2;
  
  /* Initialize MPI and PFFT */
  MPI_Init(&argc, &argv);
  pfftl_init();

  /* Create two-dimensional process grid of size np[0] x np[1], if possible */
  if( pfftl_create_procmesh_2d(MPI_COMM_WORLD, np[0], np[1], &comm_cart_2d) ){
    pfftl_fprintf(MPI_COMM_WORLD, stderr, "Error: This test file only works with %d processes.\n", np[0]*np[1]);
    MPI_Finalize();
    return 1;
  }
  
  /* Get parameters of data distribution */
  alloc_local = pfftl_local_size_dft_3d(n, comm_cart_2d, PFFT_TRANSPOSED_NONE,
      local_ni, local_i_start, local_no, local_o_start);

  /* Allocate memory */
  in  = pfftl_alloc_complex(alloc_local);
  out = pfftl_alloc_complex(alloc_local);

  /* Plan parallel forward FFT */
  plan_forw = pfftl_plan_dft_3d(
      n, in, out, comm_cart_2d, PFFT_FORWARD, PFFT_TRANSPOSED_NONE| PFFT_MEASURE| PFFT_DESTROY_INPUT);
  
  /* Plan parallel backward FFT */
  plan_back = pfftl_plan_dft_3d(
      n, out, in, comm_cart_2d, PFFT_BACKWARD, PFFT_TRANSPOSED_NONE| PFFT_MEASURE| PFFT_DESTROY_INPUT);

  /* Initialize input with random numbers */
  pfftl_init_input_complex_3d(n, local_ni, local_i_start,
      in);
  
  /* execute parallel forward FFT */
  pfftl_execute(plan_forw);
  
  /* execute parallel backward FFT */
  pfftl_execute(plan_back);
  
  /* Scale data */
  ptrdiff_t l;
  for(l=0; l < local_ni[0] * local_ni[1] * local_ni[2]; l++)
    in[l] /= (n[0]*n[1]*n[2]);

  /* Print error of back transformed data */
  err = pfftl_check_output_complex_3d(n, local_ni, local_i_start, in, comm_cart_2d);
  pfftl_printf(comm_cart_2d, "Error after one forward and backward trafo of size n=(%td, %td, %td):\n", n[0], n[1], n[2]); 
  pfftl_printf(comm_cart_2d, "maxerror = %6.2Le;\n", err);
  
  /* free mem and finalize */
  pfftl_destroy_plan(plan_forw);
  pfftl_destroy_plan(plan_back);
  MPI_Comm_free(&comm_cart_2d);
  pfftl_free(in); pfftl_free(out);
  MPI_Finalize();
  return 0;
}
