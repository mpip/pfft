#include <complex.h>
#include <pfft.h>

int main(int argc, char **argv)
{
  int np[2];
  ptrdiff_t n[3];
  ptrdiff_t alloc_local;
  ptrdiff_t local_ni[3], local_i_start[3];
  ptrdiff_t local_no[3], local_o_start[3];
  double err;
  double *in, *in2;
  pfft_complex *out, *out2;
  pfft_plan plan_forw=NULL, plan_back=NULL;
  MPI_Comm comm_cart_2d;
  
  /* Set size of FFT and process mesh */
  n[0] = 29; n[1] = 27; n[2] = 31;
  np[0] = 1; np[1] = 1;
  
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
  alloc_local = pfft_local_size_dft_r2c_3d(n, comm_cart_2d, PFFT_TRANSPOSED_OUT,
      local_ni, local_i_start, local_no, local_o_start);

  printf("%td %td %td \n", local_no[0], local_no[1], local_no[2]);


  /* Allocate memory */
  in  = pfft_alloc_real (2 * alloc_local);
  out = pfft_alloc_complex(alloc_local);

  in2  = pfft_alloc_real(2 * alloc_local);
  out2 = pfft_alloc_complex(alloc_local);

  /* Plan parallel forward FFT */
  plan_forw = pfft_plan_dft_r2c_3d(
      n, in2, out2, comm_cart_2d, PFFT_FORWARD, PFFT_TRANSPOSED_OUT| PFFT_MEASURE| PFFT_DESTROY_INPUT);
  
  /* Plan parallel backward FFT */
  plan_back = pfft_plan_dft_c2r_3d(
      n, out2, in2, comm_cart_2d, PFFT_BACKWARD, PFFT_TRANSPOSED_IN| PFFT_MEASURE| PFFT_DESTROY_INPUT);

  /* Initialize input with random numbers */
  pfft_init_input_r2c(3, n, local_ni, local_i_start,
      in);
  
  memset(in2, 0, sizeof(pfft_complex) * alloc_local);
  memset(out2, 0, sizeof(pfft_complex) * alloc_local);
  /* execute parallel forward FFT */
  pfft_execute_dft_r2c(plan_forw, in, out);
  memset(in2, 0, sizeof(pfft_complex) * alloc_local);
  memset(out2, 0, sizeof(pfft_complex) * alloc_local);
  /* execute parallel backward FFT */
  pfft_execute_dft_c2r(plan_back, out, in);
  memset(in2, 0, sizeof(pfft_complex) * alloc_local);
  memset(out2, 0, sizeof(pfft_complex) * alloc_local);
  
  /* Scale data */
  ptrdiff_t l;
  for(l=0; l < local_ni[0] * local_ni[1] * local_ni[2]; l++)
    in[l] /= (n[0]*n[1]*n[2]);

  /* Print error of back transformed data */
  err = pfft_check_output_c2r(3, n, local_ni, local_i_start, in, comm_cart_2d);
  pfft_printf(comm_cart_2d, "Error after one forward and backward trafo of size n=(%td, %td, %td):\n", n[0], n[1], n[2]); 
  pfft_printf(comm_cart_2d, "maxerror = %6.2e;\n", err);
  
  /* free mem and finalize */
  pfft_destroy_plan(plan_forw);
  pfft_destroy_plan(plan_back);
  MPI_Comm_free(&comm_cart_2d);
  pfft_free(in); pfft_free(out);
  MPI_Finalize();
  return 0;
}
