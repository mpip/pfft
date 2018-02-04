#include <complex.h>
#include <pfft.h>

int main(int argc, char **argv)
{
  int np[3];
  ptrdiff_t n[3];
  ptrdiff_t alloc_local;
  ptrdiff_t local_ni[3], local_i_start[3];
  ptrdiff_t local_no[3], local_o_start[3];
  double err, *in;
  pfft_complex *out;
  pfft_plan plan_forw=NULL, plan_back=NULL;
  MPI_Comm comm_cart_3d;
  
  /* Set size of FFT and process mesh */
  n[0] = 29; n[1] = 27; n[2] = 31;
  np[0] = 2; np[1] = 2; np[2] = 2;
  
  /* Initialize MPI and PFFT */
  MPI_Init(&argc, &argv);
  pfft_init();

  /* Create three-dimensional process grid of size np[0] x np[1] x np[2], if possible */
  if( pfft_create_procmesh(2, MPI_COMM_WORLD, np, &comm_cart_3d) ){
    pfft_fprintf(MPI_COMM_WORLD, stderr, "Error: This test file only works with %d processes.\n", np[0]*np[1]*np[2]);
    MPI_Finalize();
    return 1;
  }
  
  /* Get parameters of data distribution */
  alloc_local = pfft_local_size_dft_r2c(2, n, comm_cart_3d, PFFT_TRANSPOSED_NONE,
      local_ni, local_i_start, local_no, local_o_start);

  pfft_fprintf(MPI_COMM_WORLD, stderr, "%td.\n", alloc_local);

  /* Allocate memory */
  in  = pfft_alloc_real(2 * alloc_local);
  out = pfft_alloc_complex(alloc_local);

  /* Plan parallel forward FFT */
  plan_forw = pfft_plan_dft_r2c(2,
      n, in, out, comm_cart_3d, PFFT_FORWARD, PFFT_TRANSPOSED_NONE| PFFT_MEASURE| PFFT_DESTROY_INPUT);
  
  /* Plan parallel backward FFT */
  plan_back = pfft_plan_dft_c2r(2, 
      n, out, in, comm_cart_3d, PFFT_BACKWARD, PFFT_TRANSPOSED_NONE| PFFT_MEASURE| PFFT_DESTROY_INPUT);

  /* Initialize input with random numbers */
  pfft_init_input_real(2, n, local_ni, local_i_start,
      in);

  /* execute parallel forward FFT */
  pfft_execute(plan_forw);

  /* clear the old input */
  pfft_clear_input_real(2, n, local_ni, local_i_start,
      in);
  
  /* execute parallel backward FFT */
  pfft_execute(plan_back);
  
  /* Scale data */
  for(ptrdiff_t l=0; l < local_ni[0] * local_ni[1] * local_ni[2]; l++)
    in[l] /= (n[0]*n[1]*n[2]);
  
  /* Print error of back transformed data */
  MPI_Barrier(MPI_COMM_WORLD);
  err = pfft_check_output_real(2, n, local_ni, local_i_start, in, comm_cart_3d);
  pfft_printf(comm_cart_3d, "Error after one forward and backward trafo of size n=(%td, %td, %td):\n", n[0], n[1], n[2]); 
  pfft_printf(comm_cart_3d, "maxerror = %6.2e;\n", err);
  
  /* free mem and finalize */
  pfft_destroy_plan(plan_forw);
  pfft_destroy_plan(plan_back);
  MPI_Comm_free(&comm_cart_3d);
  pfft_free(in); pfft_free(out);
  MPI_Finalize();
  return 0;
}
