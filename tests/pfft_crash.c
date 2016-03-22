#include <complex.h>
#include <pfft.h>

int main(int argc, char **argv)
{
  int np[1];
  ptrdiff_t n[2];
  ptrdiff_t alloc_local;
  ptrdiff_t local_ni[2], local_i_start[2];
  ptrdiff_t local_no[2], local_o_start[2];
  double err;
  pfft_complex *in, *out;
  pfft_plan plan_forw=NULL, plan_back=NULL;
  MPI_Comm comm_cart_1d;

  /* Set size of FFT and process mesh */
  n[0] = 200; n[1] = 200;
  np[0] = 12;

  /* Initialize MPI and PFFT */
  MPI_Init(&argc, &argv);
  pfft_init();

  /* Create two-dimensional process grid of size np[0] x np[1], if possible */
  if( pfft_create_procmesh(1, MPI_COMM_WORLD, np, &comm_cart_1d) ){
    pfft_fprintf(MPI_COMM_WORLD, stderr, "Error: This test file only works with %d processes.\n", np[0]);
    MPI_Finalize();
    return 1;
  }

  /* Get parameters of data distribution */
  alloc_local = pfft_local_size_dft(2, n, comm_cart_1d, PFFT_TRANSPOSED_NONE,
      local_ni, local_i_start, local_no, local_o_start);

  /* Allocate memory */
  in  = pfft_alloc_complex(alloc_local);
  out = pfft_alloc_complex(alloc_local);

  /* Plan parallel forward FFT */
  plan_forw = pfft_plan_dft(
      2, n, in, out, comm_cart_1d, PFFT_FORWARD, PFFT_TRANSPOSED_NONE| PFFT_ESTIMATE| PFFT_PRESERVE_INPUT);

  /* Plan parallel backward FFT */
  plan_back = pfft_plan_dft(
      2, n, out, in, comm_cart_1d, PFFT_BACKWARD, PFFT_TRANSPOSED_NONE| PFFT_ESTIMATE| PFFT_PRESERVE_INPUT);

  /* Initialize input with random numbers */
  pfft_init_input_complex(2, n, local_ni, local_i_start,
      in);

  /* execute parallel forward FFT */
  pfft_execute(plan_forw);

  /* clear the old input */
  pfft_clear_input_complex(2, n, local_ni, local_i_start,
      in);

  /* execute parallel backward FFT */
  pfft_execute(plan_back);

  /* Scale data */
  for(ptrdiff_t l=0; l < local_ni[0] * local_ni[1]; l++)
    in[l] /= (n[0]*n[1]);

  /* Print error of back transformed data */
  err = pfft_check_output_complex(2, n, local_ni, local_i_start, in, comm_cart_1d);
  pfft_printf(comm_cart_1d, "Error after one forward and backward trafo of size n=(%td, %td):\n", n[0], n[1]); 
  pfft_printf(comm_cart_1d, "maxerror = %6.2e;\n", err);

  /* free mem and finalize */
  pfft_destroy_plan(plan_forw);
  pfft_destroy_plan(plan_back);
  MPI_Comm_free(&comm_cart_1d);
  pfft_free(in); pfft_free(out);
  MPI_Finalize();
  return 0;
}
