#include <complex.h>
#include <pfft.h>

int main(int argc, char **argv){
  int np[2];
  ptrdiff_t n[3];
  ptrdiff_t alloc_local;
  ptrdiff_t local_ni[3], local_i_start[3];
  ptrdiff_t local_no[3], local_o_start[3];
  pfft_complex *in, *out;
  pfft_plan plan=NULL;
  MPI_Comm comm_cart_2d;
  
  /* Set size of FFT and process mesh */
  n[0] = 2; n[1] = 2; n[2] = 4;
  np[0] = 2; np[1] = 2;
  
  /* Initialize MPI and PFFT */
  MPI_Init(&argc, &argv);
  pfft_init();

  /* Create two-dimensional process grid of size np[0] x np[1] */
  pfft_create_procmesh_2d(MPI_COMM_WORLD, np[0], np[1], &comm_cart_2d);

  /* Get parameters of data distribution */
  alloc_local = pfft_local_size_dft_3d(n, comm_cart_2d, PFFT_TRANSPOSED_OUT,
      local_ni, local_i_start, local_no, local_o_start);

  /* Allocate memory */
  in  = pfft_alloc_complex(alloc_local);
  out = pfft_alloc_complex(alloc_local);

  /* Plan parallel forward FFT */
  plan = pfft_plan_dft_3d(n, in, out, comm_cart_2d,
      PFFT_FORWARD, PFFT_TRANSPOSED_OUT| PFFT_MEASURE| PFFT_DESTROY_INPUT);

  /* Initialize input with random numbers */
  pfft_init_input_c2c_3d(n, local_ni, local_i_start,
      in);

  /* Execute parallel forward FFT */
  pfft_execute(plan);
  
  /* free mem and finalize MPI */
  pfft_destroy_plan(plan);
  MPI_Comm_free(&comm_cart_2d);
  pfft_free(in); pfft_free(out);
  MPI_Finalize();
  return 0;
}
