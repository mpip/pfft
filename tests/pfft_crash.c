#include <complex.h>
#include <pfft.h>

int main(int argc, char **argv)
{
  int np[1];
  ptrdiff_t n[2];
  ptrdiff_t alloc_local;
  ptrdiff_t local_ni[2], local_i_start[2];
  ptrdiff_t local_no[2], local_o_start[2];
  MPI_Comm comm_cart_1d;

  /* Set size of FFT and process mesh */
  n[0] = 3; n[1] = 3;
  np[0] = 2;

  /* Initialize MPI and PFFT */
  MPI_Init(&argc, &argv);
  pfft_init();

  MPI_Comm_size(MPI_COMM_WORLD, np);

  /* Create two-dimensional process grid of size np[0] x np[1], if possible */
  if( pfft_create_procmesh(1, MPI_COMM_WORLD, np, &comm_cart_1d) ){
    pfft_fprintf(MPI_COMM_WORLD, stderr, "Error: This test file only works with %d processes.\n", np[0]);
    MPI_Finalize();
    return 1;
  }

  /* Get parameters of data distribution */
  alloc_local = pfft_local_size_dft(2, n, comm_cart_1d, PFFT_TRANSPOSED_NONE,
      local_ni, local_i_start, local_no, local_o_start);

//   int rnk;
//   MPI_Comm_rank(comm_cart_1d, &rnk);
//   for(int r=0; r<np[0]; r++){
//     if(r==rnk){
//       printf("local_ni = [%td, %td]\n", local_ni[0], local_ni[1]);
//       printf("local_i_start = [%td, %td]\n", local_i_start[0], local_i_start[1]);
//       printf("local_no = [%td, %td]\n", local_no[0], local_no[1]);
//       printf("local_o_start = [%td, %td]\n", local_o_start[0], local_o_start[1]);
//     }
//     MPI_Barrier(comm_cart_1d);
//   }

  /* Allocate memory */
  pfft_complex *in  = pfft_alloc_complex(alloc_local);
  pfft_complex *out = pfft_alloc_complex(alloc_local);

  /* Plan parallel forward FFT */
//   const unsigned flags = PFFT_TRANSPOSED_NONE| PFFT_ESTIMATE;
//   pfft_plan plan1 = pfft_plan_dft(
//       2, n, in, in, comm_cart_1d, PFFT_FORWARD, flags);
//   pfft_plan plan2 = pfft_plan_dft(
//       2, n, in, in, comm_cart_1d, PFFT_FORWARD, flags);

  pfft_plan plan1 = pfft_plan_dft(
      2, n, in, out, comm_cart_1d, PFFT_FORWARD, PFFT_TRANSPOSED_NONE| PFFT_ESTIMATE| PFFT_PRESERVE_INPUT);
  pfft_plan plan2 = pfft_plan_dft(
      2, n, in, out, comm_cart_1d, PFFT_FORWARD, PFFT_TRANSPOSED_NONE| PFFT_ESTIMATE| PFFT_PRESERVE_INPUT);

  /* free mem and finalize */
  pfft_destroy_plan(plan1);
  pfft_destroy_plan(plan2);
  MPI_Comm_free(&comm_cart_1d);
  pfft_free(in); pfft_free(out);
  MPI_Finalize();
  return 0;
}
