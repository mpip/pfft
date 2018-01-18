#include <complex.h>
#include <pfft.h>

static void printarray(double * a, int n0, int n1)
{
  for(ptrdiff_t l=0; l < n0 * n1; l++) {
//    pfft_fprintf(MPI_COMM_WORLD, stderr, "%04.0f, ", a[l]);
    pfft_fprintf(MPI_COMM_WORLD, stderr, "%05.1f, ", a[l]);
    if ((l + 1) % n1 == 0) {
       pfft_fprintf(MPI_COMM_WORLD, stderr, "\n", a[l]);
   }
  }
}
static void gatherarray(double * a, double * in, int local_start0,
        int local_n0, int n0, int local_start1, int local_n1, int n1)
{
  int pid;
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);

  memset(a, 0, sizeof(double) * n0 * n1);
  for(ptrdiff_t l=0; l < local_n0 * local_n1; l++) {
    int x = l / local_n1, y = l % local_n1;
    a[(local_start0 + x) * n1 + local_start1 + y] = in[l];
// + pid * 1000;
  }
  MPI_Allreduce(MPI_IN_PLACE, a, n0 * n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

int main(int argc, char **argv)
{
  int np[2];
  ptrdiff_t n[2];
  ptrdiff_t alloc_local;
  ptrdiff_t local_ni[2], local_i_start[2];
  ptrdiff_t local_no[2], local_o_start[2];
  double err, *in;
  pfft_complex *out;
  pfft_plan plan_forw=NULL, plan_back=NULL;
  MPI_Comm comm_cart_2d;
  
  /* Set size of FFT and process mesh */
  n[0] = 4; n[1] = 4;
  np[0] = 2; np[1] = 2;
  
  double *truein = calloc(sizeof(double), n[0] * n[1]);
  double *trueout = calloc(sizeof(double), n[0] * (n[1] / 2 + 1) * 2);
  int pid;

  /* Initialize MPI and PFFT */
  MPI_Init(&argc, &argv);
  pfft_init();

  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  /* Create three-dimensional process grid of size np[0] x np[1], if possible */
  if( pfft_create_procmesh(2, MPI_COMM_WORLD, np, &comm_cart_2d) ){
    pfft_fprintf(MPI_COMM_WORLD, stderr, "Error: This test file only works with %d processes.\n", np[0]*np[1]);
    MPI_Finalize();
    return 1;
  }
  
  /* Get parameters of data distribution */
  alloc_local = pfft_local_size_dft_r2c(2, n, comm_cart_2d, PFFT_TRANSPOSED_NONE,
      local_ni, local_i_start, local_no, local_o_start);

  /* Allocate memory */
  in  = pfft_alloc_real(2 * alloc_local);
  out = pfft_alloc_complex(alloc_local);

  pfft_fprintf(MPI_COMM_WORLD, stderr, "planning forward\n");
  /* Plan parallel forward FFT */
  plan_forw = pfft_plan_dft_r2c(2,
      n, in, out, comm_cart_2d, PFFT_FORWARD, PFFT_TRANSPOSED_NONE| PFFT_ESTIMATE| PFFT_DESTROY_INPUT);
  
  pfft_fprintf(MPI_COMM_WORLD, stderr, "planning backward\n");
  /* Plan parallel backward FFT */
  plan_back = pfft_plan_dft_c2r(2, 
      n, out, in, comm_cart_2d, PFFT_BACKWARD, PFFT_TRANSPOSED_NONE| PFFT_ESTIMATE| PFFT_DESTROY_INPUT);

  /* Initialize input with random numbers */
  pfft_init_input_real(2, n, local_ni, local_i_start,
      in);

  for(ptrdiff_t l=0; l < local_ni[0] * local_ni[1]; l++) {
    in[l] = pid * 100 + l;
  //  if(l == 0)
  //      in[l] = 1;
  }
  gatherarray(truein, in, local_i_start[0], local_ni[0], n[0],
                          local_i_start[1], local_ni[1], n[1]);
  printarray(truein, n[0], n[1]);

  pfft_fprintf(MPI_COMM_WORLD, stderr, "running forward\n");
  /* execute parallel forward FFT */
  pfft_execute(plan_forw);


  pfft_fprintf(MPI_COMM_WORLD, stderr, "as real \n");
  gatherarray(truein, out, local_i_start[0], local_ni[0], n[0],
                          local_i_start[1], local_ni[1], n[1]);
  printarray(truein, n[0], n[1]);


  pfft_fprintf(MPI_COMM_WORLD, stderr, "as complex \n");
  gatherarray(trueout, out, local_o_start[0], local_no[0], n[0],
                            2 * local_o_start[1],
                            2 * local_no[1], 2 *(n[1] / 2 + 1));

  printarray(trueout, n[0], (n[1] / 2 + 1) * 2);

  pfft_fprintf(MPI_COMM_WORLD, stderr, "running backward\n");
  /* clear the old input */
  pfft_clear_input_real(2, n, local_ni, local_i_start,
      in);
  /* execute parallel backward FFT */
  pfft_execute(plan_back);

  /* Scale data */
  for(ptrdiff_t l=0; l < local_ni[0] * local_ni[1]; l++)
    in[l] /= (n[0]*n[1]);

  gatherarray(truein, in, local_i_start[0], local_ni[0], n[0],
                          local_i_start[1], local_ni[1], n[1]);
  printarray(truein, n[0], n[1]);

  /* Print error of back transformed data */
  MPI_Barrier(MPI_COMM_WORLD);
  err = pfft_check_output_real(2, n, local_ni, local_i_start, in, comm_cart_2d);
  pfft_printf(comm_cart_2d, "Error after one forward and backward trafo of size n=(%td, %td):\n", n[0], n[1]); 
  pfft_printf(comm_cart_2d, "maxerror = %6.2e;\n", err);

  /* free mem and finalize */
  pfft_destroy_plan(plan_forw);
  pfft_destroy_plan(plan_back);
  MPI_Comm_free(&comm_cart_2d);
  pfft_free(in); pfft_free(out);
  MPI_Finalize();
  return 0;
}
