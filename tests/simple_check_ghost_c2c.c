#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <complex.h>
/* return value 0 on success */

#include <pfft.h>

#define MAX(a,b) (((a)>(b))?(a):(b))

static void init_parameters(
    int argc, char **argv,
    ptrdiff_t *n, int *np, ptrdiff_t *gc_below, ptrdiff_t *gc_above,
    int* verbose
    )
{
  pfft_get_args(argc, argv, "-pfft_n", 3, PFFT_PTRDIFF_T, n);
  pfft_get_args(argc, argv, "-pfft_np", 3, PFFT_INT, np);
  pfft_get_args(argc, argv, "-pfft_gcb", 3, PFFT_PTRDIFF_T, gc_below);
  pfft_get_args(argc, argv, "-pfft_gca", 3, PFFT_PTRDIFF_T, gc_above);
  pfft_get_args(argc, argv, "-pfft_verb", 1, PFFT_INT, verbose);
}


int main(int argc, char **argv){
  ptrdiff_t n[3], gc_below[3], gc_above[3];
  ptrdiff_t local_ni[3], local_i_start[3];
  ptrdiff_t local_no[3], local_o_start[3];
  ptrdiff_t local_ngc[3], local_gc_start[3];
  ptrdiff_t alloc_local, alloc_local_gc;
  int np[3], rnk_self, size, verbose;
  double err;
  MPI_Comm comm_cart_2d;
  pfft_complex *data;
  pfft_gcplan ths;
  
  MPI_Init(&argc, &argv);
  pfft_init();
  MPI_Comm_rank(MPI_COMM_WORLD, &rnk_self);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  /* default values */
  n[0] = n[1] = n[2] = 8; /*  n[0] = 3; n[1] = 5; n[2] = 7;*/
  np[0]=2; np[1]=2; np[2] = 1;

  verbose = 0;
  for(int t=0; t<3; t++){
    gc_below[t] = 0;
    gc_above[t] = 0;
  }
  gc_below[0] = 0;
  gc_above[0] = 8;

  /* set values by commandline */
  init_parameters(argc, argv, n, np, gc_below, gc_above, &verbose);

  /* Create two-dimensional process grid of size np[0] x np[1], if possible */
  if( pfft_create_procmesh_2d(MPI_COMM_WORLD, np[0], np[1], &comm_cart_2d) ){
    pfft_fprintf(MPI_COMM_WORLD, stderr, "Error: This test file only works with %d processes.\n", np[0]*np[1]);
    MPI_Finalize();
    return 1;
  }

  /* Get parameters of data distribution */
  alloc_local = pfft_local_size_dft_3d(n, comm_cart_2d, PFFT_TRANSPOSED_NONE,
      local_ni, local_i_start, local_no, local_o_start);

  alloc_local_gc = pfft_local_size_gc_3d(
      local_ni, local_i_start, alloc_local, gc_below, gc_above,
      local_ngc, local_gc_start);

  /* Allocate memory */
  data = pfft_alloc_complex(alloc_local_gc);

  /* Plan parallel ghost cell send */
  ths = pfft_plan_cgc_3d(n, gc_below, gc_above,
      data, comm_cart_2d, PFFT_GC_NONTRANSPOSED);

  /* Initialize input with random numbers */
  pfft_init_input_c2c_3d(n, local_ni, local_i_start,
      data);

  /* check gcell input */
  if(verbose)
    pfft_apr_complex_3d(data, local_ni, local_i_start, "gcell input", comm_cart_2d);

  /* Execute parallel ghost cell send */
  pfft_exchange(ths);

  /* Check gcell output */
  if(verbose)
    pfft_apr_complex_3d(data, local_ngc, local_gc_start, "exchanged gcells", comm_cart_2d);
  
  /* Execute adjoint parallel ghost cell send */
  pfft_reduce(ths);

  /* check input */
  if(verbose)
    pfft_apr_complex_3d(data, local_no, local_o_start, "reduced gcells", comm_cart_2d);

  /* Scale data */
  for(ptrdiff_t l=0; l < local_ni[0] * local_ni[1] * local_ni[2]; l++)
    data[l] /= 3;

  /* Print error of back transformed data */
  MPI_Barrier(comm_cart_2d);
  err = pfft_check_output_c2c_3d(n, local_ni, local_i_start, data, comm_cart_2d);
  pfft_printf(comm_cart_2d, "Error after one gcell exchange and reduce of size n=(%td, %td, %td),\n", n[0], n[1], n[2]); 
  pfft_printf(comm_cart_2d, "gc_below = (%td, %td, %td), gc_above = (%td, %td, %td):\n", gc_below[0], gc_below[1], gc_below[2], gc_above[0], gc_above[1], gc_above[2]); 
  pfft_printf(comm_cart_2d, "maxerror = %6.2e;\n", err);


  /* free mem and finalize */
  pfft_destroy_gcplan(ths);
  MPI_Comm_free(&comm_cart_2d);
  pfft_free(data);
  MPI_Finalize();
  return 0;
}

