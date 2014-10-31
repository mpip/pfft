#include <complex.h>
#include <pfft.h>

static void init_parameters(
    int argc, char **argv,
    int *np, ptrdiff_t *n, int *iter,
    int *inplace, int* patience
    )
{
  pfft_get_args(argc, argv, "-pfft_np", 2, PFFT_INT, np);
  pfft_get_args(argc, argv, "-pfft_n", 3, PFFT_PTRDIFF_T, n);
  pfft_get_args(argc, argv, "-pfft_iter", 1, PFFT_INT, iter);
  pfft_get_args(argc, argv, "-pfft_ip", 1, PFFT_INT, inplace);
  pfft_get_args(argc, argv, "-pfft_pat", 1, PFFT_INT, patience);
  
}


int main(int argc, char **argv)
{
  int np[2];
  ptrdiff_t n[3];
  ptrdiff_t alloc_local;
  ptrdiff_t local_ni[3], local_i_start[3];
  ptrdiff_t local_no[3], local_o_start[3];
  double err;
  pfft_complex *in, *out;
  pfft_plan plan_forw=NULL, plan_back=NULL;
  MPI_Comm comm_cart_2d;
  double time;
  pfft_timer timer_forw, timer_back;
  unsigned pfft_opt_flag;

  /* setup default parameters */
  int iter = 10, inplace = 0, patience = 0;  
  
  /* Set size of FFT and process mesh */
  n[0] = n[1] = n[2] = 16;
  np[0] = 2; np[1] = 2;
 
  /* Initialize MPI and PFFT */
  MPI_Init(&argc, &argv);
  pfft_init();

  /* read parameters from command line */
  init_parameters(argc, argv, np, n, &iter, &inplace, &patience);

  /* setup FFTWs planing depth */  
  switch(patience){
    case 1: pfft_opt_flag = PFFT_MEASURE; break;
    case 2: pfft_opt_flag = PFFT_PATIENT; break;
    case 3: pfft_opt_flag = PFFT_EXHAUSTIVE; break;
    default: pfft_opt_flag = PFFT_ESTIMATE;
  }
  pfft_opt_flag |= PFFT_DESTROY_INPUT;
  
  /* Create two-dimensional process grid of size np[0] x np[1], if possible */
  if( pfft_create_procmesh_2d(MPI_COMM_WORLD, np[0], np[1], &comm_cart_2d) ){
    pfft_fprintf(MPI_COMM_WORLD, stderr, "Error: Procmesh %d x %d requires MPI launch with %d processes.\n",
        np[0], np[1], np[0]*np[1]);
    MPI_Finalize();
    MPI_Finalize();
    return 1;
  }
  
  /* Get parameters of data distribution */
  alloc_local = pfft_local_size_dft_3d(n, comm_cart_2d, PFFT_TRANSPOSED_OUT,
      local_ni, local_i_start, local_no, local_o_start);

  /* Allocate memory */
  in = pfft_alloc_complex(alloc_local);
  out = (inplace) ? in : pfft_alloc_complex(alloc_local);

  /* We often want to scale large FFTs, which do not fit on few processes. */
  if( (in == NULL) || (out == NULL)){
    fprintf(stderr, "!!! Error: Not enough memory to allocate input/output arrays !!!\n");
    MPI_Finalize();
    MPI_Finalize();
    return 1;
  }

  
  /* Plan parallel forward FFT */
  time = -MPI_Wtime();
  plan_forw = pfft_plan_dft_3d(
      n, in, out, comm_cart_2d, PFFT_FORWARD, PFFT_TRANSPOSED_OUT| pfft_opt_flag);
  time += MPI_Wtime();
//  printf("time for forw planing: %.2e\n", time);
  
  /* Plan parallel backward FFT */
  time = -MPI_Wtime();
  plan_back = pfft_plan_dft_3d(
      n, out, in, comm_cart_2d, PFFT_BACKWARD, PFFT_TRANSPOSED_IN| pfft_opt_flag);
  time += MPI_Wtime();
//  printf("time for back planing: %.2e\n", time);

  /* Initialize input with random numbers */
  pfft_init_input_complex_3d(n, local_ni, local_i_start,
      in);
  
  for(int t=0; t<iter; t++){
    /* execute parallel forward FFT */
    pfft_execute(plan_forw);

  /* clear the old input */
  pfft_clear_input_complex_3d(n, local_ni, local_i_start,
      in);
  
    /* execute parallel backward FFT */
    pfft_execute(plan_back);
  }
 
  /* check individual timers for workbalance */
  timer_forw = pfft_get_timer(plan_forw);
//    printf("timer_forw->whole = %.2e\n", timer_forw->whole);
  pfft_destroy_timer(timer_forw);
  timer_back = pfft_get_timer(plan_back);
//  printf("timer_back->whole = %.2e\n", timer_back->whole);
  pfft_destroy_timer(timer_back);

  /* read out PFFT timers */ 
  pfft_print_average_timer_adv(plan_forw, comm_cart_2d);
  pfft_print_average_timer_adv(plan_back, comm_cart_2d);
  if(inplace){
    pfft_write_average_timer_adv(plan_forw, "measure_forw_inplace.m", comm_cart_2d);
    pfft_write_average_timer_adv(plan_back, "measure_back_inplace.m", comm_cart_2d);
  } else {
    pfft_write_average_timer_adv(plan_forw, "measure_forw_outofplace.m", comm_cart_2d);
    pfft_write_average_timer_adv(plan_back, "measure_back_outofplace.m", comm_cart_2d);
  }
  
  /* Scale data */
  for(int t=0; t<iter; t++)
    for(ptrdiff_t l=0; l < local_ni[0] * local_ni[1] * local_ni[2]; l++)
      in[l] /= (n[0]*n[1]*n[2]);

  /* Print error of back transformed data */
  MPI_Barrier(MPI_COMM_WORLD);
  err = pfft_check_output_complex_3d(n, local_ni, local_i_start, in, comm_cart_2d);
  pfft_printf(comm_cart_2d, "Error after one forward and backward trafo of size n=(%td, %td, %td):\n", n[0], n[1], n[2]); 
  pfft_printf(comm_cart_2d, "maxerror = %6.2e;\n", err);
  
  /* free mem and finalize */
  pfft_destroy_plan(plan_forw);
  pfft_destroy_plan(plan_back);
  MPI_Comm_free(&comm_cart_2d);
  pfft_free(in); if(!inplace) pfft_free(out);
  MPI_Finalize();
  return 0;
}
