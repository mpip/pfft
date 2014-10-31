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
  MPI_Comm comm_cart_1d;
  double time;
  pfft_timer timer_forw, timer_back;
  unsigned fftw_flag, pfft_flag;

  /* setup default parameters */
  int iter = 10, inplace = 0, patience = 0;  
  
  /* Set size of FFT and process mesh */
  n[0] = n[1] = n[2] = 16;
  np[0] = 4; np[1] = 1;
 
  /* Initialize MPI and PFFT */
  MPI_Init(&argc, &argv);
  pfft_init();

  /* read parameters from command line */
  init_parameters(argc, argv, np, n, &iter, &inplace, &patience);

  /* setup FFTWs planing depth */  
  switch(patience){
    case 1: pfft_flag = PFFT_MEASURE; 
            fftw_flag = FFTW_MEASURE; break;
    case 2: pfft_flag = PFFT_PATIENT;
            fftw_flag = FFTW_PATIENT; break;
    case 3: pfft_flag = PFFT_EXHAUSTIVE;
            fftw_flag = FFTW_EXHAUSTIVE; break;
    default: pfft_flag = PFFT_ESTIMATE;
             fftw_flag = FFTW_ESTIMATE;
  }
  
  if(!inplace){
    pfft_flag |= PFFT_DESTROY_INPUT;
    fftw_flag |= FFTW_DESTROY_INPUT;
  }

  /* Create two-dimensional process grid of size np[0] x np[1], if possible */
  if( pfft_create_procmesh(1, MPI_COMM_WORLD, np, &comm_cart_1d) ){
    pfft_fprintf(MPI_COMM_WORLD, stderr, "Error: Procmesh %d x %d requires MPI launch with %d processes.\n",
        np[0], np[1], np[0]*np[1]);
    MPI_Finalize();
    MPI_Finalize();
    return 1;
  }
  if( np[1] != 1){
    pfft_fprintf(MPI_COMM_WORLD, stderr, "Error: Procmesh %d x %d is not one-dimensional.\n",
        np[0], np[1]);
    MPI_Finalize();
    MPI_Finalize();
    return 1;
  }
  
  /* Get parameters of data distribution */
  alloc_local = pfft_local_size_dft_3d(n, comm_cart_1d, PFFT_TRANSPOSED_OUT,
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
      n, in, out, comm_cart_1d, PFFT_FORWARD, PFFT_TRANSPOSED_OUT| pfft_flag);
  time += MPI_Wtime();
//  printf("time for forw planing: %.2e\n", time);
  
  /* Plan parallel backward FFT */
  time = -MPI_Wtime();
  plan_back = pfft_plan_dft_3d(
      n, out, in, comm_cart_1d, PFFT_BACKWARD, PFFT_TRANSPOSED_IN| pfft_flag);
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
  pfft_print_average_timer_adv(plan_forw, comm_cart_1d);
  pfft_print_average_timer_adv(plan_back, comm_cart_1d);
  if(inplace){
    pfft_write_average_timer_adv(plan_forw, "measure_forw_inplace.m", comm_cart_1d);
    pfft_write_average_timer_adv(plan_back, "measure_back_inplace.m", comm_cart_1d);
  } else {
    pfft_write_average_timer_adv(plan_forw, "measure_forw_outofplace.m", comm_cart_1d);
    pfft_write_average_timer_adv(plan_back, "measure_back_outofplace.m", comm_cart_1d);
  }
  
  /* Scale data */
  for(int t=0; t<iter; t++)
    for(ptrdiff_t l=0; l < local_ni[0] * local_ni[1] * local_ni[2]; l++)
      in[l] /= (n[0]*n[1]*n[2]);

  /* Print error of back transformed data */
  err = pfft_check_output_complex_3d(n, local_ni, local_i_start, in, comm_cart_1d);
  pfft_printf(comm_cart_1d, "Error after one forward and backward trafo of size n=(%td, %td, %td):\n", n[0], n[1], n[2]); 
  pfft_printf(comm_cart_1d, "maxerror = %6.2e;\n", err);
  


  FFTW(plan) pf, pb;
  pfft_complex *in1, *out1;
  int myrank;
  ptrdiff_t lni[3], lis[3], lno[3], los[3];
  double time_fftw[2], max_time_fftw[2];

  for(int t=0; t<3; t++){
    lni[t] = n[t]; lis[t] = 0;
    lno[t] = n[t]; los[t] = 0;
  }

  alloc_local = FFTW(mpi_local_size_3d_transposed)(n[0], n[1], n[2], comm_cart_1d,
      &lni[0], &lis[0], &lno[1], &los[1]);

  in1 = pfft_alloc_complex(alloc_local);
  out1 = (inplace) ? in1 : pfft_alloc_complex(alloc_local);
  
  pf = FFTW(mpi_plan_dft_3d)(n[0], n[1], n[2], in1, out1, comm_cart_1d,
      FFTW_FORWARD, FFTW_MPI_TRANSPOSED_OUT| fftw_flag);
  pb = FFTW(mpi_plan_dft_3d)(n[0], n[1], n[2], out1, in1, comm_cart_1d,
      FFTW_BACKWARD, FFTW_MPI_TRANSPOSED_IN| fftw_flag);
  
  pfft_init_input_complex_3d(n, lni, lis,
      in1);

  time_fftw[0] = time_fftw[1] = 0;
  for(int t=0; t<iter; t++){
    /* execute parallel forward FFT */
    time_fftw[0] -= MPI_Wtime();
    FFTW(execute)(pf);
    time_fftw[0] += MPI_Wtime();
  
    /* execute parallel backward FFT */
    time_fftw[1] -= MPI_Wtime();
    FFTW(execute)(pb);
    time_fftw[1] += MPI_Wtime();
  }
  
  /* Scale data */
  for(int t=0; t<iter; t++)
    for(ptrdiff_t l=0; l < lni[0] * lni[1] * lni[2]; l++)
      in1[l] /= (n[0]*n[1]*n[2]);

  MPI_Reduce(time_fftw, max_time_fftw, 2, MPI_DOUBLE, MPI_MAX, 0 ,comm_cart_1d);
  MPI_Comm_rank(comm_cart_1d, &myrank);
  if(myrank==0)
    printf("fftw_forw = %.2e, fftw_back = %.2e\n", max_time_fftw[0]/iter, max_time_fftw[1]/iter);
 
  err = pfft_check_output_complex_3d(n, lni, lis, in1, comm_cart_1d);
  pfft_printf(comm_cart_1d, "Error after several forward and backward FFTWs of size n=(%td, %td, %td):\n", n[0], n[1], n[2]); 
  pfft_printf(comm_cart_1d, "maxerror = %6.2e;\n", err);
  
  /* free mem and finalize */
  pfft_destroy_plan(plan_forw);
  pfft_destroy_plan(plan_back);
  MPI_Comm_free(&comm_cart_1d);
  pfft_free(in); if(!inplace) pfft_free(out);
 
  /* free mem and finalize */
  FFTW(destroy_plan)(pf);
  FFTW(destroy_plan)(pb);
  pfft_free(in1); if(!inplace) pfft_free(out1);
  
  MPI_Finalize();
  return 0;
}
