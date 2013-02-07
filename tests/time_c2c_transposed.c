#include <complex.h>
#include <pfft.h>

static void init_parameters(
    int argc, char **argv,
    ptrdiff_t *n, int *np, int *loops, int *inplace,
    unsigned *opt_flag, unsigned *tune_flag, unsigned *destroy_flag);
static void measure_pfft(
    const ptrdiff_t *n, int *np, MPI_Comm comm,
    int loops, int inplace, unsigned pfft_opt_flags);

int main(int argc, char **argv)
{
  int np[2], inplace, loops;
  ptrdiff_t n[3];
  unsigned opt, tune, destroy_input;
  
  /* Set size of FFT and process mesh */
  n[0] = 29; n[1] = 27; n[2] = 31;
  np[0] = 2; np[1] = 2;
  inplace = 0;
  opt     = PFFT_ESTIMATE;
  tune    = PFFT_NO_TUNE;
  destroy_input = PFFT_PRESERVE_INPUT;
  loops   = 1;
  
  /* Initialize MPI and PFFT */
  MPI_Init(&argc, &argv);
  pfft_init();

  /* set parameters by command line */
  init_parameters(argc, argv, n, np, &loops, &inplace, &opt, &tune, &destroy_input);

  measure_pfft(n, np, MPI_COMM_WORLD, loops, inplace, opt | tune | destroy_input);

  MPI_Finalize();
  return 0;
} 


static void measure_pfft(
    const ptrdiff_t *n, int *np, MPI_Comm comm,
    int loops, int inplace, unsigned pfft_opt_flags
    )
{
  ptrdiff_t alloc_local;
  ptrdiff_t local_ni[3], local_i_start[3];
  ptrdiff_t local_no[3], local_o_start[3];
  double err=0.0, timer[4], max_timer[4];
  pfft_complex *in, *out;
  pfft_plan plan_forw=NULL, plan_back=NULL;
  MPI_Comm comm_cart_2d;

  /* Create two-dimensional process grid of size np[0] x np[1], if possible */
  if( pfft_create_procmesh_2d(comm, np[0], np[1], &comm_cart_2d) ){
    pfft_fprintf(comm, stderr, "Error: This test file only works with %d processes.\n", np[0]*np[1]);
    return;
  }
  
  /* Get parameters of data distribution */
  alloc_local = pfft_local_size_dft_3d(n, comm_cart_2d, PFFT_TRANSPOSED_OUT,
      local_ni, local_i_start, local_no, local_o_start);

  /* Allocate memory */
  in  = pfft_alloc_complex(alloc_local);
  out = (inplace) ? in : pfft_alloc_complex(alloc_local);

  /* Plan parallel forward FFT */
  timer[0] = -MPI_Wtime();
  plan_forw = pfft_plan_dft_3d(
      n, in, out, comm_cart_2d, PFFT_FORWARD, PFFT_TRANSPOSED_OUT| pfft_opt_flags);
  timer[0] += MPI_Wtime();
  
  /* Plan parallel backward FFT */
  timer[1] = -MPI_Wtime();
  plan_back = pfft_plan_dft_3d(
      n, out, in, comm_cart_2d, PFFT_BACKWARD, PFFT_TRANSPOSED_IN| pfft_opt_flags);
  timer[1] += MPI_Wtime();

  /* Initialize input with random numbers */
  pfft_init_input_c2c_3d(n, local_ni, local_i_start,
      in);

  pfft_reset_timer(plan_forw);
  pfft_reset_timer(plan_back);

  timer[2] = timer[3] = 0;
  for(int t=0; t<loops; t++){
    /* execute parallel forward FFT */
    timer[2] -= MPI_Wtime();
    pfft_execute(plan_forw);
    timer[2] += MPI_Wtime();
    
    /* execute parallel backward FFT */
    timer[3] -= MPI_Wtime();
    pfft_execute(plan_back);
    timer[3] += MPI_Wtime();
    
    /* Scale data */
    for(ptrdiff_t l=0; l < local_ni[0] * local_ni[1] * local_ni[2]; l++)
      in[l] /= (n[0]*n[1]*n[2]);
  }
  timer[2] /= loops;
  timer[3] /= loops;

  /* Print pfft timer */
  pfft_print_average_timer_adv(plan_forw, comm_cart_2d);
  pfft_print_average_timer_adv(plan_back, comm_cart_2d);

  /* Print optimization flags */
  pfft_printf(comm_cart_2d, "\nFlags = ");
  if(pfft_opt_flags & PFFT_TUNE)
    pfft_printf(comm_cart_2d, "PFFT_TUNE");
  else
    pfft_printf(comm_cart_2d, "PFFT_NO_TUNE");

  pfft_printf(comm_cart_2d, " | ");

  if(pfft_opt_flags & PFFT_ESTIMATE)
    pfft_printf(comm_cart_2d, "PFFT_ESTIMATE");
  else if(pfft_opt_flags & PFFT_PATIENT)
    pfft_printf(comm_cart_2d, "PFFT_PATIENT");
  else if(pfft_opt_flags & PFFT_EXHAUSTIVE)
    pfft_printf(comm_cart_2d, "PFFT_EXHAUSTIVE");
  else
    pfft_printf(comm_cart_2d, "PFFT_MEASURE");

  pfft_printf(comm_cart_2d, " | ");

  if(pfft_opt_flags & PFFT_DESTROY_INPUT)
    pfft_printf(comm_cart_2d, "PFFT_DESTROY_INPUT");
  else
    pfft_printf(comm_cart_2d, "PFFT_PRESERVE_INPUT");

  pfft_printf(comm_cart_2d, "\n");


  /* Print error of back transformed data */
  err = pfft_check_output_c2c_3d(n, local_ni, local_i_start, in, comm_cart_2d);
  pfft_printf(comm_cart_2d, "Run %d loops of ", loops);
  if(inplace)
    pfft_printf(comm_cart_2d, "in-place");
  else
    pfft_printf(comm_cart_2d, "out-of-place");
  pfft_printf(comm_cart_2d, " forward and backward trafo of size n=(%td, %td, %td):\n", n[0], n[1], n[2]); 

  MPI_Reduce(&timer, &max_timer, 4, MPI_DOUBLE, MPI_MAX, 0, comm_cart_2d);
  pfft_printf(comm_cart_2d, "tune_forw = %6.2e; tune_back = %6.2e, exec_forw = %6.2e, exec_back = %6.2e, error = %6.2e\n", max_timer[0], max_timer[1], max_timer[2], max_timer[3], err);

  /* free mem and finalize */
  pfft_destroy_plan(plan_forw);
  pfft_destroy_plan(plan_back);
  MPI_Comm_free(&comm_cart_2d);
  if(in != out) pfft_free(out);
  pfft_free(in);
}

static void init_parameters(
    int argc, char **argv,
    ptrdiff_t *n, int *np, int *loops, int *inplace,
    unsigned *opt_flag, unsigned *tune_flag, unsigned *destroy_input_flag 
    )
{
  int opt=0, tune=0, destroy_input=0;

  pfft_get_args(argc, argv, "-pfft_n", 3, PFFT_PTRDIFF_T, n);
  pfft_get_args(argc, argv, "-pfft_np", 2, PFFT_INT, np);
  pfft_get_args(argc, argv, "-pfft_loops", 1, PFFT_INT, loops);
  pfft_get_args(argc, argv, "-pfft_ip", 1, PFFT_INT, inplace);
  pfft_get_args(argc, argv, "-pfft_opt", 1, PFFT_INT, &opt);
  pfft_get_args(argc, argv, "-pfft_tune", 1, PFFT_INT, &tune);
  pfft_get_args(argc, argv, "-pfft_di", 1, PFFT_INT, &destroy_input);


  switch(opt){
    case 1: *opt_flag = PFFT_MEASURE; break;
    case 2: *opt_flag = PFFT_PATIENT; break;
    case 3: *opt_flag = PFFT_EXHAUSTIVE; break;
  }

  if(destroy_input)
    *destroy_input_flag = PFFT_DESTROY_INPUT;

  if(tune)
    *tune_flag = PFFT_TUNE;
}

