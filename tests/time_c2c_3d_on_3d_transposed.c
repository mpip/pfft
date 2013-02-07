#include <complex.h>
#include <pfft.h>

static void measure_pfft(
    const ptrdiff_t *n, MPI_Comm comm_cart_3d,
    unsigned pfft_opt_flags, unsigned verbose);
static void measure_fftw(
    const ptrdiff_t *n, unsigned parallel,
    unsigned fftw_opt_flags, unsigned verbose);
static void loop_pfft_tests(
    ptrdiff_t *n, MPI_Comm comm, unsigned verbose);
static void loop_fftw_tests(
    ptrdiff_t *n, unsigned parallel, unsigned verbose);
static void init_parameters(
    int argc, char **argv,
    ptrdiff_t *n, int *np, unsigned *verbose);

int main(int argc, char **argv)
{
  int np[3];
  ptrdiff_t n[3];
  unsigned parallel, verbose;
  MPI_Comm comm_cart_1d, comm_cart_2d, comm_cart_3d;
  
  /* Set size of FFT and process mesh */
  n[0] = 128; n[1] = 128; n[2] = 128;
//   n[0] = 2; n[1] = 1; n[2] = 1;
  np[0] = 1; np[1] = 1; np[2] = 1;
  verbose = 0;

  /* Initialize MPI and PFFT */
  MPI_Init(&argc, &argv);
  pfft_init();
 
  /* set parameters by command line */
  init_parameters(argc, argv, n, np, &verbose);

  pfft_printf(MPI_COMM_WORLD, "******************************************************************************************************\n");
  pfft_printf(MPI_COMM_WORLD, "* Computation of parallel forward and backward FFT\n");
  pfft_printf(MPI_COMM_WORLD, "* for n[0] x n[1] x n[2] = %td x %td x %td Fourier coefficients (change with -pfft_n * * *)\n", n[0], n[1], n[2]);
  pfft_printf(MPI_COMM_WORLD, "* on  np[0] x np[1] x np[2] = %td x %td x %td processes (change with -pfft_np * * *)\n", np[0], np[1], np[2]);
  pfft_printf(MPI_COMM_WORLD, "*******************************************************************************************************\n\n");

  /* Create three-dimensional process grid of size np[0] x np[1] x np[2], if possible */
  if( pfft_create_procmesh(3, MPI_COMM_WORLD, np, &comm_cart_3d) ){
    pfft_fprintf(MPI_COMM_WORLD, stderr, "Error: Procmesh of size %d x %d x %d does not fit to number of allocated processes.\n", np[0], np[1], np[2]);
    pfft_fprintf(MPI_COMM_WORLD, stderr, "       Please allocate %d processes (mpiexec -np %d ...) or change the procmesh (with -pfft_np * * *).\n", np[0]*np[1]*np[2], np[0]*np[1]*np[2]);
    MPI_Finalize();
    return 1;
  }
  
  pfft_printf(MPI_COMM_WORLD, "PFFT runtimes (3d data decomposition):\n");
  loop_pfft_tests(n, comm_cart_3d, verbose);
  MPI_Comm_free(&comm_cart_3d);

  /* run 2d-data decomposition if possible */
  if( ((np[0]==1) + (np[1]==1) + (np[2]==1)) >= 1 ){
    if(np[0]==1){ np[0] = np[1]; np[1] = np[2]; np[2] = 1; }

    pfft_create_procmesh(2, MPI_COMM_WORLD, np, &comm_cart_2d);
    pfft_printf(MPI_COMM_WORLD, "\nPFFT runtimes (2d data decomposition):\n");
    loop_pfft_tests(n, comm_cart_2d, verbose);
    MPI_Comm_free(&comm_cart_2d);
  }

  /* run 1d-data decomposition if possible */
  if( ((np[0]==1) + (np[1]==1) + (np[2]==1)) >= 2 ){
    if(np[0]==1){ np[0] = np[1]; np[1] = np[2]; np[2] = 1; }

    pfft_create_procmesh(1, MPI_COMM_WORLD, np, &comm_cart_1d);
    pfft_printf(MPI_COMM_WORLD, "\nPFFT runtimes (1d data decomposition):\n");
    loop_pfft_tests(n, comm_cart_1d, verbose);
    MPI_Comm_free(&comm_cart_1d);

    pfft_printf(MPI_COMM_WORLD, "\nFFTW_MPI runtimes (1d data decomposition):\n");
    loop_fftw_tests(n, parallel=1, verbose);
  }

  /* run serial if possible */
  if( np[0]*np[1]*np[2] == 1 ){
    pfft_printf(MPI_COMM_WORLD, "\nserial FFTW runtimes (no data decomposition at all):\n");
    loop_fftw_tests(n, parallel=0, verbose);
  }

  /* free mem and finalize */
  MPI_Finalize();
  return 0;
}

static void loop_pfft_tests(
    ptrdiff_t *n, MPI_Comm comm, unsigned verbose
    )
{
  
  unsigned tune, measure, destroy;

  destroy = 0; 
  for(int k=0; k<2; k++){
    measure = PFFT_ESTIMATE;
    for(int l=0; l<2; l++){
      tune = PFFT_NO_TUNE;
      for(int m=0; m<2; m++){
        measure_pfft(n, comm, tune | measure | destroy, verbose);
        tune = PFFT_TUNE;
      }
      measure = PFFT_MEASURE;
    }
    destroy = PFFT_DESTROY_INPUT;
  }
}

static void loop_fftw_tests(
    ptrdiff_t *n, unsigned parallel, unsigned verbose
    )
{
  unsigned measure, destroy;

  destroy = 0;
  for(int k=0; k<2; k++){
    measure = FFTW_ESTIMATE;
    for(int l=0; l<2; l++){
      measure_fftw(n, parallel, measure | destroy, verbose);
      measure = FFTW_MEASURE;
    }
    destroy = FFTW_DESTROY_INPUT;
  }
}

static void measure_pfft(
    const ptrdiff_t *n, MPI_Comm comm_cart,
    unsigned pfft_opt_flags, unsigned verbose
    )
{
  ptrdiff_t alloc_local;
  ptrdiff_t local_ni[3], local_i_start[3];
  ptrdiff_t local_no[3], local_o_start[3];
  double err, timer[4];
  pfft_complex *in, *out;
  pfft_plan plan_forw=NULL, plan_back=NULL;

  /* Get parameters of data distribution */
  alloc_local = pfft_local_size_dft_3d(n, comm_cart, PFFT_TRANSPOSED_OUT,
      local_ni, local_i_start, local_no, local_o_start);

  /* Allocate memory */
  in  = pfft_alloc_complex(alloc_local);
  out = pfft_alloc_complex(alloc_local);

  /* Plan parallel forward FFT */
  timer[0] = -MPI_Wtime();
  plan_forw = pfft_plan_dft_3d(
      n, in, out, comm_cart, PFFT_FORWARD, PFFT_TRANSPOSED_OUT| pfft_opt_flags);
  timer[0] += MPI_Wtime();
  
  /* Plan parallel backward FFT */
  timer[1] = -MPI_Wtime();
  plan_back = pfft_plan_dft_3d(
      n, out, in, comm_cart, PFFT_BACKWARD, PFFT_TRANSPOSED_IN| pfft_opt_flags);
  timer[1] += MPI_Wtime();

  /* Initialize input with random numbers */
  pfft_init_input_c2c_3d(n, local_ni, local_i_start,
      in);

  if(verbose)
    pfft_apr_complex_3d(in, local_ni, local_i_start, "PFFT Input", comm_cart);

  /* execute parallel forward FFT */
  timer[2] = -MPI_Wtime();
  pfft_execute(plan_forw);
  timer[2] += MPI_Wtime();

  if(verbose)
    pfft_apr_complex_3d(out, local_no, local_o_start, "PFFT Output", comm_cart);
  
  /* execute parallel backward FFT */
  timer[3] = -MPI_Wtime();
  pfft_execute(plan_back);
  timer[3] += MPI_Wtime();
  
  /* Scale data */
  for(ptrdiff_t l=0; l < local_ni[0] * local_ni[1] * local_ni[2]; l++)
    in[l] /= (n[0]*n[1]*n[2]);

  if(verbose)
    pfft_apr_complex_3d(in, local_ni, local_i_start, "Inputs after forward and backward PFFT", comm_cart);

  /* Print error of back transformed data */
  MPI_Barrier(MPI_COMM_WORLD);
  err = pfft_check_output_c2c_3d(n, local_ni, local_i_start, in, comm_cart);

  /* Print optimization flags */
  pfft_printf(comm_cart, "Flags: ");
  if(pfft_opt_flags & PFFT_TUNE)
    pfft_printf(comm_cart, "PFFT_TUNE");
  else
    pfft_printf(comm_cart, "PFFT_NO_TUNE");

  pfft_printf(comm_cart, ", ");

  if(pfft_opt_flags & PFFT_ESTIMATE)
    pfft_printf(comm_cart, "PFFT_ESTIMATE");
  else if(pfft_opt_flags & PFFT_PATIENT)
    pfft_printf(comm_cart, "PFFT_PATIENT");
  else if(pfft_opt_flags & PFFT_EXHAUSTIVE)
    pfft_printf(comm_cart, "PFFT_EXHAUSTIVE");
  else
    pfft_printf(comm_cart, "PFFT_MEASURE");

  pfft_printf(comm_cart, ", ");

  if(pfft_opt_flags & PFFT_DESTROY_INPUT)
    pfft_printf(comm_cart, "PFFT_DESTROY_INPUT");
  else
    pfft_printf(comm_cart, "PFFT_PRESERVE_INPUT");

  pfft_printf(comm_cart, "\n");

//   pfft_printf(comm_cart, "Error after one forward and backward trafo of size n=(%td, %td, %td):\n", n[0], n[1], n[2]); 
//   pfft_printf(comm_cart, "maxerror = %6.2e;\n", err);
//   pfft_printf(comm_cart, "Tuning and execution time of forward/backward trafo of size n=(%td, %td, %td):\n", n[0], n[1], n[2]); 
  pfft_printf(comm_cart, "tune_forw = %6.2e; tune_back = %6.2e, exec_forw = %6.2e, exec_back = %6.2e, error = %6.2e\n", timer[0], timer[1], timer[2], timer[3], err);

  /* free mem and finalize */
  pfft_destroy_plan(plan_forw);
  pfft_destroy_plan(plan_back);
  pfft_free(in); pfft_free(out);
}

static void measure_fftw(
    const ptrdiff_t *n, unsigned parallel,
    unsigned fftw_opt_flags, unsigned verbose
    )
{
  ptrdiff_t alloc_local;
  ptrdiff_t local_ni[3], local_i_start[3];
  ptrdiff_t local_no[3], local_o_start[3];
  double err, timer[4];
  fftw_complex *in, *out;
  fftw_plan plan_forw=NULL, plan_back=NULL;

  for(int t=0; t<3; t++){
    local_ni[t] = local_no[t] = n[t];
    local_i_start[t] = local_o_start[t] = 0;
  }

  if(parallel){
    /* transposed output */
    local_no[0] = n[1]; local_no[1] = n[0]; local_no[2] = n[2];

    /* Get parameters of data distribution */
    alloc_local = fftw_mpi_local_size_3d_transposed(n[0], n[1], n[2], MPI_COMM_WORLD,
        &local_ni[0], &local_i_start[0], &local_no[0], &local_o_start[0]);
  } else
    alloc_local = n[0]*n[1]*n[2];

  /* Allocate memory */
  in  = fftw_alloc_complex(alloc_local);
  out = fftw_alloc_complex(alloc_local);

  /* Plan parallel forward FFT */
  if(parallel){
    timer[0] = -MPI_Wtime();
    plan_forw = fftw_mpi_plan_dft_3d(
        n[0], n[1], n[2], in, out, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MPI_TRANSPOSED_OUT| fftw_opt_flags);
    timer[0] += MPI_Wtime();
  } else {
    timer[0] = -MPI_Wtime();
    plan_forw = fftw_plan_dft_3d(
        n[0], n[1], n[2], in, out, FFTW_FORWARD, fftw_opt_flags);
    timer[0] += MPI_Wtime();
  }
  
  /* Plan parallel backward FFT */
  if(parallel){
    timer[1] = -MPI_Wtime();
    plan_back = fftw_mpi_plan_dft_3d(
        n[0], n[1], n[2], out, in, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MPI_TRANSPOSED_IN| fftw_opt_flags);
    timer[1] += MPI_Wtime();
  } else {
    timer[1] = -MPI_Wtime();
    plan_back = fftw_plan_dft_3d(
        n[0], n[1], n[2], out, in, FFTW_BACKWARD, fftw_opt_flags);
    timer[1] += MPI_Wtime();
  }

  /* Initialize input with random numbers */
  pfft_init_input_c2c_3d(n, local_ni, local_i_start,
      in);

  if(verbose)
    pfft_apr_complex_3d(in, local_ni, local_i_start, "PFFT Input", MPI_COMM_WORLD);

  /* execute parallel forward FFT */
  timer[2] = -MPI_Wtime();
  fftw_execute(plan_forw);
  timer[2] += MPI_Wtime();

  if(verbose)
    pfft_apr_complex_3d(out, local_no, local_o_start, "PFFT Output", MPI_COMM_WORLD);
  
  /* execute parallel backward FFT */
  timer[3] = -MPI_Wtime();
  fftw_execute(plan_back);
  timer[3] += MPI_Wtime();
  
  /* Scale data */
  for(ptrdiff_t l=0; l < local_ni[0] * local_ni[1] * local_ni[2]; l++)
    in[l] /= (n[0]*n[1]*n[2]);

  if(verbose)
    pfft_apr_complex_3d(in, local_ni, local_i_start, "Inputs after forward and backward PFFT", MPI_COMM_WORLD);
  
  /* Print error of back transformed data */
  MPI_Barrier(MPI_COMM_WORLD);
  err = pfft_check_output_c2c_3d(n, local_ni, local_i_start, in, MPI_COMM_WORLD);

  /* Print optimization flags */
  pfft_printf(MPI_COMM_WORLD, "Flags: ");
  if(fftw_opt_flags & FFTW_ESTIMATE)
    pfft_printf(MPI_COMM_WORLD, "FFTW_ESTIMATE");
  else if(fftw_opt_flags & FFTW_PATIENT)
    pfft_printf(MPI_COMM_WORLD, "FFTW_PATIENT");
  else if(fftw_opt_flags & FFTW_EXHAUSTIVE)
    pfft_printf(MPI_COMM_WORLD, "FFTW_EXHAUSTIVE");
  else
    pfft_printf(MPI_COMM_WORLD, "FFTW_MEASURE");

  pfft_printf(MPI_COMM_WORLD, ", ");

  if(fftw_opt_flags & FFTW_DESTROY_INPUT)
    pfft_printf(MPI_COMM_WORLD, "FFTW_DESTROY_INPUT");
  else
    pfft_printf(MPI_COMM_WORLD, "FFTW_PRESERVE_INPUT");

  pfft_printf(MPI_COMM_WORLD, "\n");

//   pfft_printf(MPI_COMM_WORLD, "Error after one forward and backward trafo of size n=(%td, %td, %td):\n", n[0], n[1], n[2]); 
//   pfft_printf(MPI_COMM_WORLD, "maxerror = %6.2e;\n", err);
//   pfft_printf(MPI_COMM_WORLD, "Tuning and execution time of forward/backward trafo of size n=(%td, %td, %td):\n", n[0], n[1], n[2]); 
  pfft_printf(MPI_COMM_WORLD, "tune_forw = %6.2e; tune_back = %6.2e, exec_forw = %6.2e, exec_back = %6.2e, error = %6.2e\n", timer[0], timer[1], timer[2], timer[3], err);

  /* free mem and finalize */
  fftw_destroy_plan(plan_forw);
  fftw_destroy_plan(plan_back);
  fftw_free(in); fftw_free(out);
}

static void init_parameters(
    int argc, char **argv,
    ptrdiff_t *n, int *np, unsigned *verbose
    )
{
  pfft_get_args(argc, argv, "-pfft_n", 3, PFFT_PTRDIFF_T, n);
  pfft_get_args(argc, argv, "-pfft_np", 3, PFFT_INT, np);
  pfft_get_args(argc, argv, "-pfft_verbose", 1, PFFT_UNSIGNED, verbose);
}

