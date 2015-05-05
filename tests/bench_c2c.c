#include <complex.h>
#include <pfft.h>

static void measure_pfft(
    const ptrdiff_t *n, MPI_Comm comm_cart_3d, int loops,
    unsigned pfft_opt_flags, int transposed, int inplace, int verbose,
    int print_timer);
static void measure_fftw(
    const ptrdiff_t *n, int parallel, int loops,
    unsigned fftw_opt_flags, int transposed, int inplace, int verbose);
static void loop_pfft_tests(
    ptrdiff_t *n, MPI_Comm comm, int loops,
    unsigned pfft_flags, int transposed, int inplace, int verbose, int cmp_flags,
    int print_timer);
static void loop_fftw_tests(
    ptrdiff_t *n, int parallel, int loops, int transposed, int inplace, int verbose);
static void init_parameters(
    int argc, char **argv,
    ptrdiff_t *n, int *np,
    unsigned *pfft_flags,
    int *loops, int *transposed, int *verbose, int *inplace,
    int *cmp_fftw, int *cmp_decomp, int *cmp_flags,
    int *print_timer);

int main(int argc, char **argv)
{
  int parallel;
  MPI_Comm comm_cart_1d, comm_cart_2d, comm_cart_3d;
  
  /* Set size of FFT and process mesh */
  ptrdiff_t n[3]        = {32,32,32};
  int       np[3]       = {1,1,1};
  int       loops       = 1;
  int       verbose     = 0;
  int       inplace     = 0;
  int       cmp_fftw    = 0;
  int       cmp_decomp  = 0;
  int       cmp_flags   = 0;
  int       transposed  = 0;
  int       print_timer = 0;
  unsigned  pfft_flags  = 0;

  /* Initialize MPI and PFFT */
  MPI_Init(&argc, &argv);
  pfft_init();
 
  /* set parameters by command line */
  init_parameters(argc, argv, n, np, &pfft_flags, &loops, &transposed, &verbose, &inplace, &cmp_fftw, &cmp_decomp, &cmp_flags, &print_timer);

  /* Create three-dimensional process grid of size np[0] x np[1] x np[2], if possible */
  if( pfft_create_procmesh(3, MPI_COMM_WORLD, np, &comm_cart_3d) ){
    pfft_fprintf(MPI_COMM_WORLD, stderr, "Error: Procmesh of size %d x %d x %d does not fit to number of allocated processes.\n", np[0], np[1], np[2]);
    pfft_fprintf(MPI_COMM_WORLD, stderr, "       Please allocate %d processes (mpiexec -np %d ...) or change the procmesh (with -pfft_np * * *).\n", np[0]*np[1]*np[2], np[0]*np[1]*np[2]);
    MPI_Finalize();
    return 1;
  }

  int num_serial_dims = (np[0]==1) + (np[1]==1) + (np[2]==1); 

  if( cmp_decomp || num_serial_dims==0){
    pfft_printf(MPI_COMM_WORLD, "* PFFT runtimes (3d data decomposition):\n");
    loop_pfft_tests(n, comm_cart_3d, loops, pfft_flags, transposed, inplace, verbose, cmp_flags, print_timer);
    pfft_printf(MPI_COMM_WORLD, "\n");
    MPI_Comm_free(&comm_cart_3d);
  }

  /* run 2d-data decomposition if possible */
  if( num_serial_dims >= 1 ){
    if( cmp_decomp || num_serial_dims==1){
      /* move serial dims to the end */ 
      if(np[1]==1){ np[1] = np[2]; np[2] = 1; }
      if(np[0]==1){ np[0] = np[1]; np[1] = 1; }
      if(np[1]==1){ np[1] = np[2]; np[2] = 1; }
  
      if( pfft_create_procmesh(2, MPI_COMM_WORLD, np, &comm_cart_2d) )
        pfft_printf(MPI_COMM_WORLD, "Error in creation of 2d procmesh of size %d x %d\n", np[0], np[1]);
      pfft_printf(MPI_COMM_WORLD, "* PFFT runtimes (2d data decomposition):\n");
      loop_pfft_tests(n, comm_cart_2d, loops, pfft_flags, transposed, inplace, verbose, cmp_flags, print_timer);
      pfft_printf(MPI_COMM_WORLD, "\n");
      MPI_Comm_free(&comm_cart_2d);
    }
  }

  /* run 1d-data decomposition if possible */
  if( num_serial_dims >= 2 ){
    /* move serial dims to the end */ 
    if(np[1]==1){ np[1] = np[2]; np[2] = 1; }
    if(np[0]==1){ np[0] = np[1]; np[1] = 1; }

    if( pfft_create_procmesh(1, MPI_COMM_WORLD, np, &comm_cart_1d) )
      pfft_printf(MPI_COMM_WORLD, "Error in creation of 2d procmesh of size %d\n", np[0]);
    pfft_printf(MPI_COMM_WORLD, "* PFFT runtimes (1d data decomposition):\n");
    loop_pfft_tests(n, comm_cart_1d, loops, pfft_flags, transposed, inplace, verbose, cmp_flags, print_timer);
    pfft_printf(MPI_COMM_WORLD, "\n");
    MPI_Comm_free(&comm_cart_1d);

    if(cmp_fftw){
      pfft_printf(MPI_COMM_WORLD, "* FFTW_MPI runtimes (1d data decomposition):\n");
      loop_fftw_tests(n, parallel=1, loops, transposed, inplace, verbose);
    }
  }

  /* run serial if possible */
  if( np[0]*np[1]*np[2] == 1 ){
    if(cmp_fftw){
      pfft_printf(MPI_COMM_WORLD, "* serial FFTW runtimes (no data decomposition at all):\n");
      loop_fftw_tests(n, parallel=0, loops, transposed, inplace, verbose);
      pfft_printf(MPI_COMM_WORLD, "\n");
    }
  }

  /* free mem and finalize */
  MPI_Finalize();
  return 0;
}

static void loop_pfft_tests(
    ptrdiff_t *n, MPI_Comm comm, int loops,
    unsigned pfft_flags, int transposed, int inplace, int verbose, int cmp_flags,
    int print_timer
    )
{
  
  unsigned tune, measure, destroy;

  if(!cmp_flags){
    measure_pfft(n, comm, loops, pfft_flags, transposed, inplace, verbose, print_timer);
    return;
  }

  destroy = 0; 
  for(int k=0; k<2; k++){
    measure = PFFT_ESTIMATE;
    for(int l=0; l<2; l++){
      tune = PFFT_NO_TUNE;
      for(int m=0; m<2; m++){
        measure_pfft(n, comm, loops, tune | measure | destroy, transposed, inplace, verbose, print_timer);
        tune = PFFT_TUNE;
      }
      measure = PFFT_MEASURE;
    }
    if(inplace) break;
    destroy = PFFT_DESTROY_INPUT;
  }
}

static void loop_fftw_tests(
    ptrdiff_t *n, int parallel, int loops, int transposed, int inplace, int verbose
    )
{
  unsigned measure, destroy;

  destroy = 0;
  for(int k=0; k<2; k++){
    measure = FFTW_ESTIMATE;
    for(int l=0; l<2; l++){
      measure_fftw(n, parallel, loops, measure | destroy, transposed, inplace, verbose);
      measure = FFTW_MEASURE;
    }
    if(inplace) break;
    destroy = FFTW_DESTROY_INPUT;
  }
}

static void measure_pfft(
    const ptrdiff_t *n, MPI_Comm comm_cart, int loops,
    unsigned pfft_opt_flags, int transposed, int inplace, int verbose,
    int print_timer
    )
{
  ptrdiff_t alloc_local;
  ptrdiff_t local_ni[3], local_i_start[3];
  ptrdiff_t local_no[3], local_o_start[3];
  double err, timer[4];
  pfft_complex *in, *out;
  pfft_plan plan_forw=NULL, plan_back=NULL;
  unsigned tr_in  = (transposed) ? PFFT_TRANSPOSED_IN  : PFFT_TRANSPOSED_NONE;
  unsigned tr_out = (transposed) ? PFFT_TRANSPOSED_OUT : PFFT_TRANSPOSED_NONE;

  /* Get parameters of data distribution */
  alloc_local = pfft_local_size_dft_3d(n, comm_cart, tr_out,
      local_ni, local_i_start, local_no, local_o_start);

  /* Allocate memory */
  in  = pfft_alloc_complex(alloc_local);
  if(inplace) out = in;
  else        out = pfft_alloc_complex(alloc_local);

  /* Plan parallel forward FFT */
  timer[0] = -MPI_Wtime();
  plan_forw = pfft_plan_dft_3d(
      n, in, out, comm_cart, PFFT_FORWARD, tr_out | pfft_opt_flags);
  timer[0] += MPI_Wtime();
  
  /* Plan parallel backward FFT */
  timer[1] = -MPI_Wtime();
  plan_back = pfft_plan_dft_3d(
      n, out, in, comm_cart, PFFT_BACKWARD, tr_in | pfft_opt_flags);
  timer[1] += MPI_Wtime();

  /* Initialize input with random numbers */
  pfft_init_input_complex_3d(n, local_ni, local_i_start,
      in);

  if(verbose)
    pfft_apr_complex_3d(in, local_ni, local_i_start, "PFFT Input", comm_cart);

  /* execute parallel forward FFT */
  timer[2] = -MPI_Wtime();
  for(int t=0; t<loops; t++)
    pfft_execute(plan_forw);
  timer[2] += MPI_Wtime();

  /* clear the old input */
  if(!inplace){
    pfft_clear_input_complex_3d(n, local_ni, local_i_start,
        in);
  }

  if(verbose)
    pfft_apr_complex_3d(out, local_no, local_o_start, "PFFT Output", comm_cart);
  
  /* execute parallel backward FFT */
  timer[3] = -MPI_Wtime();
  for(int t=0; t<loops; t++)
    pfft_execute(plan_back);
  timer[3] += MPI_Wtime();
  
  /* Scale data */
  for(int t=0; t<loops; t++)
    for(ptrdiff_t l=0; l < local_ni[0] * local_ni[1] * local_ni[2]; l++)
      in[l] /= (n[0]*n[1]*n[2]);

  if(verbose)
    pfft_apr_complex_3d(in, local_ni, local_i_start, "Inputs after forward and backward PFFT", comm_cart);

  /* Print error of back transformed data */
  MPI_Barrier(MPI_COMM_WORLD);
  err = pfft_check_output_complex_3d(n, local_ni, local_i_start, in, comm_cart);

  /* Print optimization flags */
  pfft_printf(comm_cart, "Flags: ");
  if(pfft_opt_flags & PFFT_TUNE) pfft_printf(comm_cart, "PFFT_TUNE");
  else                           pfft_printf(comm_cart, "PFFT_NO_TUNE");
  pfft_printf(comm_cart, ", ");

  if(pfft_opt_flags & PFFT_ESTIMATE)        pfft_printf(comm_cart, "PFFT_ESTIMATE");
  else if(pfft_opt_flags & PFFT_PATIENT)    pfft_printf(comm_cart, "PFFT_PATIENT");
  else if(pfft_opt_flags & PFFT_EXHAUSTIVE) pfft_printf(comm_cart, "PFFT_EXHAUSTIVE");
  else                                      pfft_printf(comm_cart, "PFFT_MEASURE");
  pfft_printf(comm_cart, ", ");

  if(pfft_opt_flags & PFFT_DESTROY_INPUT)  pfft_printf(comm_cart, "PFFT_DESTROY_INPUT, ");
  if(pfft_opt_flags & PFFT_PRESERVE_INPUT) pfft_printf(comm_cart, "PFFT_PRESERVE_INPUT, ");
  pfft_printf(comm_cart, "\n");

  pfft_printf(comm_cart, "tune_forw = %6.2e; tune_back = %6.2e, exec_forw/loops = %6.2e, exec_back/loops = %6.2e\n", timer[0], timer[1], timer[2]/loops, timer[3]/loops);
  if(loops==1) pfft_printf(MPI_COMM_WORLD, "error = %6.2e\n", err);

  if(print_timer){
    pfft_printf(MPI_COMM_WORLD, "Output of internal PFFT timers:\n");
    pfft_print_average_timer_adv(plan_forw, comm_cart);
    pfft_print_average_timer_adv(plan_back, comm_cart);
  }

  /* free mem and finalize */
  pfft_destroy_plan(plan_forw);
  pfft_destroy_plan(plan_back);
  pfft_free(in); if(!inplace) pfft_free(out);
}

static void measure_fftw(
    const ptrdiff_t *n, int parallel, int loops,
    unsigned fftw_opt_flags, int transposed, int inplace, int verbose
    )
{
  ptrdiff_t alloc_local;
  ptrdiff_t local_ni[3], local_i_start[3];
  ptrdiff_t local_no[3], local_o_start[3];
  double err, timer[4];
  fftw_complex *in, *out;
  fftw_plan plan_forw=NULL, plan_back=NULL;
  unsigned tr_in  = (transposed) ? FFTW_MPI_TRANSPOSED_IN  : 0;
  unsigned tr_out = (transposed) ? FFTW_MPI_TRANSPOSED_OUT : 0;

  for(int t=0; t<3; t++){
    local_ni[t] = local_no[t] = n[t];
    local_i_start[t] = local_o_start[t] = 0;
  }

  if(parallel){
    /* transposed output */
    if(transposed) {local_no[0] = n[1]; local_no[1] = n[0]; local_no[2] = n[2];}
    else           {local_no[0] = n[0]; local_no[1] = n[1]; local_no[2] = n[2];}

    /* Get parameters of data distribution */
    ptrdiff_t lni, lis, lno, los;
    alloc_local = fftw_mpi_local_size_3d_transposed(n[0], n[1], n[2], MPI_COMM_WORLD,
        &lni, &lis, &lno, &los);

    local_ni[0] = lni; local_i_start[0] = lis;
    if(transposed) {local_no[0] = lno; local_o_start[0] = los;}
    else           {local_no[0] = lni; local_o_start[0] = lis;}

  } else
    alloc_local = n[0]*n[1]*n[2];

  /* Allocate memory */
  in  = fftw_alloc_complex(alloc_local);
  if(inplace) out = in;
  else        out = fftw_alloc_complex(alloc_local);

  /* Plan parallel forward FFT */
  if(parallel){
    timer[0] = -MPI_Wtime();
    plan_forw = fftw_mpi_plan_dft_3d(
        n[0], n[1], n[2], in, out, MPI_COMM_WORLD, FFTW_FORWARD, tr_out | fftw_opt_flags);
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
        n[0], n[1], n[2], out, in, MPI_COMM_WORLD, FFTW_BACKWARD, tr_in | fftw_opt_flags);
    timer[1] += MPI_Wtime();
  } else {
    timer[1] = -MPI_Wtime();
    plan_back = fftw_plan_dft_3d(
        n[0], n[1], n[2], out, in, FFTW_BACKWARD, fftw_opt_flags);
    timer[1] += MPI_Wtime();
  }

  /* Initialize input with random numbers */
  pfft_init_input_complex_3d(n, local_ni, local_i_start,
      in);

  if(verbose)
    pfft_apr_complex_3d(in, local_ni, local_i_start, "PFFT Input", MPI_COMM_WORLD);

  /* execute parallel forward FFT */
  timer[2] = -MPI_Wtime();
  for(int t=0; t<loops; t++)
    fftw_execute(plan_forw);
  timer[2] += MPI_Wtime();

  if(verbose)
    pfft_apr_complex_3d(out, local_no, local_o_start, "PFFT Output", MPI_COMM_WORLD);
  
  /* execute parallel backward FFT */
  timer[3] = -MPI_Wtime();
  for(int t=0; t<loops; t++)
    fftw_execute(plan_back);
  timer[3] += MPI_Wtime();
  
  /* Scale data */
  for(ptrdiff_t l=0; l < local_ni[0] * local_ni[1] * local_ni[2]; l++)
    in[l] /= (n[0]*n[1]*n[2]);

  if(verbose)
    pfft_apr_complex_3d(in, local_ni, local_i_start, "Inputs after forward and backward PFFT", MPI_COMM_WORLD);
  
  /* Print error of back transformed data */
  MPI_Barrier(MPI_COMM_WORLD);
  err = pfft_check_output_complex_3d(n, local_ni, local_i_start, in, MPI_COMM_WORLD);

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

  pfft_printf(MPI_COMM_WORLD, "tune_forw = %6.2e; tune_back = %6.2e, exec_forw/loops = %6.2e, exec_back/loops = %6.2e\n", timer[0], timer[1], timer[2]/loops, timer[3]/loops);
  if(loops==1) pfft_printf(MPI_COMM_WORLD, "error = %6.2e\n", err);

  /* free mem and finalize */
  fftw_destroy_plan(plan_forw);
  fftw_destroy_plan(plan_back);
  fftw_free(in); if(!inplace) fftw_free(out);
}

static void init_parameters(
    int argc, char **argv,
    ptrdiff_t *n, int *np,
    unsigned *pfft_flags,
    int *loops, int *transposed, int *verbose, int *inplace,
    int *cmp_fftw, int *cmp_decomp, int *cmp_flags,
    int *print_timer
    )
{
  pfft_get_args(argc, argv, "-pfft_n", 3, PFFT_PTRDIFF_T, n);
  pfft_get_args(argc, argv, "-pfft_np", 3, PFFT_INT, np);
  pfft_get_args(argc, argv, "-pfft_loops", 1, PFFT_INT, loops);
  pfft_get_args(argc, argv, "-pfft_transposed", 0, PFFT_SWITCH, transposed);
  pfft_get_args(argc, argv, "-pfft_verbose", 0, PFFT_SWITCH, verbose);
  pfft_get_args(argc, argv, "-pfft_inplace", 0, PFFT_SWITCH, inplace);
  pfft_get_args(argc, argv, "-pfft_cmp_fftw", 0, PFFT_SWITCH, cmp_fftw);
  pfft_get_args(argc, argv, "-pfft_cmp_decomp", 0, PFFT_SWITCH, cmp_decomp);
  pfft_get_args(argc, argv, "-pfft_cmp_flags", 0, PFFT_SWITCH, cmp_flags);
  pfft_get_args(argc, argv, "-pfft_timer", 0, PFFT_SWITCH, print_timer);
  
  unsigned patience=0;
  pfft_get_args(argc, argv, "-pfft_patience", 1, PFFT_UNSIGNED, &patience);
  *pfft_flags = 0;
  switch(patience){
    case 1:  *pfft_flags |= PFFT_MEASURE;    break;
    case 2:  *pfft_flags |= PFFT_PATIENT;    break;
    case 3:  *pfft_flags |= PFFT_EXHAUSTIVE; break; 
    default: *pfft_flags |= PFFT_ESTIMATE; break;
  }
  
  int tune=0;
  pfft_get_args(argc, argv, "-pfft_tune", 0, PFFT_SWITCH, &tune);
  if(tune) *pfft_flags |= PFFT_TUNE;

  int destroy=0;
  pfft_get_args(argc, argv, "-pfft_destroy_input", 0, PFFT_SWITCH, &destroy);
  if(destroy) *pfft_flags |= PFFT_DESTROY_INPUT;

  pfft_printf(MPI_COMM_WORLD, "******************************************************************************************************\n");
  pfft_printf(MPI_COMM_WORLD, "* Computation of loops=%d parallel forward and backward FFTs (change with -pfft_loops *)\n", *loops);
#ifdef _OPENMP
  pfft_printf(MPI_COMM_WORLD, "* of with %d threads per process (change with OMP_NUM_THREADS=* before the command line) \n", pfft_get_nthreads());
#endif
  pfft_printf(MPI_COMM_WORLD, "* for n[0] x n[1] x n[2] = %td x %td x %td Fourier coefficients (change with -pfft_n * * *)\n", n[0], n[1], n[2]);
  pfft_printf(MPI_COMM_WORLD, "* on  np[0] x np[1] x np[2] = %td x %td x %td processes (change with -pfft_np * * *)\n", np[0], np[1], np[2]);
  
  pfft_printf(MPI_COMM_WORLD, "* with:\n");

  pfft_printf(MPI_COMM_WORLD, "*      - ");
  if(*transposed) pfft_printf(MPI_COMM_WORLD, "transposed data layout");
  else            pfft_printf(MPI_COMM_WORLD, "non-transposed data layout");
  pfft_printf(MPI_COMM_WORLD, " (change with -pfft_transposed)\n");
  
  pfft_printf(MPI_COMM_WORLD, "*      - ");
  if(*verbose) pfft_printf(MPI_COMM_WORLD, "verbose output");
  else         pfft_printf(MPI_COMM_WORLD, "non-verbose output");
  pfft_printf(MPI_COMM_WORLD, " (change with -pfft_verbose)\n");
  
  pfft_printf(MPI_COMM_WORLD, "*      - ");
  if(*inplace) pfft_printf(MPI_COMM_WORLD, "in-place transforms");
  else         pfft_printf(MPI_COMM_WORLD, "out-of-place transforms");
  pfft_printf(MPI_COMM_WORLD, " (change with -pfft_inplace)\n");
  
  pfft_printf(MPI_COMM_WORLD, "*      - ");
  if(*cmp_decomp) pfft_printf(MPI_COMM_WORLD, "enabled decomposition comparison");
  else            pfft_printf(MPI_COMM_WORLD, "disabled decomposition comparison");
  pfft_printf(MPI_COMM_WORLD, " (change with -pfft_cmp_decomp)\n");

  pfft_printf(MPI_COMM_WORLD, "*      - ");
  if(*cmp_fftw) pfft_printf(MPI_COMM_WORLD, "enabled FFTW comparison");
  else          pfft_printf(MPI_COMM_WORLD, "disabled FFTW comparison");
  pfft_printf(MPI_COMM_WORLD, " (change with -pfft_cmp_fftw)\n");

  pfft_printf(MPI_COMM_WORLD, "*      - ");
  if(*cmp_flags) pfft_printf(MPI_COMM_WORLD, "enabled comparison of all planner flags");
  else           pfft_printf(MPI_COMM_WORLD, "disabled comparison of all planner flags");
  pfft_printf(MPI_COMM_WORLD, " (change with -pfft_cmp_flags)\n");

  pfft_printf(MPI_COMM_WORLD, "*      - ");
  if(*cmp_flags) pfft_printf(MPI_COMM_WORLD, "enabled output of internal PFFT timer");
  else           pfft_printf(MPI_COMM_WORLD, "disabled output of internal PFFT timer");
  pfft_printf(MPI_COMM_WORLD, " (change with -pfft_timer)\n");

  if(!*cmp_flags) {
    pfft_printf(MPI_COMM_WORLD, "*      - ");
    pfft_printf(MPI_COMM_WORLD, "pfft_flags = ");
    if(*pfft_flags & PFFT_ESTIMATE)        pfft_printf(MPI_COMM_WORLD, "PFFT_ESTIMATE");
    else if(*pfft_flags & PFFT_PATIENT)    pfft_printf(MPI_COMM_WORLD, "PFFT_PATIENT");
    else if(*pfft_flags & PFFT_EXHAUSTIVE) pfft_printf(MPI_COMM_WORLD, "PFFT_EXHAUSTIVE");
    else                                   pfft_printf(MPI_COMM_WORLD, "PFFT_MEASURE");
    if(*pfft_flags & PFFT_TUNE)            pfft_printf(MPI_COMM_WORLD, " | PFFT_TUNE");
    else                                   pfft_printf(MPI_COMM_WORLD, " | PFFT_NO_TUNE");
    if(*pfft_flags & PFFT_DESTROY_INPUT)   pfft_printf(MPI_COMM_WORLD, " | PFFT_DESTROY_INPUT");
    pfft_printf(MPI_COMM_WORLD, "\n");
    pfft_printf(MPI_COMM_WORLD, "*      ");
    pfft_printf(MPI_COMM_WORLD, "  (change with [-pfft_patience  0|1|2|3] [-pfft_tune] [-pfft_destroy_input])\n");
  }

  pfft_printf(MPI_COMM_WORLD, "*******************************************************************************************************\n\n");

  if(*loops!=1) pfft_printf(MPI_COMM_WORLD, "!!! Warning: error checks are only available for loops=1 !!!\n");
  if(*inplace)  pfft_printf(MPI_COMM_WORLD, "!!! Warning: inplace transforms do not support DESTROY_INPUT flag !!!\n");
}

