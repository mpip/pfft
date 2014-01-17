#include <complex.h>
#include <pfft.h>

static void init_parameters(
    int argc, char **argv,
    int *n, int *iter,
    int *inplace, int* patience
    )
{
  pfft_get_args(argc, argv, "-pfft_n", 3, PFFT_INT, n);
  pfft_get_args(argc, argv, "-pfft_iter", 1, PFFT_INT, iter); 
  pfft_get_args(argc, argv, "-pfft_ip", 1, PFFT_INT, inplace);
  pfft_get_args(argc, argv, "-pfft_pat", 1, PFFT_INT, patience);
}


int main(int argc, char **argv)
{
  int n[3];
  pfft_complex *in, *out;
  FFTW(plan) plan_forw=NULL, plan_back=NULL;
  double err, time, time_fftw[2], max_time_fftw[2];
  unsigned fftw_flag;

  /* setup default parameters */
  int iter = 10, inplace = 0, patience = 0;  
  
  /* Set size of FFT and process mesh */
  n[0] = n[1] = n[2] = 16;
 
  /* Initialize MPI and PFFT */
  MPI_Init(&argc, &argv);
  pfft_init();

  /* read parameters from command line */
  init_parameters(argc, argv, n, &iter, &inplace, &patience);

  /* setup FFTWs planing depth */  
  switch(patience){
    case 1: fftw_flag = FFTW_MEASURE; break;
    case 2: fftw_flag = FFTW_PATIENT; break;
    case 3: fftw_flag = FFTW_EXHAUSTIVE; break;
    default: fftw_flag = FFTW_ESTIMATE;
  }
  
  if(!inplace)
    fftw_flag |= FFTW_DESTROY_INPUT;

  /* Allocate memory */
  in = pfft_alloc_complex(n[0]*n[1]*n[2]);
  out = (inplace) ? in : pfft_alloc_complex(n[0]*n[1]*n[2]);

  /* We often want to scale large FFTs, which do not fit on few processes. */
  if( (in == NULL) || (out == NULL)){
    fprintf(stderr, "!!! Error: Not enough memory to allocate input/output arrays !!!\n");
    MPI_Finalize();
    MPI_Finalize();
    return 1;
  }

  ptrdiff_t local_ni[3], local_i_start[3], n_ptr[3];
  for(int t=0; t<3; t++){
    local_i_start[t] = 0;
    n_ptr[t] = local_ni[t] = (ptrdiff_t) n[t];
  }
  
  plan_forw = FFTW(plan_dft_3d)(n[0], n[1], n[2], in, out, FFTW_FORWARD, fftw_flag);
  plan_back = FFTW(plan_dft_3d)(n[0], n[1], n[2], out, in, FFTW_BACKWARD, fftw_flag);
  
  /* Initialize input with random numbers */
  pfft_init_input_complex_3d(n_ptr, local_ni, local_i_start,
      in);

  time_fftw[0] = time_fftw[1] = 0;
  for(int t=0; t<iter; t++){
    /* execute parallel forward FFT */
    time_fftw[0] -= MPI_Wtime();
    FFTW(execute)(plan_forw);
    time_fftw[0] += MPI_Wtime();
  
    /* execute parallel backward FFT */
    time_fftw[1] -= MPI_Wtime();
    FFTW(execute)(plan_back);
    time_fftw[1] += MPI_Wtime();
  }
  
  /* Scale data */
  for(int t=0; t<iter; t++)
    for(ptrdiff_t l=0; l < n[0] * n[1] * n[2]; l++)
      in[l] /= (n[0]*n[1]*n[2]);

  printf("fftw_forw = %.2e, fftw_back = %.2e\n", time_fftw[0]/iter, time_fftw[1]/iter);
 
  err = pfft_check_output_complex_3d(n_ptr, local_ni, local_i_start, in, MPI_COMM_WORLD);
  printf("Error after several forward and backward FFTWs of size n=(%td, %td, %td):\n", n[0], n[1], n[2]); 
  printf("maxerror = %6.2e;\n", err);
  
  /* free mem and finalize */
  FFTW(destroy_plan)(plan_forw);
  FFTW(destroy_plan)(plan_back);
  pfft_free(in); if(!inplace) pfft_free(out);
 
  MPI_Finalize();
  return 0;
}
