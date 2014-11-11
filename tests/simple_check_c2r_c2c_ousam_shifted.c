#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <pfft.h>

static pfft_complex semirandom(const ptrdiff_t *N, ptrdiff_t k0, ptrdiff_t k1, ptrdiff_t k2);
static void init_input(const ptrdiff_t *N, const ptrdiff_t *local_N, const ptrdiff_t *local_N_start, pfft_complex *data);
static double compare_c2c_c2r(const ptrdiff_t *local_N_c, const ptrdiff_t *local_N_r, const pfft_complex *data_c, const pfft_complex *data_r, MPI_Comm comm);

int main(int argc, char **argv)
{
  int np[2];
  ptrdiff_t n[3], ni[3], howmany;
  double err;

  ptrdiff_t alloc_local_c;
  ptrdiff_t local_ni_c[3], local_i_start_c[3];
  ptrdiff_t local_no_c[3], local_o_start_c[3];
  pfft_complex *in_c, *out_c;
  pfft_plan plan_forw_c=NULL, plan_back_c=NULL;

  ptrdiff_t alloc_local_r, alloc_local_forw, alloc_local_back;
  ptrdiff_t local_ni_r[3], local_i_start_r[3];
  ptrdiff_t local_no_r[3], local_o_start_r[3];
  pfft_complex *in_r;
  double *out_r;
  pfft_plan plan_forw_r=NULL, plan_back_r=NULL;

  MPI_Comm comm_cart_2d;

  /* Set size of FFT and process mesh */
  ni[0] = 4; ni[1] = 4; ni[2] = 4;
  n[0] = 6; n[1] = 6; n[2] = 6;
  np[0] = 2; np[1] = 2;
  howmany = 1;

  /* Initialize MPI and PFFT */
  MPI_Init(&argc, &argv);
  pfft_init();

  /* Create two-dimensional process grid of size np[0] x np[1], if possible */
  if( pfft_create_procmesh_2d(MPI_COMM_WORLD, np[0], np[1], &comm_cart_2d) ){
    pfft_fprintf(MPI_COMM_WORLD, stderr, "Error: This test file only works with %d processes.\n", np[0]*np[1]);
    MPI_Finalize();
    return 1;
  }

  /* Get parameters of data distribution */
  alloc_local_forw = pfft_local_size_many_dft(3, n, ni, n, howmany, PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, comm_cart_2d, PFFT_TRANSPOSED_NONE| PFFT_SHIFTED_IN | PFFT_SHIFTED_OUT,
      local_ni_c, local_i_start_c, local_no_c, local_o_start_c);
  alloc_local_back = pfft_local_size_many_dft(3, n, n, ni, howmany, PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, comm_cart_2d, PFFT_TRANSPOSED_NONE| PFFT_SHIFTED_IN | PFFT_SHIFTED_OUT,
      local_no_c, local_o_start_c, local_ni_c, local_i_start_c);

  alloc_local_c = (alloc_local_forw > alloc_local_back) ? alloc_local_forw : alloc_local_back;

  alloc_local_forw = pfft_local_size_many_dft_c2r(3, n, ni, n, howmany, PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, comm_cart_2d, PFFT_TRANSPOSED_NONE| PFFT_SHIFTED_IN | PFFT_SHIFTED_OUT,
      local_ni_r, local_i_start_r, local_no_r, local_o_start_r);
  alloc_local_back = pfft_local_size_many_dft_r2c(3, n, n, ni, howmany, PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, comm_cart_2d, PFFT_TRANSPOSED_NONE| PFFT_SHIFTED_IN | PFFT_SHIFTED_OUT,
      local_no_r, local_o_start_r, local_ni_r, local_i_start_r);

  alloc_local_r = (alloc_local_forw > alloc_local_back) ? alloc_local_forw : alloc_local_back;


  /* Allocate memory */
  in_c  = pfft_alloc_complex(alloc_local_c);
  out_c = pfft_alloc_complex(alloc_local_c);
  in_r  = pfft_alloc_complex(alloc_local_r);
  out_r = pfft_alloc_real(2*alloc_local_r);


  /* Plan parallel forward FFT */
  plan_forw_c = pfft_plan_many_dft(3,
      n, ni, n, howmany, PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, in_c, out_c, comm_cart_2d, PFFT_FORWARD, PFFT_TRANSPOSED_NONE| PFFT_MEASURE| PFFT_DESTROY_INPUT| PFFT_SHIFTED_IN | PFFT_SHIFTED_OUT);
  plan_forw_r = pfft_plan_many_dft_c2r(3,
      n, ni, n, howmany, PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, in_r, out_r, comm_cart_2d, PFFT_FORWARD, PFFT_TRANSPOSED_NONE| PFFT_MEASURE| PFFT_DESTROY_INPUT| PFFT_SHIFTED_IN | PFFT_SHIFTED_OUT);

  /* Plan parallel backward FFT */
  plan_back_c = pfft_plan_many_dft(3,
      n, n, ni, howmany, PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, out_c, in_c, comm_cart_2d, PFFT_BACKWARD, PFFT_TRANSPOSED_NONE| PFFT_MEASURE| PFFT_DESTROY_INPUT| PFFT_SHIFTED_IN | PFFT_SHIFTED_OUT);
  plan_back_r = pfft_plan_many_dft_r2c(3,
      n, n, ni, howmany, PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, out_r, in_r, comm_cart_2d, PFFT_BACKWARD, PFFT_TRANSPOSED_NONE| PFFT_MEASURE| PFFT_DESTROY_INPUT| PFFT_SHIFTED_IN | PFFT_SHIFTED_OUT);

  /* Initialize input with random numbers */
  init_input(ni, local_ni_c, local_i_start_c, in_c);
  init_input(ni, local_ni_r, local_i_start_r, in_r);

//   pfft_apr_complex_3d(in_c, local_ni_c, local_i_start_c, "c2c input:\n", comm_cart_2d);
//   pfft_apr_complex_3d(in_r, local_ni_r, local_i_start_r, "c2r input:\n", comm_cart_2d);

  /* execute parallel forward FFT */
  pfft_execute(plan_forw_c);

  /* clear the old input */

  pfft_execute(plan_forw_r);

//   pfft_apr_complex_3d(out_c, local_no_c, local_o_start_c, "c2c output:\n", comm_cart_2d);
//   pfft_apr_real_3d(out_r, local_no_r, local_o_start_r, "c2r output:\n", comm_cart_2d);

  /* execute parallel backward FFT */
  pfft_execute(plan_back_c);
  pfft_execute(plan_back_r);

//   pfft_apr_complex_3d(in_c, local_ni_c, local_i_start_c, "c2c^ output:\n", comm_cart_2d);
//   pfft_apr_complex_3d(in_r, local_ni_r, local_i_start_r, "c2r^ output:\n", comm_cart_2d);

  /* Scale data */
  for(ptrdiff_t l=0; l < local_ni_c[0] * local_ni_c[1] * local_ni_c[2]; l++)
    in_c[l] /= (n[0]*n[1]*n[2]);
  for(ptrdiff_t l=0; l < local_ni_r[0] * local_ni_r[1] * local_ni_r[2]; l++)
    in_r[l] /= (n[0]*n[1]*n[2]);

  /* Print error of back transformed data */
  err = compare_c2c_c2r(local_ni_c, local_ni_r, in_c, in_r, comm_cart_2d);

  pfft_printf(comm_cart_2d, "Error after one forward and backward trafo of size n=(%td, %td, %td):\n", n[0], n[1], n[2]);
  pfft_printf(comm_cart_2d, "maxerror = %6.2e;\n", err);
  
  /* free mem and finalize */
  pfft_destroy_plan(plan_forw_c);
  pfft_destroy_plan(plan_back_c);
  pfft_destroy_plan(plan_forw_r);
  pfft_destroy_plan(plan_back_r);
  MPI_Comm_free(&comm_cart_2d);
  pfft_free(in_c); pfft_free(out_c);
  pfft_free(in_r); pfft_free(out_r);
  MPI_Finalize();
  return 0;
}

#define DATA_INIT(i) (( (double)1000 ) / ( (double)( (i) == 0 ? 1 : i) ))

static pfft_complex semirandom(const ptrdiff_t *N, ptrdiff_t k0, ptrdiff_t k1, ptrdiff_t k2) {
  if ((k0 == -N[0]/2) || (k1 == -N[1]/2) || (k2 == -N[2]/2))
    return 0 + 0*I;
  if (k0 == (N[0]+1)/2)
    return semirandom(N, -k0, k1, k2);
  if (k1 == (N[1]+1)/2)
    return semirandom(N, k0, -k1, k2);
  if (k2 == (N[2]+1)/2)
    return semirandom(N, k0, k1, -k2);

  ptrdiff_t l0 = (k0) ? k0+N[0] : 1;
  ptrdiff_t l1 = (k1) ? k1-2*N[1] : 1;
  ptrdiff_t l2 = (k2) ? k2+N[2] : 1;

  double re = DATA_INIT(l0+l1*l2) + 3*DATA_INIT(l0*l1+l2) + DATA_INIT(l2*l0+l1);
  double im = 3*DATA_INIT(l0+l1*l1) + DATA_INIT(l1+l2*l2) - DATA_INIT(l2+l0*l0);

  return re + im*I;
}

static void init_input(const ptrdiff_t *N, const ptrdiff_t *local_N, const ptrdiff_t *local_N_start, pfft_complex *data)
{
  int m = 0;
  for(ptrdiff_t k0=local_N_start[0]; k0<local_N_start[0]+local_N[0]; k0++)
    for(ptrdiff_t k1=local_N_start[1]; k1<local_N_start[1]+local_N[1]; k1++)
      for(ptrdiff_t k2=local_N_start[2]; k2<local_N_start[2]+local_N[2]; k2++, m++)
        data[m] = semirandom(N, k0, k1, k2) + conj(semirandom(N, -k0, -k1, -k2));
}

static double compare_c2c_c2r(const ptrdiff_t *local_N_c, const ptrdiff_t *local_N_r, const pfft_complex *data_c, const pfft_complex *data_r, MPI_Comm comm)
{
  double err = 0;
  double glob_max_err = 0;

  for (int k0=0; k0 < local_N_r[0]; k0++)
    for (int k1=1; k1 < local_N_r[1]; k1++)
      for (int k2=2; k2 < local_N_r[2]; k2++) {
        double complex r = data_r[k2 + k1*local_N_r[2] + k0*local_N_r[1]*local_N_r[2]];
        double complex c = data_c[k2 + k1*local_N_c[2] + k0*local_N_c[1]*local_N_c[2]];
        double re = creal(r) - creal(c);
        double im = cimag(r) - cimag(c);
        double tmp = sqrt(re*re + im*im);
        if (tmp > err)
          err = tmp;
      }
  MPI_Allreduce(&err, &glob_max_err, 1, MPI_DOUBLE, MPI_MAX, comm);
  return glob_max_err;
}

