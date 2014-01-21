#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <pfft.h>

void init_input_c2r(const int rnk_n, const ptrdiff_t *n, const ptrdiff_t *local_ni, const ptrdiff_t *local_i_start, pfft_complex *data);

int main(int argc, char **argv)
{
  int np[2];
  ptrdiff_t n[3], ni[3], no[3];
  ptrdiff_t alloc_local_forw, alloc_local_back, alloc_local, howmany;
  ptrdiff_t local_ni[3], local_i_start[3];
  ptrdiff_t local_n[3], local_start[3];
  ptrdiff_t local_no[3], local_o_start[3];
  double err, *out;
  pfft_complex *in, *in_check;
  pfft_plan plan_forw=NULL, plan_back=NULL;
  MPI_Comm comm_cart_2d;
  
  /* Set size of FFT and process mesh */
  ni[0] = ni[1] = ni[2] = no[0] = no[1] = no[2] = 16;
  n[0] = 29; n[1] = 27; n[2] = 31;
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
  alloc_local_forw = pfft_local_size_many_dft_c2r(3, n, ni, n, howmany,
      PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS,
      comm_cart_2d, PFFT_TRANSPOSED_NONE,
      local_ni, local_i_start, local_n, local_start);

  alloc_local_back = pfft_local_size_many_dft_r2c(3, n, n, no, howmany,
      PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS,
      comm_cart_2d, PFFT_TRANSPOSED_NONE,
      local_n, local_start, local_no, local_o_start);
  
  alloc_local = (alloc_local_forw > alloc_local_back) ?
    alloc_local_forw : alloc_local_back;

  /* Allocate memory */
  in  = pfft_alloc_complex(alloc_local);
  in_check = pfft_alloc_complex(alloc_local);
  out = pfft_alloc_real(2 * alloc_local);

  /* Plan parallel forward FFT */
  plan_forw = pfft_plan_many_dft_c2r(
      3, n, ni, n, howmany, PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS,
      in, out, comm_cart_2d, PFFT_FORWARD, PFFT_TRANSPOSED_NONE| PFFT_MEASURE| PFFT_DESTROY_INPUT);
  
  /* Plan parallel backward FFT */
  plan_back = pfft_plan_many_dft_r2c(
      3, n, n, no, howmany, PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS,
      out, in, comm_cart_2d, PFFT_BACKWARD, PFFT_TRANSPOSED_NONE| PFFT_MEASURE| PFFT_DESTROY_INPUT);

  /* Initialize input with random numbers */
//  m = 0;
//  for(int k0=0; k0<n[0]; k0++){
//    for(int k1=0; k1<n[1]; k1++){
//      for(int k2=0; k2<n[2]; k2++, m++)
//        in[m] = 1000/(m+1);
//      for(int k2=n[2]; k2<n[2]/2+1; k2++, m++)
//        in[m] = 0;
//    }
//  }
  
  init_input_c2r(3, ni, local_ni, local_i_start, in_check);
  init_input_c2r(3, ni, local_ni, local_i_start, in);

  /* execute parallel forward FFT */
  pfft_execute(plan_forw);
  
  /* execute parallel backward FFT */
  pfft_execute(plan_back);
  
  /* Scale data */
  for(ptrdiff_t l=0; l < local_ni[0] * local_ni[1] * local_ni[2]; l++)
    in[l] /= (n[0]*n[1]*n[2]);
  
  /* Print error of back transformed data */
  MPI_Barrier(MPI_COMM_WORLD);

  err = 0;
  for(ptrdiff_t l=0; l < local_ni[0] * local_ni[1] * local_ni[2]; l++) {
    double re = creal(in[l]) - creal(in_check[l]);
    double im = cimag(in[l]) - cimag(in_check[l]);
    double tmp = sqrt(re*re + im*im);
    if (tmp > err)
      err = tmp;
  }
    
//  err = pfft_check_output_c2c(3, pn, local_ni, local_i_start, in, comm_cart_2d);
  pfft_printf(comm_cart_2d, "Error after one forward and backward trafo of size n=(%td, %td, %td):\n", n[0], n[1], n[2]); 
  pfft_printf(comm_cart_2d, "maxerror = %6.2e;\n", err);
  
  /* free mem and finalize */
  pfft_destroy_plan(plan_forw);
  pfft_destroy_plan(plan_back);
  MPI_Comm_free(&comm_cart_2d);
  pfft_free(in); pfft_free(out);
  MPI_Finalize();
  return 0;
}

void init_input_c2r(const int rnk_n, const ptrdiff_t *n, const ptrdiff_t *local_ni, const ptrdiff_t *local_i_start, pfft_complex *data)
{
  ptrdiff_t ln_tot = pfft_prod_INT(rnk_n, local_ni);
  ptrdiff_t *pn = malloc(sizeof(ptrdiff_t) * rnk_n);

  for (int k=0; k<rnk_n; k++)
    pn[k] = n[k];
  pn[rnk_n-1] = pn[rnk_n-1]/2 + 1;

  pfft_init_input_complex(rnk_n, pn, local_ni, local_i_start, data);

  if (local_i_start[rnk_n-1] == 0)
    for(int k=0; k<ln_tot; k += local_ni[rnk_n-1])
      data[k] = 0;
  if ((n[rnk_n-1]%2 == 0) && (local_i_start[rnk_n-1] + local_ni[rnk_n-1] == pn[rnk_n-1]))
    for(int k=pn[rnk_n-1]-1 - local_i_start[rnk_n-1]; k<ln_tot; k+= local_ni[rnk_n-1])
      data[k] = creal((double complex)data[k]);
}
