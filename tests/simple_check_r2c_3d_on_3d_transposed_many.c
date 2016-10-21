#include <complex.h>
#include <pfft.h>
#include <math.h>

int main(int argc, char **argv)
{
  int np[3];
  ptrdiff_t n[3];
  ptrdiff_t alloc_local;
  ptrdiff_t local_ni[3], local_i_start[3];
  ptrdiff_t local_no[3], local_o_start[3];
  ptrdiff_t howmany;

  double err, *in, *check;
  double dx,dy,dz, normfac;
  pfft_complex *out;
  pfft_plan plan_forw=NULL, plan_back=NULL;
  MPI_Comm comm_cart_3d;
  int rank;

  howmany = 2;

  /* Set size of FFT and process mesh */
  n[0] = 8; n[1] = 8; n[2] = 8;

  np[0] = 1; np[1] = 1; np[2] = 2; // IF I CHOOSE np[2] > 1, COMBINED WITH howmany = 2, THEN THE TEST FAILS

  dx = 1.0/n[0];
  dy = 1.0/n[1];
  dz = 1.0/n[2];

  /* Initialize MPI and PFFT */
  MPI_Init(&argc, &argv);
  pfft_init();

  /* Create three-dimensional process grid of size np[0] x np[1] x np[2], if possible */
  if( pfft_create_procmesh(3, MPI_COMM_WORLD, np, &comm_cart_3d) ){
    pfft_fprintf(MPI_COMM_WORLD, stderr, "Error: This test file only works with %d processes.\n", np[0]*np[1]*np[2]);
    MPI_Finalize();
    return 1;
  }
  MPI_Comm_rank (comm_cart_3d, &rank);

  /* Get parameters of data distribution */
  alloc_local = pfft_local_size_many_dft_r2c(3, n, n, n, howmany, PFFT_DEFAULT_BLOCKS,PFFT_DEFAULT_BLOCKS, comm_cart_3d, PFFT_TRANSPOSED_OUT,
      local_ni, local_i_start, local_no, local_o_start);

  printf("rnk %d : in start %td %td %td, size %td %td %td, total %td\n",rank,local_i_start[0],local_i_start[1],local_i_start[2],local_ni[0],local_ni[1],local_ni[2],2*alloc_local);
  printf("\t\t\t\t\t\t\trnk %d : out start %td %td %td, size %td %td %td, total %td\n",rank,local_o_start[0],local_o_start[1],local_o_start[2],local_no[0],local_no[1],local_no[2],alloc_local);

  // allocate memory
  check = pfft_alloc_real(local_ni[0]*local_ni[1]*local_ni[2]*howmany);

  //  /* Allocate memory */
  in  = pfft_alloc_real(2 * alloc_local);
  // in has twice as much as we need but this is required for in-place transform. 
  out = pfft_alloc_complex(alloc_local);
  // out is allocated like this, should be the same as out = pfft_alloc_complex(local_no[0]*local_no[1]*local_no[2]);

  /* Plan parallel forward FFT */
  plan_forw = pfft_plan_many_dft_r2c(3, n, n, n, howmany,PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, in, out, comm_cart_3d, PFFT_FORWARD, PFFT_TRANSPOSED_OUT | PFFT_MEASURE| PFFT_DESTROY_INPUT);

  /* Plan parallel backward FFT */
  plan_back = pfft_plan_many_dft_c2r(3, n, n, n, howmany,PFFT_DEFAULT_BLOCKS, PFFT_DEFAULT_BLOCKS, out, in, comm_cart_3d, PFFT_BACKWARD, PFFT_TRANSPOSED_IN | PFFT_MEASURE| PFFT_DESTROY_INPUT);

  // initialize input with sine wave
  for(ptrdiff_t i=0; i<local_ni[0]; i++)
    for(ptrdiff_t j=0; j<local_ni[1]; j++)
      for(ptrdiff_t k=0; k<local_ni[2]; k++)
      {
        ptrdiff_t idx = (i*local_ni[1] + j)*local_ni[2] + k;
        idx = idx*howmany; // correct linear idx for the components

        double myx = (i + local_i_start[0])*dx;
        double myy = (j + local_i_start[1])*dy;
        double myz = (k + local_i_start[2])*dz;

        in[idx] = -3.0*(2.0*M_PI)*(2.0*M_PI)*(sin(2.0*M_PI*myx) * sin(2.0*M_PI*myy) * sin(2.0*M_PI*myz));
        check[idx] = sin(2.0*M_PI*myx) * sin(2.0*M_PI*myy) * sin(2.0*M_PI*myz);

        // second component has a different wave (with 100 x scale factor)
        if(howmany==2)
        {
          in[idx+1] =  -3.0*(2.0*M_PI)*(2.0*M_PI)*100*(cos(2.0*M_PI*myx) * cos(2.0*M_PI*myy) * cos(2.0*M_PI*myz));
          check[idx+1] =100*(cos(2.0*M_PI*myx) * cos(2.0*M_PI*myy) * cos(2.0*M_PI*myz));            
        }
      }

  /* execute parallel forward FFT */
  pfft_execute(plan_forw);

  normfac = 1.0/n[0];
  normfac /= n[1];
  normfac /= n[2];

  for(ptrdiff_t j=0; j<local_no[1]; j++)
    for(ptrdiff_t k=0; k<local_no[2]; k++)
      for(ptrdiff_t i=0; i<local_no[0]; i++)
      {

        // get the wave number (we have even n[i])
        ptrdiff_t kx,ky,kz;
        kx = (i + local_o_start[0] <= n[0]/2) ? i + local_o_start[0] : -(n[0] - i - local_o_start[0]);
        ky = (j + local_o_start[1] <= n[1]/2) ? j + local_o_start[1] : -(n[1] - j - local_o_start[1]);
        kz = (k + local_o_start[2] <= n[2]/2) ? k + local_o_start[2] : -(n[2] - k - local_o_start[2]);

        double rkx = kx*2.0*M_PI; 
        double rky = ky*2.0*M_PI;
        double rkz = kz*2.0*M_PI;

        // invert the laplacian
        double inv_denom;
        if(kx==0 && ky==0 && kz==0)
        {
          inv_denom = 0.0;
        }
        else
        {
          inv_denom = -1.0 / (rkx*rkx + rky*rky + rkz*rkz);
        }

        // rescale the output
        inv_denom *= normfac;

        ptrdiff_t idx = (j*local_no[2] + k)*local_no[0] + i;
        idx = idx*howmany; // correct linear idx for the components

        for(ptrdiff_t c=0; c<howmany; ++c)
          out[idx+c] = inv_denom*creal(out[idx+c]) + 1i*inv_denom*cimag(out[idx+c]);
      }


  /* execute parallel backward FFT */
  pfft_execute(plan_back);

  /* Print error of back transformed data */

  MPI_Barrier(MPI_COMM_WORLD);

  for(ptrdiff_t c=0; c<howmany; ++c)
  {
    err = 0.0;


    pfft_printf(comm_cart_3d, "Error after one forward and backward trafo of size n=(%td, %td, %td), component %d:\n", n[0], n[1], n[2], c);
    pfft_printf(comm_cart_3d, "maxerror = %6.2e;\n", err);
  }

  /* free mem and finalize */
  pfft_destroy_plan(plan_forw);
  pfft_destroy_plan(plan_back);
  MPI_Comm_free(&comm_cart_3d);
  pfft_free(in); pfft_free(out); pfft_free(check);

  MPI_Finalize();
  return 0;
}
