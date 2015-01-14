/* This file is a test for the transpose algorithms for fftw with OpenMP */

/* TODO/ FIXME: fix segmentation fault somewhere */


/*#include <complex.h>*/
#include <fftw3.h>
/*#include <complex.h>*/
/*#include <pfft.h> */
#include <omp.h>

#define MSIZE 4
#define NSIZE 4


int initialize_matrix(fftw_complex*out,int m,int n)
{
  int i;int j;
  for(i=0;i<m;i++)
    for(j=0;j<n;j++)
    {
      out[i*m+j][0]=2*i+j;
    }
  return 0;
}


int printmatrix(fftw_complex*matrix,int m,int n)
{
  /*only print real part*/
  int i;int j;
  for(i=0;i<m;i++)
  {  
    for(j=0;j<n;j++)
    {
        printf(" %7.3f",matrix[i*m+j][0]);
    }
    printf("\n");
  }  
  printf("\nend of print\n");
  return 0;
}


int main(int argc, char **argv)
{
  int rows=MSIZE,columns=NSIZE;
  fftw_complex*m1 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*rows*columns);
  fftw_iodim howmany_dims[2];
  howmany_dims[0].n=rows;
  howmany_dims[0].is=columns;
  howmany_dims[0].os=1;
  howmany_dims[1].n=columns;
  howmany_dims[1].is=1;
  howmany_dims[1].is=rows;
  const int howmany_rank = sizeof(howmany_dims)/sizeof(howmany_dims[0]);

  fftw_plan plan_transpose = fftw_plan_guru_dft(/*rank*/ 0,NULL,howmany_rank, howmany_dims,m1,m1,FFTW_FORWARD,FFTW_ESTIMATE); /* this causes the error most likely */


  initialize_matrix(m1,rows,columns);
  printmatrix(m1,rows,columns);
  printf("... segfault after this line:\n");
  
  fftw_execute(plan_transpose);  /* FIXME: this line leads to a segfault. possibly wrong call of fftw_plan_guru_dft */
  /*
  TODO: TASK1: use the fftw transpose algorithms
  */
  printmatrix(m1,rows,columns);

  return 0;


}
