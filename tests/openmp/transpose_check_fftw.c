/* This file is a test for the transpose algorithms for fftw with OpenMP */

/* however, most of the code is just crap right now */


/*#include <complex.h>*/
#include <fftw3.h>
/*#include <complex.h>*/
/*#include <pfft.h> */
#include <omp.h>




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
  return 0;
}


int main(int argc, char **argv)
{
  int rows=4,columns=4;
  fftw_complex*m1 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*rows*columns);
  initialize_matrix(m1,rows,columns);

  /*
  TODO: TASK1: use the fftw transpose algorithms
  */
  printmatrix(m1,rows,columns);

  return 0;


}
