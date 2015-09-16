/* TODO: make this file about a self written matrix transpose algorithm */

/*#include <complex.h>*/
/*#include <fftw3.h>*/
#include <stdlib.h>
/*#include <time.h>*/
/*#include <complex.h>*/
/*#include <pfft.h>  */
#include <omp.h>
#include <stdio.h>

typedef double fftw_complex[2];


// out of place transposition for square matrix
int transpose_2(fftw_complex*in,fftw_complex*out,int rows,int columns)
{
  // rows and columns correspondent to in
  if(rows!=columns) return -1;
  if(in==out) return -1;
  int i,j;
  #pragma omp parallel for private(j)
  for(i=0;i<rows;i++)
  {
    for(j=0;j<rows;j++)
    {
      out[j*columns+i][0]=in[i*columns+j][0];
      out[j*columns+i][1]=in[i*columns+j][1];
    }
  }
  return 0;
}

// in_place transposition for square matrix
int transpose_1(fftw_complex*in,int rows, int columns)
{
  if(rows!=columns) return -1;
  int i,j;
  /*register*/ int x1,x2;
  double swap1,swap0;
// TODO: add parralelization
  #pragma omp parallel for private(j,x1,x2,swap0,swap1) 
  for(i=0;i<rows;i++)
  {
    for(j=0;j<i;j++)
    {
      x1=i*columns+j;
      x2=j*columns+i;
      swap0=in[x1][0];
      swap1=in[x1][1];
      
      in[x1][0]=in[x2][0];
      in[x1][1]=in[x2][1];
      in[x2][0]=swap0;
      in[x2][1]=swap1;
    }
  }
  return 0;
}

int transpose_3(fftw_complex*in, int rows,int columns)
{
  if(rows!=columns) return -1;
  int i,j;
  double*in2=(double*)in;
  /*register*/ int x1,x2;
  double swap1,swap0,swap2,swap3;
// TODO: add parralelization
  #pragma omp parallel for private(j,x1,x2,swap0,swap1) 
  for(i=0;i<rows;i++)
  {
    for(j=i+1;j<rows;j++)
    {
      x1=2*(i*columns+j);
      x2=2*(j*columns+i);
      swap0=in2[x1];
      swap1=in2[x1+1];     
      swap2=in2[x2];
      swap3=in2[x2+1];
      in2[x2]=swap0;
      in2[x2+1]=swap1;
      in2[x1]=swap2;
      in2[x1+1]=swap3;
    }
  }
  return 0;
}

int initialize_matrix(fftw_complex*out,int m,int n)
{
  int i;int j;
  for(i=0;i<m;i++)
    for(j=0;j<n;j++)
    {
      out[i*m+j][0]=2*i+j;
      out[i*m+j][1]=2*i-j;
    }
  return 0;
}

int printmatrix(fftw_complex*matrix,int m,int n)
{
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
  int rows=1,columns=1;
  int inplace=0;
  int nthreads=1;
  int runs=1;
  int print=0;
  if(argc==1)
  {
    printf("# usage: <program name> <rows> <columns> <inplace>\n");
    printf("# inplace should be 0 or 1 \n");
    return 0;
  }
  if(argc>1)
  {
    sscanf(argv[1],"%d",&rows);
  }
  if(argc>2)
  {
    sscanf(argv[2],"%d",&columns);
  }
  if(argc>3)
  {
    sscanf(argv[3],"%d",&inplace);
  }
  double wtime,avgtime;
  nthreads=omp_get_max_threads();
  printf("# Test program for transposing a Matrix with OpenMP\n");
  printf("# Number of avaliable Processors: %d\n",omp_get_num_procs());
/*  printf("# Number of Threads being used: %d\n",omp_get_max_threads());
  printf("Enter a number:");
  scanf("%d",&nthreads);
  printf("\n");*/
  fftw_complex*m1 = (fftw_complex *) malloc(sizeof(fftw_complex)*rows*columns);
  fftw_complex*m2 = (fftw_complex *) malloc(sizeof(fftw_complex)*rows*columns);
  initialize_matrix(m1,rows,columns);
  if(print==1) printmatrix(m1,rows,columns);
  if(inplace==1) 
  {
    int i=0;
    wtime=omp_get_wtime();
    for(i=0;i<runs;i++)
    {
      transpose_1(m1,rows,columns);
      /* this computation does not really make sense*/
    }
    wtime=omp_get_wtime()-wtime;
  }
  if(inplace==0) 
  {
    int i=0;
    wtime=omp_get_wtime();
    for(i=0;i<runs;i++)
    {
      transpose_2(m1,m2,rows,columns);
    }
    wtime=omp_get_wtime()-wtime;
  }
  if(inplace==3) 
  {
    int i=0;
    wtime=omp_get_wtime();
    for(i=0;i<runs;i++)
    {
      transpose_3(m1,rows,columns);
      /* this computation does not really make sense*/
    }
    wtime=omp_get_wtime()-wtime;
  }
  avgtime=wtime/runs;
  if(print==1) printmatrix(m1,rows,columns);
  printf("# Time for calculation used: %.6f\n",wtime);
  printf("# rows,columns,nthreads,inplace,runs,total time,avg time\n");
  printf("%d %d %d %d %d %.6f %.6f\n",rows,columns,nthreads,inplace,runs,wtime,avgtime);
 
  //  printmatrix(m1,rows,columns);
  /*struct timeval t1,t2;
  gettimeofday(&t1,NULL);
  gettimeofday(&t2,NULL);
  float seconds=t2.tv_sec-t1.tv_sec+(t2.tv_usec-t1.tv_usec)/1000000.0;
  printf("The measured time interval is %.6f seconds\n\n",seconds);*/
  free(m1);
  return 0;
}
