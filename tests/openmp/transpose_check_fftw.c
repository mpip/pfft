/* This file is a test for the transpose algorithms for fftw with OpenMP */


/*#include <complex.h>*/
#include <fftw3.h>
#include <stdlib.h>
#include <sys/time.h>
/*#include <complex.h>*/
#include <pfft.h> 
#include <omp.h>

//#define MSIZE 8192
//#define NSIZE 8192


int initialize_matrix(fftw_complex*out,int m,int n)
{
  int i;int j;
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
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
  int rows=1,columns=1;
  int nthreads=1;
  int inplace=0;
  int measure=0;
  int runs=1;
  pfft_get_args(argc,argv,"-pfft_test_runs",1,PFFT_INT,&runs); 
  pfft_get_args(argc,argv,"-pfft_omp_threads",1,PFFT_INT,&nthreads); 
  pfft_get_args(argc,argv,"-pfft_omp_rows",1,PFFT_INT,&rows); 
  pfft_get_args(argc,argv,"-pfft_omp_columns",1,PFFT_INT,&columns); 
  pfft_get_args(argc,argv,"-pfft_omp_inplace",1,PFFT_SWITCH,&inplace); 
  pfft_get_args(argc,argv,"-pfft_omp_measure",1,PFFT_SWITCH,&measure); 
  if(inplace!=1) inplace=0;

/*  printf("Enter a number:");
  scanf("%d",&nthreads);
  printf("\n");*/
  fftw_complex*m1 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*rows*columns);
  fftw_complex*m2 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*rows*columns);
  fftw_iodim howmany_dims[2];
  howmany_dims[0].n=rows;
  howmany_dims[0].is=columns;
  howmany_dims[0].os=1;
  howmany_dims[1].n=columns;
  howmany_dims[1].is=1;
  howmany_dims[1].os=rows;
  const int howmany_rank = 2;
  // init OMP threads for fftw
  fftw_init_threads();
  fftw_plan_with_nthreads(nthreads);
  
  fftw_plan plan_transpose_outplace,plan_transpose_inplace,plansimple2d;

  int fftw_flags=FFTW_ESTIMATE;
  if(measure==1) fftw_flags=FFTW_MEASURE;


  plan_transpose_outplace = fftw_plan_guru_dft(/*rank*/ 0,NULL,howmany_rank, howmany_dims,m1,m2,FFTW_FORWARD,fftw_flags);
  plan_transpose_inplace = fftw_plan_guru_dft(/*rank*/ 0,NULL,howmany_rank, howmany_dims,m1,m1,FFTW_FORWARD,fftw_flags);
/*  plansimple2d=fftw_plan_2d()*/
  /*plansimple2d=fftw_plan_dft_2d(rows,columns,m1,m1,FFTW_FORWARD,FFTW_ESTIMATE);*/
  
  /*printf("# rows=%d, columns=%d, threads=%d,inplace=%d, measure=%d\n",rows,columns,nthreads,inplace,measure);
  initialize_matrix(m1,rows,columns);
  printf("# start calculation...\n");*/
  /*  printf("# check");*/
  /* printmatrix(m1,rows,columns);*/
  
  struct timeval t1,t2;
  gettimeofday(&t1,NULL);
   
  int i;
  for(i=0;i<runs;i++)
  {
    if(inplace==1) fftw_execute(plan_transpose_inplace);  
    if(inplace==0) fftw_execute(plan_transpose_outplace);  
  }
  /*if(selectplan==3) fftw_execute(plansimple2d);*/
  
  gettimeofday(&t2,NULL);
  float seconds=t2.tv_sec-t1.tv_sec+(t2.tv_usec-t1.tv_usec)/1000000.0;
  printf("# output: rows,cols,inplace,,measure,threads,runs,total time,average time\n");
  printf("%d %d %d %d %d %d %.6f %6f\n",rows,columns,inplace,measure,nthreads,runs,seconds,seconds/runs);
 /* printmatrix(m1,rows,columns);*/
  free(m1);
  free(m2);
  return 0;

}
