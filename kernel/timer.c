/*
 * Copyright (c) 2011-2013 Michael Pippig
 *
 * This file is part of PFFT.
 *
 * PFFT is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PFFT is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PFFT.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "pfft.h"
#include "ipfft.h"
#include <string.h>     /* strcmp */


static void reset_timer(
    PX(timer) ths);
static size_t length(
    const PX(timer) ths);
static void fprint_average_timer(
    MPI_Comm comm, FILE *file, const PX(plan) ths, unsigned flags);
static void write_info_header(
    MPI_Comm comm, FILE *file);
static void write_run_specific_infos(
    MPI_Comm comm, FILE *file, const PX(plan) ths);
static int get_matlab_index(
    MPI_Comm comm);

static void fprint_average_timer_prefixed(
    MPI_Comm comm, FILE *file, const char *prefix,
    const PX(timer) ths, unsigned flags);

static int file_exists(
    const char *name);
static FILE* open_or_create_file_to_append(
    MPI_Comm comm, const char *name);


/*********************************************
 * Functions to operate on timer in PX(plan)
 ********************************************/

PX(timer) PX(get_timer)(
    const PX(plan) ths
    )
{
  return PX(copy_timer)(ths->timer);
}


void PX(reset_timer)(
    PX(plan) ths
    )
{
  reset_timer(ths->timer);
}

void PX(print_average_timer)(
    const PX(plan) ths, MPI_Comm comm
    )
{
  write_info_header(comm, stdout);
  write_run_specific_infos(comm, stdout, ths);
  fprint_average_timer(comm, stdout, ths, PFFTI_PRINT_TIMER_BASIC);
}

void PX(print_average_timer_adv)(
    const PX(plan) ths, MPI_Comm comm
    )
{
  PX(print_average_timer)(ths, comm);
  fprint_average_timer(comm, stdout, ths, PFFTI_PRINT_TIMER_ADV);
}

void PX(write_average_timer)(
    const PX(plan) ths, const char *name, MPI_Comm comm
    )
{
  int newfile;
  FILE *f;
  
  newfile = !file_exists(name);
  f = open_or_create_file_to_append(comm, name);
  if( newfile )
    write_info_header(comm, f);
  
  write_run_specific_infos(comm, f, ths);
  fprint_average_timer(comm, f, ths, PFFTI_PRINT_TIMER_BASIC);
  
  fclose(f);
}

void PX(write_average_timer_adv)(
    const PX(plan) ths, const char *name, MPI_Comm comm
    )
{
  FILE *f;
  PX(write_average_timer)(ths, name, comm);
  f = open_or_create_file_to_append(comm, name);
  fprint_average_timer(comm, f, ths, PFFTI_PRINT_TIMER_ADV);
  fclose(f);
}

static void fprint_average_timer(
    MPI_Comm comm, FILE *file, const PX(plan) ths, unsigned flags
    )
{
  if(ths->sign == FFTW_FORWARD)
    fprint_average_timer_prefixed(comm, file, "pfft_forw", ths->timer, flags);
  else
    fprint_average_timer_prefixed(comm, file, "pfft_back", ths->timer, flags);
}

static void write_info_header(
    MPI_Comm comm, FILE *file
    )
{
  PX(fprintf)(comm, file, "\n%% n - FFT size\n");
  PX(fprintf)(comm, file, "%% np - process grid\n");
  PX(fprintf)(comm, file, "%% procs - number of processes\n");
  PX(fprintf)(comm, file, "%% pfft - PFFT runtime\n");
  PX(fprintf)(comm, file, "%% index(i) = log(procs(i)) + 1\n");
}

static void write_run_specific_infos(
    MPI_Comm comm, FILE *file, const PX(plan) ths
    )
{
  int size, idx = get_matlab_index(comm);
  MPI_Comm_size(comm, &size);
  
  if(ths->pfft_flags & PFFT_ESTIMATE)
    PX(fprintf)(comm, file, "\n%% pfft_flags == PFFT_ESTIMATE");
  else if(ths->pfft_flags & PFFT_PATIENT)
    PX(fprintf)(comm, file, "\n%% pfft_flags == PFFT_PATIENT");
  else if(ths->pfft_flags & PFFT_EXHAUSTIVE)
    PX(fprintf)(comm, file, "\n%% pfft_flags == PFFT_EXHAUSTIVE");
  else
    PX(fprintf)(comm, file, "\n%% pfft_flags == PFFT_MEASURE");
  if(ths->pfft_flags & PFFT_TRANSPOSED_IN)
    PX(fprintf)(comm, file, " | PFFT_TRANSPOSED_IN");
  if(ths->pfft_flags & PFFT_TRANSPOSED_OUT)
    PX(fprintf)(comm, file, " | PFFT_TRANSPOSED_OUT");
  if(ths->pfft_flags & PFFT_SHIFTED_IN)
    PX(fprintf)(comm, file, " | PFFT_SHIFTED_IN");
  if(ths->pfft_flags & PFFT_SHIFTED_OUT)
    PX(fprintf)(comm, file, " | PFFT_SHIFTED_OUT");
  if(ths->pfft_flags & PFFT_TUNE)
    PX(fprintf)(comm, file, " | PFFT_TUNE");
  else
    PX(fprintf)(comm, file, " | PFFT_NO_TUNE");
  if(ths->pfft_flags & PFFT_PRESERVE_INPUT)
    PX(fprintf)(comm, file, " | PFFT_PRESERVE_INPUT");
  if(ths->pfft_flags & PFFT_DESTROY_INPUT)
    PX(fprintf)(comm, file, " | PFFT_DESTROY_INPUT");
  if(ths->pfft_flags & PFFT_BUFFERED_INPLACE)
    PX(fprintf)(comm, file, " | PFFT_BUFFERED_INPLACE");
  PX(fprintf)(comm, file, "\n");

  if(ths->fftw_flags & FFTW_ESTIMATE)
    PX(fprintf)(comm, file, "\n%% fftw_flags == FFTW_ESTIMATE\n");
  else if(ths->fftw_flags & FFTW_PATIENT)
    PX(fprintf)(comm, file, "\n%% fftw_flags == FFTW_PATIENT\n");
  else if(ths->fftw_flags & FFTW_EXHAUSTIVE)
    PX(fprintf)(comm, file, "\n%% fftw_flags == FFTW_EXHAUSTIVE\n");
  else
    PX(fprintf)(comm, file, "\n%% fftw_flags == FFTW_MEASURE\n");
  
  PX(fprintf)(comm, file, "index(%d) = %d;  ", idx, idx);
  PX(fprintf)(comm, file, "procs(%d) = %d;  ", idx, size);

  PX(fprintf)(comm, file, "np_pfft(%d, 1:%d) = [", idx, ths->rnk_pm);
  for(int t=0; t<ths->rnk_pm; t++)
    PX(fprintf)(comm, file, "%td ", ths->np[t]);
  PX(fprintf)(comm, file, "];\n");

  PX(fprintf)(comm, file, "n_pfft(%d, 1:%d) = [", idx, ths->rnk_n);
  for(int t=0; t<ths->rnk_n; t++)
    PX(fprintf)(comm, file, "%td ", ths->n[t]);
  PX(fprintf)(comm, file, "];  ");

  PX(fprintf)(comm, file, "ni_pfft(%d, 1:%d) = [", idx, ths->rnk_n);
  for(int t=0; t<ths->rnk_n; t++)
    PX(fprintf)(comm, file, "%td ", ths->ni[t]);
  PX(fprintf)(comm, file, "];  ");

  PX(fprintf)(comm, file, "no_pfft(%d, 1:%d) = [", idx, ths->rnk_n);
  for(int t=0; t<ths->rnk_n; t++)
    PX(fprintf)(comm, file, "%td ", ths->no[t]);
  PX(fprintf)(comm, file, "];\n");
}

static int get_matlab_index(
    MPI_Comm comm
    )
{
  int size;
  R idx;

  MPI_Comm_size(comm, &size);
  idx = pfft_log2((R) size)+1;
  return (int) idx;
}

  


/*************************************
 * Functions to operate on PX(timer)
 ************************************/

PX(timer) PX(copy_timer)(
    const PX(timer) orig
    )
{
  PX(timer) copy = PX(mktimer)(orig->rnk_pm);
  
  copy->iter = orig->iter;
  copy->whole = orig->whole;
  for(int t=0; t<orig->rnk_trafo; t++)
    copy->trafo[t] = orig->trafo[t];
  for(int t=0; t<orig->rnk_remap; t++)
    copy->remap[t] = orig->remap[t];
  for(int t=0; t<2; t++)
    copy->remap_3dto2d[t] = orig->remap_3dto2d[t];

  return copy;
}

void PX(average_timer)(
     PX(timer) ths
     )
{
  if(ths->iter <= 0)
    return;
  
  ths->whole /= ths->iter;
  for(int t=0; t<ths->rnk_trafo; t++)
    ths->trafo[t] /= ths->iter;
  for(int t=0; t<ths->rnk_remap; t++)
    ths->remap[t] /= ths->iter;
  for(int t=0; t<2; t++)
    ths->remap_3dto2d[t] /= ths->iter;
  
  ths->iter = 1;
}

PX(timer) PX(add_timers)(
    const PX(timer) sum1, const PX(timer) sum2
    )
{
  PX(timer) res = PX(copy_timer)(sum1);
  
  res->iter += sum2->iter;
  res->whole += sum2->whole;
  for(int t=0; t<res->rnk_trafo; t++)
    res->trafo[t] += sum2->trafo[t];
  for(int t=0; t<res->rnk_remap; t++)
    res->remap[t] += sum2->remap[t];
  for(int t=0; t<2; t++)
    res->remap_3dto2d[t] += sum2->remap_3dto2d[t];

  return res;
}

PX(timer) PX(reduce_max_timer)(
    const PX(timer) ths, MPI_Comm comm
    )
{
  PX(timer) ths_max;
  double *times = PX(convert_timer2vec)(ths);
  double *times_max = (double*) malloc(sizeof(double) * length(ths));

  MPI_Reduce(times, times_max, length(ths), MPI_DOUBLE, MPI_MAX, 0, comm);
  ths_max = PX(convert_vec2timer)(times_max);

  free(times); free(times_max);
  return ths_max;
}

double* PX(convert_timer2vec)(
    const PX(timer) ths
    )
{
  int m=0;
  double *times = (double*) malloc(sizeof(double) * length(ths));

  times[m++] = (double) ths->rnk_pm;
  times[m++] = (double) ths->rnk_trafo;
  times[m++] = (double) ths->rnk_remap;
  times[m++] = (double) ths->iter;
  times[m++] = ths->whole;
  for(int t=0; t<ths->rnk_trafo; t++)
    times[m++] = ths->trafo[t];
  for(int t=0; t<ths->rnk_remap; t++)
    times[m++] = ths->remap[t];
  for(int t=0; t<2; t++)
    times[m++] = ths->remap_3dto2d[t];

  return times;
}

PX(timer) PX(convert_vec2timer)(
    const double *times
    )
{
  int m=3;
  PX(timer) ths = PX(mktimer)((int) times[0]);
  
  ths->iter  = (int) times[m++];
  ths->whole = times[m++];
  for(int t=0; t<ths->rnk_trafo; t++)
    ths->trafo[t] = times[m++];
  for(int t=0; t<ths->rnk_remap; t++)
    ths->remap[t] = times[m++];
  for(int t=0; t<2; t++)
    ths->remap_3dto2d[t] = times[m++];

  return ths;
}

static void reset_timer(
    PX(timer) ths
    )
{
  ths->iter = 0;
  ths->whole = 0;
  for(int t=0; t<ths->rnk_trafo; t++)
    ths->trafo[t] = 0;
  for(int t=0; t<ths->rnk_remap; t++)
    ths->remap[t] = 0;
  for(int t=0; t<2; t++)
    ths->remap_3dto2d[t] = 0;
}

PX(timer) PX(mktimer)(
    int rnk_pm
    )
{ 
  PX(timer) ths = (PX(timer)) malloc(sizeof(PX(timer_s)));
  
  ths->rnk_pm = rnk_pm;
  ths->rnk_trafo = 2*rnk_pm+2;
  ths->rnk_remap = 2*rnk_pm;
  ths->trafo = (double*) malloc(sizeof(double) * (size_t) ths->rnk_trafo);
  ths->remap = (double*) malloc(sizeof(double) * (size_t) ths->rnk_remap);

  reset_timer(ths);

  return ths;
}

static size_t length(
    const PX(timer) ths
    )
{
  /* +3 for rnk_pm, rnk_trafo, rnk_remap */
  /* +1 for number of iterations */
  /* +1 for whole trafo timer */
  /* +2 for remap_3dto2d[2] */
  return (size_t) (ths->rnk_trafo + ths->rnk_remap + 7);
}

void PX(destroy_timer)(
    PX(timer) ths
    )
{
  if(ths==NULL)
    return;
  
  if(ths->trafo != NULL)
    free(ths->trafo);
  if(ths->remap != NULL)
    free(ths->remap);

  free(ths);
}

static void fprint_average_timer_prefixed(
    MPI_Comm comm, FILE *file, const char *prefix,
    const PX(timer) ths, unsigned flags
    )
{
  int k=0, l=0;
  int idx = get_matlab_index(comm);
  PX(timer) mt = PX(reduce_max_timer)(ths, comm);
  PX(average_timer)(mt);
 
  if(flags & PFFTI_PRINT_TIMER_BASIC){
    PX(fprintf)(comm, file, "%s_iter(%d)    = %d;  ", prefix, idx, ths->iter);
    PX(fprintf)(comm, file, "%s(%d)   = %.3e;\n", prefix, idx, mt->whole);
  } else if(flags & PFFTI_PRINT_TIMER_ADV){
    /* print times of transposed out step */
    PX(fprintf)(comm, file, "%s_remap_3dto2d(%d, 2)   = %.3e;\n",
        prefix, idx, mt->remap_3dto2d[0]);
    for(int t=0; t<mt->rnk_pm; t++, k++, l++){
      PX(fprintf)(comm, file, "%s_trafo%d(%d, 2)   = %.3e;  ",
          prefix, k+1, idx, mt->trafo[k]);
      PX(fprintf)(comm, file, "%s_remap%d(%d, 1) = %.3e;\n",
          prefix, l+1, idx, mt->remap[l]);
    }
    PX(fprintf)(comm, file, "%s_trafo%d(%d, 2)   = %.3e;\n",
        prefix, k+1, idx, mt->trafo[k]);
    k++;

    /* print times of transposed in step */
    for(int t=0; t<mt->rnk_pm; t++, k++, l++){
      PX(fprintf)(comm, file, "%s_trafo%d(%d, 2)   = %.3e;  ",
          prefix, k+1, idx, mt->trafo[k]);
      PX(fprintf)(comm, file, "%s_remap%d(%d, 1) = %.3e;\n",
          prefix, l+1, idx, mt->remap[l]);
    }
    PX(fprintf)(comm, file, "%s_trafo%d(%d, 2)   = %.3e;\n",
        prefix, k+1, idx, mt->trafo[k]);
    PX(fprintf)(comm, file, "%s_remap_2dto3d(%d, 2)   = %.3e;\n",
        prefix, idx, mt->remap_3dto2d[1]);
  }

  PX(destroy_timer)(mt);
}





/****************************************
 * Functions to operate with files
 ***************************************/

static int file_exists(
    const char *name
    )
{
  FILE *f;
  
  f=fopen(name, "r");
  if(f != NULL)
    fclose(f);
  
  return (f != NULL);
}

static FILE* open_or_create_file_to_append(
    MPI_Comm comm, const char *name
    )
{
  FILE *f=fopen(name, "a+");
  if(f==NULL){
    PX(fprintf)(comm, stderr, "Error: Cannot open file %s.\n", name);
    exit(1);
  }
  return f;
}


