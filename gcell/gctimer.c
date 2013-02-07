/*
 * Copyright (c) 2010-2013 Michael Pippig
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

static void reset_gctimer(
    PX(gctimer) ths);
static void fprint_average_gctimer(
    MPI_Comm comm, FILE *file, PX(gcplan) ths, unsigned flags);
static void write_info_header(
    MPI_Comm comm, FILE *file);
static void write_run_specific_infos(
    MPI_Comm comm, FILE *file, PX(gcplan) ths);
static int get_matlab_index(
    MPI_Comm comm);
static void fprint_average_gctimer_prefixed(
    MPI_Comm comm, FILE *file, const char *prefix, PX(gctimer) ths,
    unsigned flags);
static int file_exists(
    const char *name);
static FILE* open_or_create_file_to_append(
    MPI_Comm comm, const char *name);

/*********************************************
 * Functions to operate on timer in PX(gcplan)
 ********************************************/

PX(gctimer) PX(get_gctimer_exg)(
    const PX(gcplan) ths
    )
{
  return PX(copy_gctimer)(ths->timer_exg);
}

PX(gctimer) PX(get_gctimer_red)(
    const PX(gcplan) ths
    )
{
  return PX(copy_gctimer)(ths->timer_red);
}

void PX(reset_gctimers)(
    PX(gcplan) ths
    )
{
  reset_gctimer(ths->timer_exg);
  reset_gctimer(ths->timer_red);
}



/*************************************
 * Functions to operate on PX(gctimer)
 ************************************/

void PX(print_average_gctimer)(
    const PX(gcplan) ths, MPI_Comm comm
    )
{
  write_info_header(comm, stdout);
  write_run_specific_infos(comm, stdout, ths);
  fprint_average_gctimer(comm, stdout, ths, PFFTI_PRINT_TIMER_BASIC);
}

void PX(print_average_gctimer_adv)(
    const PX(gcplan) ths, MPI_Comm comm
    )
{
  PX(print_average_gctimer)(ths, comm);
  fprint_average_gctimer(comm, stdout, ths, PFFTI_PRINT_TIMER_ADV);
}

void PX(write_average_gctimer)(
    const PX(gcplan) ths, const char *name, MPI_Comm comm
    )
{
  int newfile;
  FILE *f;
  
  newfile = !file_exists(name);
  f = open_or_create_file_to_append(comm, name);
  if( newfile )
    write_info_header(comm, f);
  
  write_run_specific_infos(comm, f, ths);
  fprint_average_gctimer(comm, f, ths, PFFTI_PRINT_TIMER_BASIC);
  
  fclose(f);
}

void PX(write_average_gctimer_adv)(
    const PX(gcplan) ths, const char *name, MPI_Comm comm
    )
{
  FILE *f;
  PX(write_average_gctimer)(ths, name, comm);
  f = open_or_create_file_to_append(comm, name);
  fprint_average_gctimer(comm, f, ths, PFFTI_PRINT_TIMER_ADV);
  fclose(f);
}

static void fprint_average_gctimer(
    MPI_Comm comm, FILE *file, PX(gcplan) ths, unsigned flags
    )
{
  fprint_average_gctimer_prefixed(comm, file, "gcells_exg", ths->timer_exg, flags);
  fprint_average_gctimer_prefixed(comm, file, "gcells_red", ths->timer_red, flags);
}

static void write_info_header(
    MPI_Comm comm, FILE *file
    )
{
  PX(fprintf)(comm, file, "\n%% n - FFT size\n");
  PX(fprintf)(comm, file, "%% np - process grid\n");
  PX(fprintf)(comm, file, "%% procs - number of processes\n");
  PX(fprintf)(comm, file, "%% index(i) = log(procs(i)) + 1\n");
}

static void write_run_specific_infos(
    MPI_Comm comm, FILE *file, PX(gcplan) ths
    )
{
  
  int size;
  int idx = get_matlab_index(comm);
  MPI_Comm_size(comm, &size);
  
  PX(fprintf)(comm, file, "\nindex(%d) = %d;  ", idx, idx);
  PX(fprintf)(comm, file, "procs(%d) = %d;  ", idx, size);
  PX(fprintf)(comm, file, "np_gc(%d, 1:3) = [%d %d %d];\n",
      idx, ths->np[0], ths->np[1], ths->np[2]);
  PX(fprintf)(comm, file, "n_gc(%d, 1:3) = [%td %td %td];  ",
      idx, ths->n[0], ths->n[1], ths->n[2]);
  PX(fprintf)(comm, file, "gc_below(%d, 1:3) = [%td %td %td];  ",
      idx, ths->gc_below[0], ths->gc_below[1], ths->gc_below[2]);
  PX(fprintf)(comm, file, "gc_above(%d, 1:3) = [%td %td %td];\n",
      idx, ths->gc_above[0], ths->gc_above[1], ths->gc_above[2]);
}

static int get_matlab_index(
    MPI_Comm comm
    )
{
  int size;
  R idx;

  MPI_Comm_size(comm, &size);
  idx = pfft_log2((R)size)+1;
  return (int) idx;
}



PX(gctimer) PX(copy_gctimer)(
    const PX(gctimer) orig
    )
{
  PX(gctimer) ths = PX(gc_mktimer)();
  
  ths->iter  = orig->iter;
  ths->whole = orig->whole;
  ths->pad_zeros = orig->pad_zeros;
  ths->exchange  = orig->exchange;
  
  return ths;
}

void PX(average_gctimer)(
     PX(gctimer) ths
     )
{
  if(ths->iter <= 0)
    return;

  ths->whole /= ths->iter;
  ths->pad_zeros /= ths->iter;
  ths->exchange /= ths->iter;
  ths->iter = 1;
}

PX(gctimer) PX(add_gctimers)(
    const PX(gctimer) sum1, const PX(gctimer) sum2
    )
{
  PX(gctimer) ths = PX(copy_gctimer)(sum1);
  
  ths->iter  += sum2->iter;
  ths->whole += sum2->whole;
  ths->pad_zeros += sum2->pad_zeros;
  ths->exchange  += sum2->exchange;

  return ths;
}

PX(gctimer) PX(reduce_max_gctimer)(
    const PX(gctimer) ths, MPI_Comm comm
    )
{
  double times[4], times_max[4];

  PX(convert_gctimer2vec)(ths, times);
  MPI_Reduce(times, times_max, 4, MPI_DOUBLE, MPI_MAX, 0, comm);

  return PX(convert_vec2gctimer)(times_max);
}

void PX(convert_gctimer2vec)(
    const PX(gctimer) gctimer,
    double *times
    )
{
  times[0] = (double) gctimer->iter;
  times[1] = gctimer->whole;
  times[2] = gctimer->pad_zeros;
  times[3] = gctimer->exchange;
}

PX(gctimer) PX(convert_vec2gctimer)(
    const double *times
    )
{
  PX(gctimer) ths = PX(gc_mktimer)();

  ths->iter      = (int) times[0];
  ths->whole     = times[1];
  ths->pad_zeros = times[2];
  ths->exchange  = times[3];

  return ths;
}

static void reset_gctimer(
    PX(gctimer) ths
    )
{
  ths->iter = 0;
  ths->whole = 0;
  ths->pad_zeros = 0;
  ths->exchange = 0;
}

PX(gctimer) PX(gc_mktimer)(
    void
    )
{
  PX(gctimer) ths = (PX(gctimer)) malloc(sizeof(PX(gctimer_s)));
  reset_gctimer(ths);
  return ths;  
}

void PX(destroy_gctimer)(
    PX(gctimer) ths
    )
{
  if(ths==NULL)
    return;
  
  free(ths);
}


static void fprint_average_gctimer_prefixed(
    MPI_Comm comm, FILE *file, const char *prefix, PX(gctimer) ths,
    unsigned flags
    )
{
  int idx = get_matlab_index(comm);
  PX(gctimer) mt;
  mt = PX(reduce_max_gctimer)(ths, comm);
  PX(average_gctimer)(mt);
  
  if(flags & PFFTI_PRINT_TIMER_BASIC){
    PX(fprintf)(comm, file, "%s_gc_iter(%d)      = %d;  ", prefix, idx, ths->iter);
    PX(fprintf)(comm, file, "%s_gcells(%d)       = %.3e;\n", prefix, idx, mt->whole);
  } else if(flags & PFFTI_PRINT_TIMER_ADV){
    PX(fprintf)(comm, file, "%s_gc_pad_zeros(%d) = %.3e;  ", prefix, idx, mt->pad_zeros);
    PX(fprintf)(comm, file, "%s_gc_exchange(%d)  = %.3e;\n", prefix, idx, mt->exchange);
  }
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

