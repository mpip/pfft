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
#include "util.h"

/* use MPI_Datatype to send ghostcells */
/* This option is more comfortable for programming and avoids extra buffers for send/recv.
 * However the buffered send/recv is slightly faster. */
#define PFFT_GC_USE_MPI_DATATYPE 0

/* time parts of the ghostcell send */
#define PFFT_ENABLE_GC_TIMER 0

/* definition of timing macros */
#if PFFT_ENABLE_GC_TIMER 
#define PFFT_GC_INIT_TIMING(comm) \
  int tm_rank; \
  MPI_Comm_rank(comm, &tm_rank); \
  double tm_timer, tm_global_timer;
#define PFFT_GC_START_TIMING() \
  tm_timer = -MPI_Wtime();
#define PFFT_GC_FINISH_TIMING(comm, str) \
  tm_timer += MPI_Wtime(); \
  MPI_Reduce(&tm_timer, &tm_global_timer, 1, MPI_DOUBLE, MPI_MAX, 0, comm); \
  if(!tm_rank) printf("PFFT_GC_TIMING: %s takes %e s\n", str, tm_global_timer);
#else
#define PFFT_GC_INIT_TIMING(comm)
#define PFFT_GC_START_TIMING()
#define PFFT_GC_FINISH_TIMING(comm, str)
#endif


static void exchange_gcells_along_one_dim(
    PX(gcplan) ths, int dim);
static void reduce_gcells_along_one_dim(
    PX(gcplan) ths, int dim);
static void exchange_gcells_above_along_one_dim(
    PX(gcplan) ths, int dim);
static void exchange_gcells_below_along_one_dim(
    PX(gcplan) ths, int dim);
static void reduce_gcells_above_along_one_dim(
    PX(gcplan) ths, int dim);
static void reduce_gcells_below_along_one_dim(
    PX(gcplan) ths, int dim);
static void sendrecv_gcells_along_one_dim(
    PX(gcplan) ths, int dim, int dir,
    INT numSend, INT sendOffset, INT numRecv, INT recvOffset);
static void addsendrecv_gcells_along_one_dim(
    PX(gcplan) ths, int dim, int dir,
    INT numSend, INT sendOffset, INT numRecv, INT recvOffset);
#if PFFT_GC_USE_MPI_DATATYPE
static void isend_slices(
    R *data, int rnk_n, const INT *localArraySize, int dim, int dir, INT tupleSize,
    INT numSend, INT sendOffset, MPI_Comm commCart1d, 
    MPI_Request *request);
static void irecv_slices(
    R *data, int rnk_n, const INT *localArraySize, int dim, int dir, INT tupleSize,
    INT numRecv, INT recvOffset, MPI_Comm commCart1d, MPI_Request *request);
static void addrecv_slices(
    R *data, int rnk_n, const INT *localArraySize, int dim, int dir, INT tupleSize,
    INT numRecv, INT recvOffset, MPI_Comm commCart1d);

static void create_mpi_datatype_slices(
    int rnk_n, const INT *n, INT tuple,
    int dim, INT num_slices, INT offset,
    MPI_Datatype *newtype);

#endif
/* FIXME: Generalize to arbitrary dimesions
 * Idea: canonicalize to three dims: n0 * slices * n1 */
static void add_buffer_to_slices(
    R *data, R *buffer,
    const INT *localArraySize3d, INT tupleSize,
    int dim, INT numSlices, INT offset);

#if !PFFT_GC_USE_MPI_DATATYPE
static void copy_slices_to_buffer(
    R *data, R *buffer,
    const INT *localArraySize3d, INT tupleSize,
    int dim, INT numSlices, INT offset);
static void copy_buffer_to_slices(
    R *data, R *buffer,
    const INT *localArraySize3d, INT tupleSize,
    int dim, INT numSlices, INT offset);
static INT calculate_buffer_size(
    int rnk_n, const INT *n, INT tuple, INT num_slices, int dim);
static R* allocate_buffer(
    int rnk_n, const INT *n, INT tuple,
    int dim, INT num_slices,
    INT *buffer_size);
static void isend_packed_slices(
    int dir, INT buf_size, R *buffer, MPI_Comm commCart1d,
    MPI_Request *request);
static void irecv_packed_slices(
    int dir, INT buf_size, MPI_Comm commCart1d,
    MPI_Request *request, R *buffer);
static void pack_slices(
    INT buffer_size, int rnk_n, const INT *n, INT tuple,
    int dim, INT num_slices, INT offset, R *data,
    R *buffer); 
static void unpack_slices(
    int rnk_n, const INT *n, INT tuple,
    int dim, INT num_slices, INT offset,
    INT buffer_size, R *buffer,
    unsigned add_buffer,
    R *data);
#endif





void PX(exchange_gc_sendrecv)(
    PX(gcplan) ths
    )
{
  for(int dim=0; dim<ths->rnk_n; dim++)
    exchange_gcells_along_one_dim(ths, dim);
}


void PX(reduce_gc_sendrecv)(
    PX(gcplan) ths
    )
{
  for(int dim=0; dim<ths->rnk_n; dim++)
    reduce_gcells_along_one_dim(ths, dim);
}


static void exchange_gcells_along_one_dim(
    PX(gcplan) ths, int dim
    )
{
  exchange_gcells_below_along_one_dim(ths, dim);
  exchange_gcells_above_along_one_dim(ths, dim);
}


static void reduce_gcells_along_one_dim(
    PX(gcplan) ths, int dim
    )
{
  reduce_gcells_below_along_one_dim(ths, dim);
  reduce_gcells_above_along_one_dim(ths, dim);
}


static void exchange_gcells_above_along_one_dim(
    PX(gcplan) ths, int dim
    )
{
  INT numSendLeftOver, numRecvLeftOver, sendOffset, recvOffset;
  INT numCurrentSend, numCurrentRecv, numSendAvail, numRecvAvail;
  INT localArrayStart, localArrayEnd;
  INT globalArraySize = ths->n[dim], blockSize = ths->blk[dim];

  localArrayStart  = localArrayEnd = ths->gc_below[dim];
  localArrayEnd   += ths->loc_n[dim];

#if PFFT_DEBUG_GHOSTCELLS
  int myrank; MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  R rsum, grsum;
  INT start[3], end[3], ngc[3], k;
  INT l0, l1, l2, l, m0, m1, m2, m, n0, n1, n2, n;
  for(int t=0; t<3; t++)
    ngc[t] = ths->gc_below[t] + ths->loc_n[t] + ths->gc_above[t];
#endif

#if PFFT_DEBUG_GHOSTCELLS
  for(int t=0; t<3; t++){
    start[t] = ths->gc_below[t];
    end[t] = ths->gc_below[t] + ths->loc_n[t];
  }

  rsum = 0.0;
  for(INT k0=start[0]; k0<end[0]; k0++){
    for(INT k1=start[1]; k1<end[1]; k1++){
      for(INT k2=start[2]; k2<end[2]; k2++){
        for(INT t=0; t<ths->tuple; t++){
          k = t + ths->tuple*(k2 + ngc[2]*(k1 + ngc[1]*k0));
          rsum += fabs(ths->data[k]);
        }
      }
    }
  }
  MPI_Reduce(&rsum, &grsum, 1, PFFT_MPI_REAL_TYPE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(!myrank) fprintf(stderr, "PFFT GC-Send Above: dim = %d, Sum of initial data: %e\n", dim, grsum);
#endif

#if PFFT_DEBUG_GHOSTCELLS
  for(int t=0; t<3; t++){
    start[t] = (t<=dim) ? 0 : ths->gc_below[t];
    end[t] = ths->gc_below[t] + ths->loc_n[t];
    end[t] += (t<dim) ? ths->gc_above[t] : 0;
  }

  rsum = 0.0;
  for(INT k0=start[0]; k0<end[0]; k0++){
    for(INT k1=start[1]; k1<end[1]; k1++){
      for(INT k2=start[2]; k2<end[2]; k2++){
        for(INT t=0; t<ths->tuple; t++){
          k = t + ths->tuple*(k2 + ngc[2]*(k1 + ngc[1]*k0));
          rsum += fabs(ths->data[k]);
        }
      }
    }
  }
  MPI_Reduce(&rsum, &grsum, 1, PFFT_MPI_REAL_TYPE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(!myrank) fprintf(stderr, "PFFT GC-Send Above: dim = %d, Sum of data before gcsend: %e\n", dim, grsum);
#endif
#if PFFT_DEBUG_GHOSTCELLS
  if(!myrank) fprintf(stderr, "PFFT GC-Send Above: dim = %d, Before gcsend: loc_n = [%td, %td, %td], start = [%td, %td, %td], end = [%td, %td, %td], ngc = [%td, %td, %td], tuple = %td\n",
      dim, ths->loc_n[0], ths->loc_n[1], ths->loc_n[2], start[0], start[1], start[2], end[0], end[1], end[2], ngc[0], ngc[1], ngc[2], ths->tuple);
#endif
#if PFFT_DEBUG_GHOSTCELLS
  l0=ths->gc_below[0]; l1=ths->gc_below[1]; l2=ths->gc_below[2];
  m0=l0+ths->gc_above[0]-1, m1=l1; m2=l2;
  n0=0, n1=l1; n2=l2;
  l=l2 + ngc[2]*( l1+ngc[1]*l0 );
  m=m2 + ngc[2]*( m1+ngc[1]*m0 );
  n=n2 + ngc[2]*( n1+ngc[1]*n0 );
  if(!myrank) fprintf(stderr, "PFFT GC-Send Above: dim = %d, Before gcsend: data[%td, %td, %td] = %e + I* %e, data[%td, %td, %td] = %e + I* %e, data[%td, %td, %td] = %e + I* %e\n",
      dim, n0, n1, n2, ths->data[2*n], ths->data[2*n+1], l0, l1, l2, ths->data[2*l], ths->data[2*l+1], m0, m1, m2, ths->data[2*m], ths->data[2*m+1]); 
#endif

  sendOffset = localArrayStart;
  recvOffset = localArrayEnd;
  numSendLeftOver = numRecvLeftOver = ths->gc_above[dim];
  for(int shift=0; (numSendLeftOver > 0) || (numRecvLeftOver > 0); shift++){
    numSendAvail = PX(local_block_size_shifted)(
        globalArraySize, blockSize, shift,   ths->comms_pm[dim]);
    numRecvAvail = PX(local_block_size_shifted)(
        globalArraySize, blockSize, shift+1, ths->comms_pm[dim]);

    numCurrentSend = MIN(numSendLeftOver, numSendAvail);
    numCurrentRecv = MIN(numRecvLeftOver, numRecvAvail);

    sendrecv_gcells_along_one_dim(
        ths, dim, -1, numCurrentSend, sendOffset, numCurrentRecv, recvOffset);

    sendOffset += numCurrentSend;
    recvOffset += numCurrentRecv;

    numSendLeftOver -= numCurrentSend;
    numRecvLeftOver -= numCurrentRecv;
  }

#if PFFT_DEBUG_GHOSTCELLS
  for(int t=0; t<3; t++){
    start[t] = (t<=dim) ? 0 : ths->gc_below[t];
    end[t] = ths->gc_below[t] + ths->loc_n[t];
    end[t] += (t<=dim) ? ths->gc_above[t] : 0;
  }

  rsum = 0.0;
  for(INT k0=start[0]; k0<end[0]; k0++){
    for(INT k1=start[1]; k1<end[1]; k1++){
      for(INT k2=start[2]; k2<end[2]; k2++){
        for(INT t=0; t<ths->tuple; t++){
          k = t + ths->tuple*(k2 + ngc[2]*(k1 + ngc[1]*k0));
          rsum += fabs(ths->data[k]);
        }
      }
    }
  }
  MPI_Reduce(&rsum, &grsum, 1, PFFT_MPI_REAL_TYPE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(!myrank) fprintf(stderr, "PFFT GC-Send Above: dim = %d, Sum of data after gcsend: %e\n", dim, grsum);
#endif
#if PFFT_DEBUG_GHOSTCELLS
  if(!myrank) fprintf(stderr, "PFFT GC-Send Above: dim = %d, After gcsend: loc_n = [%td, %td, %td], start = [%td, %td, %td], end = [%td, %td, %td], tuple = %td\n",
      dim, ths->loc_n[0], ths->loc_n[1], ths->loc_n[2], start[0], start[1], start[2], end[0], end[1], end[2], ths->tuple);
#endif
#if PFFT_DEBUG_GHOSTCELLS
  l0=ths->gc_above[0]-1; l1=ths->gc_below[1]; l2=ths->gc_below[2];
  m0=ths->gc_below[0] + ths->loc_n[0], m1=l1; m2=l2;
  n0=0, n1=l1; n2=l2;
  l=l2 + ngc[2]*( l1+ngc[1]*l0 );
  m=m2 + ngc[2]*( m1+ngc[1]*m0 );
  n=n2 + ngc[2]*( n1+ngc[1]*n0 );
  if(!myrank) fprintf(stderr, "PFFT GC-Send Above: dim = %d, After gcsend: data[%td, %td, %td] = %e + I* %e, data[%td, %td, %td] = %e + I* %e, data[%td, %td, %td] = %e + I* %e\n",
      dim, n0, n1, n2, ths->data[2*n], ths->data[2*n+1], l0, l1, l2, ths->data[2*l], ths->data[2*l+1], m0, m1, m2, ths->data[2*m], ths->data[2*m+1]); 
#endif
}


static void exchange_gcells_below_along_one_dim(
    PX(gcplan) ths, int dim
    )
{
  INT numSendLeftOver, numRecvLeftOver, sendOffset, recvOffset;
  INT numCurrentSend, numCurrentRecv, numSendAvail, numRecvAvail;
  INT localArrayStart, localArrayEnd;
  INT globalArraySize = ths->n[dim], blockSize = ths->blk[dim];

  localArrayStart  = localArrayEnd = ths->gc_below[dim];
  localArrayEnd   += ths->loc_n[dim];

#if PFFT_DEBUG_GHOSTCELLS
  int myrank; MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  R rsum, grsum;
  INT start[3], end[3], ngc[3], k;
  for(int t=0; t<3; t++)
    ngc[t] = ths->gc_below[t] + ths->loc_n[t] + ths->gc_above[t];
#endif

#if PFFT_DEBUG_GHOSTCELLS
  for(int t=0; t<3; t++){
    start[t] = ths->gc_below[t];
    end[t] = ths->gc_below[t] + ths->loc_n[t];
  }

  rsum = 0.0;
  for(INT k0=start[0]; k0<end[0]; k0++){
    for(INT k1=start[1]; k1<end[1]; k1++){
      for(INT k2=start[2]; k2<end[2]; k2++){
        for(INT t=0; t<ths->tuple; t++){
          k = t + ths->tuple*(k2 + ngc[2]*(k1 + ngc[1]*k0));
          rsum += fabs(ths->data[k]);
        }
      }
    }
  }
  MPI_Reduce(&rsum, &grsum, 1, PFFT_MPI_REAL_TYPE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(!myrank) fprintf(stderr, "PFFT GC-Send Below: dim = %d, Sum of initial data: %e\n", dim, grsum);
#endif

#if PFFT_DEBUG_GHOSTCELLS
  for(int t=0; t<3; t++){
    start[t] = (t<dim) ? 0 : ths->gc_below[t];
    end[t] = ths->gc_below[t] + ths->loc_n[t];
    end[t] += (t<dim) ? ths->gc_above[t] : 0;
  }

  rsum = 0.0;
  for(INT k0=start[0]; k0<end[0]; k0++){
    for(INT k1=start[1]; k1<end[1]; k1++){
      for(INT k2=start[2]; k2<end[2]; k2++){
        for(INT t=0; t<ths->tuple; t++){
          k = t + ths->tuple*(k2 + ngc[2]*(k1 + ngc[1]*k0));
          rsum += fabs(ths->data[k]);
        }
      }
    }
  }
  MPI_Reduce(&rsum, &grsum, 1, PFFT_MPI_REAL_TYPE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(!myrank) fprintf(stderr, "PFFT GC-Send Below: dim = %d, Sum of data before gcsend: %e\n", dim, grsum);
#endif
#if PFFT_DEBUG_GHOSTCELLS
  if(!myrank) fprintf(stderr, "PFFT GC-Send Below: dim = %d, Before gcsend: loc_n = [%td, %td, %td], start = [%td, %td, %td], end = [%td, %td, %td], ngc = [%td, %td, %td], tuple = %td\n",
      dim, ths->loc_n[0], ths->loc_n[1], ths->loc_n[2], start[0], start[1], start[2], end[0], end[1], end[2], ngc[0], ngc[1], ngc[2], ths->tuple);
#endif

  sendOffset = localArrayEnd;
  recvOffset = localArrayStart;
  numSendLeftOver = numRecvLeftOver = ths->gc_below[dim];
  for(int shift=0; (numSendLeftOver > 0) || (numRecvLeftOver > 0); shift-- ){
    numSendAvail = PX(local_block_size_shifted)(
        globalArraySize, blockSize, shift,   ths->comms_pm[dim]);
    numRecvAvail = PX(local_block_size_shifted)(
        globalArraySize, blockSize, shift-1, ths->comms_pm[dim]);

    numCurrentSend = MIN(numSendLeftOver, numSendAvail);
    numCurrentRecv = MIN(numRecvLeftOver, numRecvAvail);

    sendOffset -= numCurrentSend;
    recvOffset -= numCurrentRecv;

    sendrecv_gcells_along_one_dim(
        ths, dim, +1, numCurrentSend, sendOffset, numCurrentRecv, recvOffset);

    numSendLeftOver -= numCurrentSend;
    numRecvLeftOver -= numCurrentRecv;
  }

#if PFFT_DEBUG_GHOSTCELLS
  for(int t=0; t<3; t++){
    start[t] = (t<=dim) ? 0 : ths->gc_below[t];
    end[t] = ths->gc_below[t] + ths->loc_n[t];
    end[t] += (t<dim) ? ths->gc_above[t] : 0;
  }

  rsum = 0.0;
  for(INT k0=start[0]; k0<end[0]; k0++){
    for(INT k1=start[1]; k1<end[1]; k1++){
      for(INT k2=start[2]; k2<end[2]; k2++){
        for(INT t=0; t<ths->tuple; t++){
          k = t + ths->tuple*(k2 + ngc[2]*(k1 + ngc[1]*k0));
          rsum += fabs(ths->data[k]);
        }
      }
    }
  }
  MPI_Reduce(&rsum, &grsum, 1, PFFT_MPI_REAL_TYPE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(!myrank) fprintf(stderr, "PFFT GC-Send Below: dim = %d, Sum of data after gcsend: %e\n", dim, grsum);
#endif
#if PFFT_DEBUG_GHOSTCELLS
  if(!myrank) fprintf(stderr, "PFFT GC-Send Below: dim = %d, After gcsend: loc_n = [%td, %td, %td], start = [%td, %td, %td], end = [%td, %td, %td], tuple = %td\n",
      dim, ths->loc_n[0], ths->loc_n[1], ths->loc_n[2], start[0], start[1], start[2], end[0], end[1], end[2], ths->tuple);
#endif
}


static void reduce_gcells_above_along_one_dim(
    PX(gcplan) ths, int dim
    )
{
  INT numSendLeftOver, numRecvLeftOver, sendOffset, recvOffset;
  INT numCurrentSend, numCurrentRecv, numSendAvail, numRecvAvail;
  INT localArrayStart, localArrayEnd, numSendPossible, numRecvPossible;
  INT globalArraySize = ths->n[dim], blockSize = ths->blk[dim];
  int numShifts;

  localArrayStart  = localArrayEnd = ths->gc_below[dim];
  localArrayEnd   += ths->loc_n[dim];

  /* corresponds to loop in exchange_gcells_above_along_one_dim 
   * with switched sender and receiver */
  numShifts = 0;
  numSendPossible = numRecvPossible = 0;
  numSendLeftOver = numRecvLeftOver = ths->gc_above[dim];
  for(int shift=0; (numSendLeftOver > 0) || (numRecvLeftOver > 0); shift++){
    numSendAvail = PX(local_block_size_shifted)(
        globalArraySize, blockSize, shift+1, ths->comms_pm[dim]);
    numRecvAvail = PX(local_block_size_shifted)(
        globalArraySize, blockSize, shift,   ths->comms_pm[dim]);

    numSendPossible += numSendAvail;
    numRecvPossible += numRecvAvail;

    numCurrentSend = MIN(numSendLeftOver, numSendAvail);
    numCurrentRecv = MIN(numRecvLeftOver, numRecvAvail);

    numSendLeftOver -= numCurrentSend;
    numRecvLeftOver -= numCurrentRecv;

    numShifts++;
  }

  sendOffset = localArrayEnd + ths->gc_above[dim];
  recvOffset = localArrayStart + ths->gc_above[dim];
  numSendLeftOver = numRecvLeftOver = ths->gc_above[dim];
  for(int shift = numShifts-1; shift >= 0; shift-- ){
    numSendPossible -= PX(local_block_size_shifted)(
        globalArraySize, blockSize, shift+1, ths->comms_pm[dim]);
    numRecvPossible -= PX(local_block_size_shifted)(
        globalArraySize, blockSize, shift,   ths->comms_pm[dim]);

    numCurrentSend = numSendLeftOver - MIN(numSendLeftOver, numSendPossible);
    numCurrentRecv = numRecvLeftOver - MIN(numRecvLeftOver, numRecvPossible);

    sendOffset -= numCurrentSend;
    recvOffset -= numCurrentRecv;

    addsendrecv_gcells_along_one_dim(
        ths, dim, +1, numCurrentSend, sendOffset, numCurrentRecv, recvOffset);

    numSendLeftOver -= numCurrentSend;
    numRecvLeftOver -= numCurrentRecv;

  }
}


static void reduce_gcells_below_along_one_dim(
    PX(gcplan) ths, int dim
    )
{
  INT numSendLeftOver, numRecvLeftOver, sendOffset, recvOffset;
  INT numCurrentSend, numCurrentRecv, numSendAvail, numRecvAvail;
  INT localArrayEnd, numSendPossible, numRecvPossible;
  INT globalArraySize = ths->n[dim], blockSize = ths->blk[dim];
  int numShifts;

  localArrayEnd  = ths->gc_below[dim];
  localArrayEnd += ths->loc_n[dim];

  /* corresponds to loop in exchange_gcells_above_along_one_dim 
   * with switched sender and receiver */
  numShifts = 0;
  numSendPossible = numRecvPossible = 0;
  numSendLeftOver = numRecvLeftOver = ths->gc_below[dim];
  for(int shift=0; (numSendLeftOver > 0) || (numRecvLeftOver > 0); shift--){
    numSendAvail = PX(local_block_size_shifted)(
        globalArraySize, blockSize, shift-1,   ths->comms_pm[dim]);
    numRecvAvail = PX(local_block_size_shifted)(
        globalArraySize, blockSize, shift, ths->comms_pm[dim]);

    numSendPossible += numSendAvail;
    numRecvPossible += numRecvAvail;

    numCurrentSend = MIN(numSendLeftOver, numSendAvail);
    numCurrentRecv = MIN(numRecvLeftOver, numRecvAvail);

    numSendLeftOver -= numCurrentSend;
    numRecvLeftOver -= numCurrentRecv;

    numShifts--;
  }

  sendOffset = 0;
  recvOffset = localArrayEnd - ths->gc_below[dim];
  numSendLeftOver = numRecvLeftOver = ths->gc_below[dim];
  for(int shift = numShifts+1; shift <= 0; shift++ ){
    numSendPossible -= PX(local_block_size_shifted)(
        globalArraySize, blockSize, shift-1,   ths->comms_pm[dim]);
    numRecvPossible -= PX(local_block_size_shifted)(
        globalArraySize, blockSize, shift, ths->comms_pm[dim]);

    numCurrentSend = numSendLeftOver - MIN(numSendLeftOver, numSendPossible);
    numCurrentRecv = numRecvLeftOver - MIN(numRecvLeftOver, numRecvPossible);

    addsendrecv_gcells_along_one_dim(
        ths, dim, -1, numCurrentSend, sendOffset, numCurrentRecv, recvOffset);

    sendOffset += numCurrentSend;
    recvOffset += numCurrentRecv;

    numSendLeftOver -= numCurrentSend;
    numRecvLeftOver -= numCurrentRecv;
  }
}


#if PFFT_GC_USE_MPI_DATATYPE
static void sendrecv_gcells_along_one_dim(
    PX(gcplan) ths, int dim, int dir,
    INT numSend, INT sendOffset, INT numRecv, INT recvOffset
    )
{
  MPI_Request mpi_req[2];

  PFFT_GC_INIT_TIMING(MPI_COMM_WORLD);
  PFFT_GC_START_TIMING();
  isend_slices(ths->data, ths->rnk_n, ths->ngc, dim, dir,
      ths->tuple, numSend, sendOffset, ths->comms_pm[dim], &mpi_req[0]);
  irecv_slices(ths->data, ths->rnk_n, ths->ngc, dim, dir,
      ths->tuple, numRecv, recvOffset, ths->comms_pm[dim], &mpi_req[1]);
  MPI_Waitall(2, mpi_req, MPI_STATUSES_IGNORE);
  PFFT_GC_FINISH_TIMING(MPI_COMM_WORLD, "send/recv");
}
#else
static void sendrecv_gcells_along_one_dim(
    PX(gcplan) ths, int dim, int dir,
    INT numSend, INT sendOffset, INT numRecv, INT recvOffset
    )
{
  INT sbuf_size=0, rbuf_size=0;
  R *sbuf=NULL, *rbuf=NULL;
  MPI_Request mpi_req[2];
  unsigned add_buffer=0;

  PFFT_GC_INIT_TIMING(MPI_COMM_WORLD);
  PFFT_GC_START_TIMING();
  sbuf = allocate_buffer(ths->rnk_n, ths->ngc, ths->tuple, dim, numSend,
      &sbuf_size);
  rbuf = allocate_buffer(ths->rnk_n, ths->ngc, ths->tuple, dim, numRecv,
      &rbuf_size);
  PFFT_GC_FINISH_TIMING(MPI_COMM_WORLD, "allocate buffers");
  
  PFFT_GC_START_TIMING();
  pack_slices(sbuf_size, ths->rnk_n, ths->ngc, ths->tuple, dim, numSend, sendOffset, ths->data,
      sbuf);
  PFFT_GC_FINISH_TIMING(MPI_COMM_WORLD, "pack_slices");
 
#if PFFT_ENABLE_GC_TIMER
  if(!tm_rank) printf("sbuf_size = %td, rbuf_size = %td\n", sbuf_size, rbuf_size); 
#endif
  PFFT_GC_START_TIMING();
  isend_packed_slices(dir, sbuf_size, sbuf, ths->comms_pm[dim],
      &mpi_req[0]);
  irecv_packed_slices(dir, rbuf_size, ths->comms_pm[dim],
      &mpi_req[1], rbuf);
  MPI_Waitall(2, mpi_req, MPI_STATUSES_IGNORE);
  PFFT_GC_FINISH_TIMING(MPI_COMM_WORLD, "send/recv");

  PFFT_GC_START_TIMING();
  unpack_slices(ths->rnk_n, ths->ngc, ths->tuple, dim, numRecv, recvOffset, rbuf_size, rbuf, add_buffer,
      ths->data);
  PFFT_GC_FINISH_TIMING(MPI_COMM_WORLD, "unpack_slices");

  if(sbuf != NULL) free(sbuf);
  if(rbuf != NULL) free(rbuf);
}
#endif

#if PFFT_GC_USE_MPI_DATATYPE
static void addsendrecv_gcells_along_one_dim(
    PX(gcplan) ths, int dim, int dir,
    INT numSend, INT sendOffset, INT numRecv, INT recvOffset
    )
{
  MPI_Request request;

  PFFT_GC_INIT_TIMING(MPI_COMM_WORLD);
  PFFT_GC_START_TIMING();
  isend_slices(ths->data, ths->rnk_n, ths->ngc, dim, dir,
      ths->tuple, numSend, sendOffset, ths->comms_pm[dim], &request);
  addrecv_slices(ths->data, ths->rnk_n, ths->ngc, dim, dir,
      ths->tuple, numRecv, recvOffset, ths->comms_pm[dim]);
  MPI_Wait(&request, MPI_STATUS_IGNORE);
  PFFT_GC_FINISH_TIMING(MPI_COMM_WORLD, "send/addrecv");
}
#else
static void addsendrecv_gcells_along_one_dim(
    PX(gcplan) ths, int dim, int dir,
    INT numSend, INT sendOffset, INT numRecv, INT recvOffset
    )
{
  INT sbuf_size=0, rbuf_size=0;
  R *sbuf=NULL, *rbuf=NULL;
  MPI_Request mpi_req[2];
  unsigned add_buffer=1;

  sbuf = allocate_buffer(ths->rnk_n, ths->ngc, ths->tuple, dim, numSend,
      &sbuf_size);
  rbuf = allocate_buffer(ths->rnk_n, ths->ngc, ths->tuple, dim, numRecv,
      &rbuf_size);
  
  pack_slices(sbuf_size, ths->rnk_n, ths->ngc, ths->tuple, dim, numSend, sendOffset, ths->data,
      sbuf);

  isend_packed_slices(dir, sbuf_size, sbuf, ths->comms_pm[dim],
      &mpi_req[0]);
  irecv_packed_slices(dir, rbuf_size, ths->comms_pm[dim],
      &mpi_req[1], rbuf);
  MPI_Waitall(2, mpi_req, MPI_STATUSES_IGNORE);

  unpack_slices(ths->rnk_n, ths->ngc, ths->tuple, dim, numRecv, recvOffset, rbuf_size, rbuf, add_buffer,
      ths->data);

  if(sbuf != NULL) free(sbuf);
  if(rbuf != NULL) free(rbuf);
}
#endif

#if PFFT_GC_USE_MPI_DATATYPE
static void isend_slices(
    R *data, int rnk_n, const INT *localArraySize, int dim, int dir, INT tupleSize,
    INT numSend, INT sendOffset, MPI_Comm commCart1d, 
    MPI_Request *request
    )
{
  int from, to, disp;
  INT numSendTotal=tupleSize;
  MPI_Datatype sendtype;

  for(int t=0; t<rnk_n; t++)
    numSendTotal *= (dim==t) ? numSend : localArraySize[t];
  
  /* Return immediately if nothing needs to be send */
  if(numSendTotal == 0){
    *request = MPI_REQUEST_NULL;
    return;
  }

  PFFT_GC_INIT_TIMING(MPI_COMM_WORLD);
  PFFT_GC_START_TIMING();
  create_mpi_datatype_slices(rnk_n, localArraySize, tupleSize, dim, numSend, sendOffset,
      &sendtype);
  MPI_Type_commit(&sendtype);
  PFFT_GC_FINISH_TIMING(MPI_COMM_WORLD, "create_mpi_datatype");

  disp = (dir<0) ? -1 : +1;
  MPI_Cart_shift(commCart1d, 0, disp, &from, &to);

  MPI_Isend(data, 1, sendtype, to, 7, commCart1d, request);

  MPI_Type_free(&sendtype);
}
#else
static void isend_packed_slices(
    int dir, INT buf_size, R *buffer, MPI_Comm commCart1d,
    MPI_Request *request
    )
{
  int from, to, disp;
  
  /* Return immediately if nothing needs to be send */
  if(buf_size == 0){
    *request = MPI_REQUEST_NULL;
    return;
  }
 
  disp = (dir<0) ? -1 : +1;
  MPI_Cart_shift(commCart1d, 0, disp, &from, &to);

  MPI_Isend(buffer, buf_size, PFFT_MPI_REAL_TYPE, to, 7, commCart1d, request);
}
#endif


#if PFFT_GC_USE_MPI_DATATYPE
static void irecv_slices(
    R *data, int rnk_n, const INT *localArraySize, int dim, int dir, INT tupleSize,
    INT numRecv, INT recvOffset, MPI_Comm commCart1d, MPI_Request *request
    )
{
  int from, to, disp;
  INT numRecvTotal=tupleSize;
  MPI_Datatype recvtype;
  
  for(int t=0; t<rnk_n; t++)
    numRecvTotal *= (dim==t) ? numRecv : localArraySize[t];
  
  /* Return immediately if nothing needs to be received */
  if(numRecvTotal == 0){
    *request = MPI_REQUEST_NULL;
    return;
  }

  PFFT_GC_INIT_TIMING(MPI_COMM_WORLD);
  PFFT_GC_START_TIMING();
  create_mpi_datatype_slices(rnk_n, localArraySize, tupleSize, dim, numRecv, recvOffset,
      &recvtype);
  MPI_Type_commit(&recvtype);
  PFFT_GC_FINISH_TIMING(MPI_COMM_WORLD, "create_mpi_datatype");

  disp = (dir<0) ? -1 : +1;
  MPI_Cart_shift(commCart1d, 0, -disp, &to, &from);

  MPI_Irecv(data, 1, recvtype, from, 7, commCart1d, request);

  MPI_Type_free(&recvtype);
}
#else
static void irecv_packed_slices(
    int dir, INT buf_size, MPI_Comm commCart1d,
    MPI_Request *request, R *buffer
    )
{
  int from, to, disp;
  
  /* Return immediately if nothing needs to be received */
  if(buf_size == 0){
    *request = MPI_REQUEST_NULL;
    return;
  }

  disp = (dir<0) ? -1 : +1;
  MPI_Cart_shift(commCart1d, 0, -disp, &to, &from);

  MPI_Irecv(buffer, buf_size, PFFT_MPI_REAL_TYPE, from, 7, commCart1d, request);
}
#endif



#if PFFT_GC_USE_MPI_DATATYPE
static void addrecv_slices(
    R *data, int rnk_n, const INT *localArraySize, int dim, int dir, INT tupleSize,
    INT numRecv, INT recvOffset, MPI_Comm commCart1d
    )
{
  int from, to, disp;
  INT numRecvTotal=tupleSize;
  R *buffer;

  for(int t=0; t<rnk_n; t++)
    numRecvTotal *= (dim==t) ? numRecv : localArraySize[t];

  /* Return immediately if nothing needs to be received */
  if(numRecvTotal == 0)
    return;

  disp = (dir<0) ? -1 : +1;
  MPI_Cart_shift(commCart1d, 0, -disp, &to, &from);

  buffer = (R*) malloc(sizeof(R) * (size_t) numRecvTotal);
  MPI_Recv(buffer, (int) numRecvTotal, PFFT_MPI_REAL_TYPE, from, 7, commCart1d, MPI_STATUS_IGNORE);

  add_buffer_to_slices(data, buffer, localArraySize, tupleSize, dim, numRecv, recvOffset);
  free(buffer);
}

static void create_mpi_datatype_slices(
    int rnk_n, const INT *n, INT tuple,
    int dim, INT num_slices, INT offset,
    MPI_Datatype *newtype
    )
{
  int t;
  int *sizes, *subsizes, *starts;
  MPI_Datatype oldtype=PFFT_MPI_REAL_TYPE;

  sizes = PX(malloc_int)(rnk_n+1);
  subsizes = PX(malloc_int)(rnk_n+1);
  starts = PX(malloc_int)(rnk_n+1);

  for(t=0; t<rnk_n; t++){
    sizes[t]    = (int) ( n[t] );
    subsizes[t] = (int) ( (t==dim) ? num_slices : n[t] );
    starts[t]   = (int) ( (t==dim) ? offset    : 0 );
  }
  sizes[t] = subsizes[t] = (int) (tuple);
  starts[t]= 0;

  /* FIXME: array_of_sizes, array_of_subsizes, array_of_starts are integer */
  /* TODO:  implement MPI_Type_create_subarray with INT similar MPI_Datatypes*/
  MPI_Type_create_subarray(rnk_n+1, sizes, subsizes, starts, MPI_ORDER_C, oldtype, newtype);

  free(sizes); free(subsizes); free(starts);
}

#endif

static void add_buffer_to_slices(
    R *data, R *buffer,
    const INT *localArraySize3d, INT tupleSize,
    int dim, INT numSlices, INT offset
    )
{
  const int d=3;
  INT n[d], nSub[d], nStart[d];
  INT k0, k1, k2, h, ind;

  for(int t=0; t<3; t++){
    n[t]      = localArraySize3d[t];
    nSub[t]   = (t==dim) ? numSlices : localArraySize3d[t];
    nStart[t] = (t==dim) ? offset    : 0;
  }

  ind = 0;
  for(k0 = nStart[0]; k0 < nStart[0] + nSub[0]; k0++){
    for(k1 = nStart[1]; k1 < nStart[1] + nSub[1]; k1++){
      for(k2 = nStart[2]; k2 < nStart[2] + nSub[2]; k2++){
        for(h = 0; h < tupleSize; h++){
          data[h+tupleSize*(k2+n[2]*(k1+n[1]*k0))] += buffer[ind++];
        }
      }
    }
  }
}


#if !PFFT_GC_USE_MPI_DATATYPE
static void copy_slices_to_buffer(
    R *data, R *buffer,
    const INT *localArraySize3d, INT tupleSize,
    int dim, INT numSlices, INT offset
    )
{
  const int d=3;
  INT n[d], nSub[d], nStart[d];
  INT k0, k1, k2, h, ind;

  for(int t=0; t<3; t++){
    n[t]      = localArraySize3d[t];
    nSub[t]   = (t==dim) ? numSlices : localArraySize3d[t];
    nStart[t] = (t==dim) ? offset    : 0;
  }

  ind = 0;
  for(k0 = nStart[0]; k0 < nStart[0] + nSub[0]; k0++){
    for(k1 = nStart[1]; k1 < nStart[1] + nSub[1]; k1++){
      for(k2 = nStart[2]; k2 < nStart[2] + nSub[2]; k2++){
        for(h = 0; h < tupleSize; h++, ind++){
          buffer[ind] = data[h+tupleSize*(k2+n[2]*(k1+n[1]*k0))];
        }
      }
    }
  }
}

static void copy_buffer_to_slices(
    R *data, R *buffer,
    const INT *localArraySize3d, INT tupleSize,
    int dim, INT numSlices, INT offset
    )
{
  const int d=3;
  INT n[d], nSub[d], nStart[d];
  INT k0, k1, k2, h, ind;

  for(int t=0; t<3; t++){
    n[t]      = localArraySize3d[t];
    nSub[t]   = (t==dim) ? numSlices : localArraySize3d[t];
    nStart[t] = (t==dim) ? offset    : 0;
  }

  ind = 0;
  for(k0 = nStart[0]; k0 < nStart[0] + nSub[0]; k0++){
    for(k1 = nStart[1]; k1 < nStart[1] + nSub[1]; k1++){
      for(k2 = nStart[2]; k2 < nStart[2] + nSub[2]; k2++){
        for(h = 0; h < tupleSize; h++){
          data[h+tupleSize*(k2+n[2]*(k1+n[1]*k0))] += buffer[ind++];
        }
      }
    }
  }
}

static INT calculate_buffer_size(
    int rnk_n, const INT *n, INT tuple, INT num_slices, int dim
    )
{
  INT buf_size=tuple;

  for(int t=0; t<rnk_n; t++)
    buf_size *= (dim==t) ? num_slices : n[t];

  return buf_size;
}

static R* allocate_buffer(
    int rnk_n, const INT *n, INT tuple,
    int dim, INT num_slices,
    INT *buffer_size
    )
{
  *buffer_size = calculate_buffer_size(rnk_n, n, tuple, num_slices, dim);
  if(*buffer_size == 0)
    return NULL;
  else
    return (R*) malloc(sizeof(R) * (*buffer_size));
}


/* Copy a block of size n[rnk_n] starting at offset[rnk_n] into
 * a contiguous buffer.  */
static void pack_slices(
    INT buffer_size, int rnk_n, const INT *n, INT tuple,
    int dim, INT num_slices, INT offset, R *data,
    R *buffer    
    )
{
  int t;
  INT *sizes, *subsizes, *starts;
 
  if(buffer_size == 0)
    return;

  if(rnk_n==3){
    copy_slices_to_buffer(data, buffer, n, tuple, dim, num_slices, offset);
  } else {
    /* TODO: this implementation is quite slow */

    /* calculate size and offset of data block */ 
    sizes = PX(malloc_INT)(rnk_n+1);
    subsizes = PX(malloc_INT)(rnk_n+1);
    starts = PX(malloc_INT)(rnk_n+1);
  
    for(t=0; t<rnk_n; t++){
      sizes[t]    = n[t];
      subsizes[t] = (t==dim) ? num_slices : n[t];
      starts[t]   = (t==dim) ? offset    : 0;
    }
    /* last dim corresponds to tuple size */
    sizes[t] = subsizes[t] = (tuple);
    starts[t]= 0;
  
    INT subind_lin, ind_lin;
    INT *subind_vec = PX(malloc_INT)(rnk_n+1);
    
    /* loop over linearized index in contiguous buffer */
    for(subind_lin=0; subind_lin<buffer_size; subind_lin++){
      
      /* calculate vectorized index in contiguous buffer */
      INT tmp = subind_lin;
      for(t=(rnk_n+1)-1; t>=0; t--){
        subind_vec[t] = tmp % subsizes[t];
        tmp /= subsizes[t];
      }
  
      /* calculate vectorized index in strided buffer */
      ind_lin = subind_vec[0] + starts[0];
      for(t=1; t<rnk_n+1; t++){
        ind_lin *= sizes[t];
        ind_lin += subind_vec[t] + starts[t];
      }
  
      /* copy strided data into contiguous buffer */
      buffer[subind_lin] = data[ind_lin];
    }
    
    free(subind_vec);
    free(sizes); free(subsizes); free(starts);
  }
}


/* Copy/Add a contiguous buffer into a block of size n[rnk_n] starting at offset[rnk_n].  */
static void unpack_slices(
    int rnk_n, const INT *n, INT tuple,
    int dim, INT num_slices, INT offset,
    INT buffer_size, R *buffer,
    unsigned add_buffer,
    R *data
    )
{
  int t;
  INT *sizes, *subsizes, *starts;
 
  if(buffer_size==0)
    return;

  if(rnk_n==3){
    if(add_buffer)
      add_buffer_to_slices(data, buffer, n, tuple, dim, num_slices, offset);
    else
      copy_buffer_to_slices(data, buffer, n, tuple, dim, num_slices, offset);
  } else {
    /* TODO: this implementation is quite slow */
    
    /* calculate size and offset of data block */ 
    sizes = PX(malloc_INT)(rnk_n+1);
    subsizes = PX(malloc_INT)(rnk_n+1);
    starts = PX(malloc_INT)(rnk_n+1);
  
    for(t=0; t<rnk_n; t++){
      sizes[t]    = n[t];
      subsizes[t] = (t==dim) ? num_slices : n[t];
      starts[t]   = (t==dim) ? offset    : 0;
    }
    /* last dim corresponds to tuple size */
    sizes[t] = subsizes[t] = (tuple);
    starts[t]= 0;
  
    INT subind_lin, ind_lin;
    INT *subind_vec = PX(malloc_INT)(rnk_n+1);

    /* loop over linearized index in contiguous buffer */
    for(subind_lin=0; subind_lin<buffer_size; subind_lin++){
      
      /* calculate vectorized index in contiguous buffer */
      INT tmp = subind_lin;
      for(t=(rnk_n+1)-1; t>=0; t--){
        subind_vec[t] = tmp % subsizes[t];
        tmp /= subsizes[t];
      }
  
      /* calculate vectorized index in strided buffer */
      ind_lin=subind_vec[0] + starts[0];
      for(t=1; t<rnk_n+1; t++){
        ind_lin *= sizes[t];
        ind_lin += subind_vec[t] + starts[t];
      }
  
      /* add/write contiguous buffer into strided data */
      if(add_buffer)
        data[ind_lin] += buffer[subind_lin];
      else
        data[ind_lin] = buffer[subind_lin];
    }

    free(subind_vec);
    free(sizes); free(subsizes); free(starts);
  }
}
#endif


