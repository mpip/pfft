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

/* Warning: This ghostcell implementation only supports data distributions where the first two dimensions
   are distributed among a two-dimensional processor grid. */

/* get corners from the following neighbors
      |
      v
-->oxxo
   x  x
   x  x
   oxxo<--
   ^
   |
*/

static void exchange_gcells(
    PX(gcplan) ths, unsigned gcflags);
static void reduce_gcells(
    PX(gcplan) ths, unsigned gcflags);
static void exchange_gcells_along_n0(
    PX(gcplan) ths, unsigned gcflags, int direction,
    INT gc_xstart, INT gc_target_xstart, INT gc_xget);
static void exchange_gcells_along_n1(
    PX(gcplan) ths, unsigned gcflag, int direction,
    INT gc_ystart, INT gc_target_ystart, INT gc_yget);
static void execute_gcell_exchange(
    INT gc_xstart, INT gc_ystart, INT gc_target_xstart, INT gc_target_ystart,
    INT gc_xget, INT gc_yget, const INT *ngc, const INT* target_ngc,
    INT tuple_size, int target_rank, MPI_Win win, unsigned gcflags,
    R *data);
static int sync_communication_finished(
    const INT *gc_below_remain, const INT *gc_above_remain, MPI_Comm comm);


void PX(exchange_gc_RMA)(
    PX(gcplan) ths
    )
{
  exchange_gcells(ths, PFFTI_GC_BORDERS);
  exchange_gcells(ths, PFFTI_GC_CORNERS);
}


void PX(reduce_gc_RMA)(
    PX(gcplan) ths
    )
{
  reduce_gcells(ths, PFFTI_GC_BORDERS);
  reduce_gcells(ths, PFFTI_GC_CORNERS);
}


static void exchange_gcells(
    PX(gcplan) ths, unsigned gcflags
    )
{
  const int d=2;
  INT gc_below_remain[d], gc_above_remain[d];
  INT gc_below_start[d], gc_above_start[d];
  INT gc_below_avail[d], gc_above_avail[d];
  INT gc_prev_start[d], gc_next_start[d];
  INT gc_below_get[d], gc_above_get[d];
  
  for(int t=0; t<2; t++){
    gc_below_remain[t] = ths->gc_below[t];
    gc_above_remain[t] = ths->gc_above[t];
    gc_below_start[t] = ths->gc_below[t];
    gc_above_start[t] = ths->gc_below[t] + ths->loc_n[t];
    gc_prev_start[t] = ths->gc_below[t] + PX(local_block_size_shifted)(
        ths->n[t], ths->blk[t], -1, ths->comms_pm[t]);
    gc_next_start[t] = ths->gc_below[t];
  }
  
  for(int shift = 1; 1; shift++ ){
    for(int t=0; t<2; t++){
      gc_below_avail[t] = PX(local_block_size_shifted)(
          ths->n[t], ths->blk[t], -shift, ths->comms_pm[t]);
      gc_above_avail[t] = PX(local_block_size_shifted)(
          ths->n[t], ths->blk[t], +shift, ths->comms_pm[t]);
      
      gc_below_get[t]   = MIN(gc_below_remain[t], gc_below_avail[t]);
      gc_above_get[t]   = MIN(gc_above_remain[t], gc_above_avail[t]);
  
      gc_below_start[t] -= gc_below_get[t];
      gc_prev_start[t]  -= gc_below_get[t];
    }
    
    /* exchange is implemented with MPI_Get */
    MPI_Win_fence(MPI_MODE_NOPRECEDE | MPI_MODE_NOPUT, ths->win);

//    MPI_Win_post(ths->grp, 0, ths->win);
//    MPI_Win_start(ths->grp, 0, ths->win);
    
    exchange_gcells_along_n0(
        ths, gcflags| PFFTI_GC_TRAFO, -1, gc_below_start[0], gc_prev_start[0], gc_below_get[0]);
    exchange_gcells_along_n0(
        ths, gcflags| PFFTI_GC_TRAFO, +1, gc_above_start[0], gc_next_start[0], gc_above_get[0]);
    exchange_gcells_along_n1(
        ths, gcflags| PFFTI_GC_TRAFO, -1, gc_below_start[1], gc_prev_start[1], gc_below_get[1]);
    exchange_gcells_along_n1(
        ths, gcflags| PFFTI_GC_TRAFO, +1, gc_above_start[1], gc_next_start[1], gc_above_get[1]);
    
    for(int t=0; t<2; t++){
      gc_above_start[t] += gc_above_get[t];
      gc_next_start[t]  += gc_above_get[t];
      
      gc_below_remain[t] -= gc_below_avail[t];
      gc_above_remain[t] -= gc_above_avail[t];
    }

    MPI_Win_fence(MPI_MODE_NOSTORE | MPI_MODE_NOSUCCEED , ths->win);

//    MPI_Win_complete(ths->win);
//    MPI_Win_wait(ths->win);
    
//    if( (gc_below_remain[0] <= 0) && (gc_above_remain[0] <= 0) )
//      if( (gc_below_remain[1] <= 0) && (gc_above_remain[1] <= 0) )
//        break;

    /* synchronized end of communication */
    if( sync_communication_finished(gc_below_remain, gc_above_remain, ths->comm_cart) )
      break;
  }
}


static void reduce_gcells(
    PX(gcplan) ths, unsigned gcflags
    )
{
  const int d=2;
  int shifts;
  INT gc_below_remain[d], gc_above_remain[d];
  INT gc_below_start[d], gc_above_start[d];
  INT gc_below_avail[d], gc_above_avail[d];
  INT gc_prev_start[d], gc_next_start[d];
  INT gc_below_put[d], gc_above_put[d];
  
  for(int t=0; t<2; t++){
    gc_below_remain[t] = ths->gc_below[t];
    gc_above_remain[t] = ths->gc_above[t];
  }
  
  for(shifts = 0; 1; shifts++){
//    if( (gc_below_remain[0] <= 0) && (gc_above_remain[0] <= 0) )
//      if( (gc_below_remain[1] <= 0) && (gc_above_remain[1] <= 0) )
//        break;

    /* synchronized end of communication */
    if( sync_communication_finished(gc_below_remain, gc_above_remain, ths->comm_cart) )
      break;
    
    for(int t=0; t<2; t++){
      gc_below_avail[t] = PX(local_block_size_shifted)(
          ths->n[t], ths->blk[t], -shifts, ths->comms_pm[t]);
      gc_above_avail[t] = PX(local_block_size_shifted)(
          ths->n[t], ths->blk[t], +shifts, ths->comms_pm[t]);
      
      gc_below_remain[t] -= gc_below_avail[t];
      gc_above_remain[t] -= gc_above_avail[t];
    }
  }
  
  for(int t=0; t<2; t++){
    gc_below_start[t] = 0;
    gc_above_start[t] = ths->ngc[t];
    gc_prev_start[t] = PX(local_block_size_shifted)(
        ths->n[t], ths->blk[t], -1, ths->comms_pm[t]);
    gc_next_start[t] = ths->gc_below[t] + ths->gc_above[t];
  }

  for(int shift = shifts; shift > 0; shift-- ){
    for(int t=0; t<2; t++){
      gc_below_avail[t] = PX(local_block_size_shifted)(
          ths->n[t], ths->blk[t], -shift, ths->comms_pm[t]);
      gc_above_avail[t] = PX(local_block_size_shifted)(
          ths->n[t], ths->blk[t], +shift, ths->comms_pm[t]);
      
      gc_below_remain[t] += gc_below_avail[t];
      gc_above_remain[t] += gc_above_avail[t];
      
      gc_below_put[t] = MAX(0, MIN(gc_below_remain[t], gc_below_avail[t]));
      gc_above_put[t] = MAX(0, MIN(gc_above_remain[t], gc_above_avail[t]));
      
      gc_above_start[t] -= gc_above_put[t];
      gc_next_start[t]  -= gc_above_put[t];
    }
   
    MPI_Win_fence(MPI_MODE_NOPRECEDE, ths->win);
    
//    MPI_Win_post(ths->grp, 0, ths->win);
//    MPI_Win_start(ths->grp, 0, ths->win);
    
    exchange_gcells_along_n0(
        ths, gcflags| PFFTI_GC_ADJOINT, -1, gc_below_start[0], gc_prev_start[0], gc_below_put[0]);
    exchange_gcells_along_n0(
        ths, gcflags| PFFTI_GC_ADJOINT, +1, gc_above_start[0], gc_next_start[0], gc_above_put[0]);
    exchange_gcells_along_n1(
        ths, gcflags| PFFTI_GC_ADJOINT, -1, gc_below_start[1], gc_prev_start[1], gc_below_put[1]);
    exchange_gcells_along_n1(
        ths, gcflags| PFFTI_GC_ADJOINT, +1, gc_above_start[1], gc_next_start[1], gc_above_put[1]);
    
    for(int t=0; t<2; t++){
      gc_below_start[t] += gc_below_put[t];
      gc_prev_start[t]  += gc_below_put[t];
    }
  
    MPI_Win_fence(MPI_MODE_NOSTORE | MPI_MODE_NOSUCCEED , ths->win);
 
//    MPI_Win_complete(ths->win);
//    MPI_Win_wait(ths->win);
  }
}

/* The end of communication shifts must be synchronized on all processes.
 * Otherwise deadlocks occur because ofunmatched window sychronization calls. */
static int sync_communication_finished(
    const INT *gc_below_remain, const INT *gc_above_remain, MPI_Comm comm
    )
{
  int finished = 0, all_finished;
  
  if( (gc_below_remain[0] <= 0) && (gc_above_remain[0] <= 0) )
    if( (gc_below_remain[1] <= 0) && (gc_above_remain[1] <= 0) )
      finished = 1;

  MPI_Allreduce(&finished, &all_finished, 1, MPI_INT, MPI_MIN, comm);
  return all_finished;
}




static void exchange_gcells_along_n0(
    PX(gcplan) ths, unsigned gcflags, int direction,
    INT gc_xstart, INT gc_target_xstart, INT gc_xget
    )
{
  int target_rank, dim=0;
  INT *target_ngc, gc_ystart, gc_target_ystart, gc_yget;
  
  if(gc_xget <= 0)
    return;
  
  target_ngc = PX(malloc_INT)(ths->rnk_n);
  for(int t=0; t<ths->rnk_n; t++)
    target_ngc[t] = ths->loc_n[t];

  target_rank     = (direction == -1) ? ths->rnk_prec[dim] : ths->rnk_succ[dim];
  target_ngc[dim] = (direction == -1) ? ths->ngc_prec[dim] : ths->ngc_succ[dim];
  
  if(gcflags & PFFTI_GC_CORNERS){
    gc_ystart = gc_target_ystart = (direction == -1) ? ths->gc_below[1] + ths->loc_n[1] : 0;
    gc_yget = (direction == -1) ? ths->gc_above[1] : ths->gc_below[1];
  } else {
    gc_ystart = gc_target_ystart = ths->gc_below[1];
    gc_yget = ths->loc_n[1];
  }
  
  execute_gcell_exchange(
    gc_xstart, gc_ystart, gc_target_xstart, gc_target_ystart, gc_xget, gc_yget,
    ths->ngc, target_ngc, ths->tuple, target_rank, ths->win, gcflags,
    ths->data);
}


static void exchange_gcells_along_n1(
    PX(gcplan) ths, unsigned gcflags, int direction,
    INT gc_ystart, INT gc_target_ystart, INT gc_yget
    )
{
  int target_rank, dim=1;
  INT *target_ngc, gc_xstart, gc_target_xstart, gc_xget;
  
  if(gc_yget <= 0)
    return;

  target_ngc = PX(malloc_INT)(ths->rnk_n);
  for(int t=0; t<ths->rnk_n; t++)
    target_ngc[t] = ths->loc_n[t];
  
  target_rank     = (direction == -1) ? ths->rnk_prec[dim] : ths->rnk_succ[dim];
  target_ngc[dim] = (direction == -1) ? ths->ngc_prec[dim] : ths->ngc_succ[dim];

  if(gcflags & PFFTI_GC_CORNERS){
    gc_xstart = gc_target_xstart = (direction == -1) ? 0 : ths->gc_below[0] + ths->loc_n[0];
    gc_xget = (direction == -1) ? ths->gc_below[0] : ths->gc_above[0];
  } else {
    gc_xstart = gc_target_xstart = ths->gc_below[0];
    gc_xget = ths->loc_n[0];
  }
  
  execute_gcell_exchange(
    gc_xstart, gc_ystart, gc_target_xstart, gc_target_ystart, gc_xget, gc_yget,
    ths->ngc, target_ngc, ths->tuple, target_rank, ths->win, gcflags,
    ths->data);

  free(target_ngc);
}


static void execute_gcell_exchange(
    INT gc_xstart, INT gc_ystart, INT gc_target_xstart, INT gc_target_ystart,
    INT gc_xget, INT gc_yget, const INT *ngc, const INT* target_ngc,
    INT tuple_size, int target_rank, MPI_Win win, unsigned gcflags,
    R *data
    )
{
  int num_recv, num_send; /* FIXME: 64-bit portability issue, but MPI forces to use int */
  INT offset_recv, offset_send;
  
  for(INT k0 = 0; k0 < gc_xget; k0++){
    offset_recv = MACRO_PLAIN_INDEX_3D(
        k0 + gc_xstart,        gc_ystart,        0, ngc) * tuple_size;
    offset_send = MACRO_PLAIN_INDEX_3D(
        k0 + gc_target_xstart, gc_target_ystart, 0, target_ngc) * tuple_size;
    num_recv = num_send = (int) (gc_yget * ngc[2] * tuple_size);
    
    if(gcflags & PFFTI_GC_ADJOINT)
      MPI_Accumulate(data + offset_recv, num_recv, PFFT_MPI_REAL_TYPE,
          target_rank, offset_send, num_send, PFFT_MPI_REAL_TYPE, MPI_SUM, win);
    else
      MPI_Get(data + offset_recv, num_recv, PFFT_MPI_REAL_TYPE,
          target_rank, offset_send, num_send, PFFT_MPI_REAL_TYPE, win);
  }
}

