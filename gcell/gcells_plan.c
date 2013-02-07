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


static PX(gcplan) gc_mkplan(
    int rnk_n);
static void decompose(
    const INT *pn, const INT *block,
    int rnk_pm, const MPI_Comm *comms_pm,
    INT *local_n, INT *local_start);
static void get_coords(
    int rnk_pm, const MPI_Comm *comms_pm,
    int *coords_pm);
static int gc_RMA_applicable(
    int rnk_n, int rnk_pm, const INT *gc_below, const INT *gc_above);
static void vcopy_INT_transposed(
    int rnk_n, const INT* n,
    int rnk_pm, unsigned transp_flag,
    INT *cp);
static void create_neighbor_group(
    MPI_Comm comm_cart, MPI_Group *group_nbr);

/* FIXME: Generalize these functions to arbitray dimensions */
static void pad_array_with_zeros_outside_3d(
    const INT *n, const INT *zeros_below, const INT *zeros_above,
    INT tuple_size, R *real_data);
static void truncate_array_outside_3d(
    const INT *n, const INT *trunc_below, const INT *trunc_above,
    INT tuple_size, R *real_data);



INT PX(local_size_gc_internal)(
    int rnk_n, const INT *local_n, const INT *local_start,
    INT tuple, const INT *gc_below_user, const INT *gc_above_user,
    INT *local_ngc, INT *local_gc_start
    )
{
  INT mem = tuple;
  INT *gc_below, *gc_above;

  gc_below = PX(malloc_INT)(rnk_n);
  gc_above = PX(malloc_INT)(rnk_n);

  PX(evaluate_user_gcells)(
      rnk_n, gc_below_user, gc_above_user,
      gc_below, gc_above);

  PX(vadd_INT)(rnk_n, gc_below, local_n, local_ngc);
  PX(vadd_INT)(rnk_n, gc_above, local_ngc, local_ngc);
  PX(vsub_INT)(rnk_n, local_start, gc_below, local_gc_start);

  for(int t = 0; t < rnk_n; t++)
    mem *= gc_below[t] + local_n[t] + gc_above[t];
  
  free(gc_below); free(gc_above);
  return mem;
}


/* pfft_flag in {0, PFFT_TRANSPOSED_IN, PFFT_TRANSPOSED_OUT} */
PX(gcplan) PX(plan_rgc_internal)(
    int rnk_n, const INT *n, INT tuple_size, const INT *block_user,
    const INT *gc_below_user, const INT *gc_above_user,
    R *data, int rnk_pm, MPI_Comm *comms_pm, MPI_Comm comm_cart,
    unsigned gc_flags
    )
{
  int dims = 1, periods = 1, reorder = 0, rnk_self;
  INT ngc_total, *dummy;
  INT *gc_below, *gc_above;
  PX(gcplan) ths = NULL;

  MPI_Comm_rank(comm_cart, &rnk_self);

  gc_below = PX(malloc_INT)(rnk_n);
  gc_above = PX(malloc_INT)(rnk_n);

  PX(evaluate_user_gcells)(
      rnk_n, gc_below_user, gc_above_user,
      gc_below, gc_above);

  /* return NULL for zero gcells */
  int nothing_to_do = 1;
  for(int t=0; t<rnk_n; t++)
    if( (gc_below[t] != 0) || (gc_above[t] != 0) )
      nothing_to_do = 0;

  if(nothing_to_do){
    free(gc_below); free(gc_above);
    return NULL;
  }

  /* FIXME: generalize to arbitray dimensions */
  if(rnk_n != 3){
    PX(fprintf)(comm_cart, stderr,
        "Error: Gcell send only available for three dimensions.\n");
    return NULL;
  }

  ths = gc_mkplan(rnk_n);

  ths->data = data;
  ths->tuple = tuple_size;

  /* save arrays in transposed order */
  vcopy_INT_transposed(rnk_n, n, rnk_pm, gc_flags, ths->n);
  vcopy_INT_transposed(rnk_n, gc_below, rnk_pm, gc_flags, ths->gc_below);
  vcopy_INT_transposed(rnk_n, gc_above, rnk_pm, gc_flags, ths->gc_above);
 
  /* calculate 'blk' from transposed array */ 
  PX(evaluate_user_block_size)(rnk_pm, ths->n, block_user, comms_pm, ths->blk);
  for(int t=rnk_pm; t<rnk_n; t++)
    ths->blk[t] = ths->n[t];

  /* calculate 'loc_n' from transposed array */ 
  dummy = PX(malloc_INT)(rnk_pm);
  decompose(ths->n, ths->blk, rnk_pm, comms_pm,
      ths->loc_n, dummy);
  for(int t=rnk_pm; t<rnk_n; t++)
    ths->loc_n[t] = ths->n[t];
  free(dummy);

  /* calculate 'ngc' from transposed arrays */
  for(int t=0; t<rnk_n; t++)
    ths->ngc[t] = ths->gc_below[t] + ths->loc_n[t] + ths->gc_above[t];

  /* procmesh remains the same for transposed layout */
  for(int t=0; t<rnk_pm; t++){
    MPI_Cart_shift(comm_cart, t, 1, &ths->rnk_prec[t], &ths->rnk_succ[t]);
    MPI_Comm_dup(comms_pm[t], &ths->comms_pm[t]);
    MPI_Comm_size(comms_pm[t], &ths->np[t]);
  }
  for(int t=rnk_pm; t<rnk_n; t++){
    ths->rnk_prec[t] = ths->rnk_succ[t] = rnk_self;
    MPI_Cart_create(MPI_COMM_SELF, 1, &dims, &periods, reorder, &ths->comms_pm[t]);
    ths->np[t] = 1;
  }

  /* create memory window for RMA algorithm */
  ngc_total = tuple_size;
  for(int t=0; t<rnk_n; t++)
    ngc_total *= ths->ngc[t];

  MPI_Win_create(data, (MPI_Aint) ((size_t) ngc_total * sizeof(R)), sizeof(R),
      MPI_INFO_NULL, comm_cart, &ths->win);
  MPI_Comm_dup(comm_cart, &ths->comm_cart);

  /* group all neighbor prcesses for RMA algorithm */
  create_neighbor_group(comm_cart, &ths->grp);
  
  if(gc_flags & PFFT_GC_RMA)
    if( !gc_RMA_applicable(rnk_n, rnk_pm, ths->gc_below, ths->gc_above) )
      PX(fprintf)(comm_cart, stderr,
          "Error: RMA Gcell algorithm doesn't support this data decompostion.\n");

  /* default: call RMA ghost cell send if possible */
  ths->alg_flag = PFFT_GC_RMA;
  if( !gc_RMA_applicable(rnk_n, rnk_pm, ths->gc_below, ths->gc_above) )
    ths->alg_flag = PFFT_GC_SENDRECV;
  if( gc_flags & PFFT_GC_SENDRECV )
    ths->alg_flag = PFFT_GC_SENDRECV;

  /* calculate neighboring local size of distributed dimensions */  
  for(int t=0; t<rnk_n; t++){
    ths->ngc_prec[t] = ths->gc_below[t] + ths->gc_above[t];
    ths->ngc_prec[t] += PX(local_block_size_shifted)(
        ths->n[t], ths->blk[t], -1, ths->comms_pm[t]);
    ths->ngc_succ[t] = ths->gc_below[t] + ths->gc_above[t];
    ths->ngc_succ[t] += PX(local_block_size_shifted)(
        ths->n[t], ths->blk[t], +1, ths->comms_pm[t]);
  }

  free(gc_below); free(gc_above);
  return ths;
}

static int gc_RMA_applicable(
    int rnk_n, int rnk_pm, const INT *gc_below, const INT *gc_above
    )
{
  /* The RMA ghostcell implementation only supports data
   * distributions where the first two dimensions
   * are distributed among a two-dimensional processor grid. */

  if(rnk_n != 3)
    return 0;
  if(rnk_pm != 2)
    return 0;
  if(gc_below[2] != 0)
    return 0;
  if(gc_above[2] != 0)
    return 0;

  return 1;
}


static void vcopy_INT_transposed(
    int rnk_n, const INT* n,
    int rnk_pm, unsigned transp_flag,
    INT *cp
    )
{
  for(int t=0; t<rnk_pm; t++)
    cp[t] = (transp_flag & PFFT_GC_TRANSPOSED) ?
      n[(t+1)%rnk_pm] : n[t];

  for(int t=rnk_pm; t<rnk_n; t++)
    cp[t] = n[t];
}

/* calculate block sizes from physical array size
 * for the 'rnk_pm' distributed dimensions */
static void decompose(
    const INT *pn, const INT *block,
    int rnk_pm, const MPI_Comm *comms_pm,
    INT *local_n, INT *local_start
    )
{
  int *coords_pm = PX(malloc_int)(rnk_pm);

  get_coords(rnk_pm, comms_pm,
      coords_pm);

  PX(decompose)(pn, block, rnk_pm, coords_pm,
      local_n, local_start);

  free(coords_pm);
}


static void get_coords(
    int rnk_pm, const MPI_Comm *comms_pm,
    int *coords_pm
    )
{
  int rnk;
  
  for(int t=0; t<rnk_pm; t++){
    MPI_Comm_rank(comms_pm[t], &rnk);
    MPI_Cart_coords(comms_pm[t], rnk, 1, &coords_pm[t]); 
  }
}


static PX(gcplan) gc_mkplan(
    int rnk_n
    )
{
  PX(gcplan) ths;

  ths = (PX(gcplan)) malloc(sizeof(gcplan_s));
  
  ths->rnk_n = rnk_n;

  /* allocate arrays */
  ths->n        = PX(malloc_INT)(rnk_n);
  ths->loc_n    = PX(malloc_INT)(rnk_n);
  ths->gc_below = PX(malloc_INT)(rnk_n);
  ths->gc_above = PX(malloc_INT)(rnk_n);
  ths->ngc      = PX(malloc_INT)(rnk_n);
  ths->ngc_prec = PX(malloc_INT)(rnk_n);
  ths->ngc_succ = PX(malloc_INT)(rnk_n);

  ths->np       = PX(malloc_int)(rnk_n);
  ths->rnk_prec = PX(malloc_int)(rnk_n);
  ths->rnk_succ = PX(malloc_int)(rnk_n);

  /* allocate mem for rnk_n comms, instead of only rnk_pm */
  ths->blk      = PX(malloc_INT)(rnk_n);
  ths->comms_pm = (MPI_Comm*) malloc(sizeof(MPI_Comm) * (size_t) rnk_n);

  ths->timer_exg = PX(gc_mktimer)();
  ths->timer_red = PX(gc_mktimer)();
  
  return ths;
}



void PX(rmplan_gc)(
    PX(gcplan) ths
    )
{
  /* plan was already destroyed or never initialized */
  if(ths == NULL)
    return;

  free(ths->n);
  free(ths->loc_n);
  free(ths->gc_below);
  free(ths->gc_above);
  free(ths->ngc);
  free(ths->ngc_prec);
  free(ths->ngc_succ);

  free(ths->np);
  free(ths->rnk_prec);
  free(ths->rnk_succ);

  free(ths->blk);
  for(int t=0; t < ths->rnk_n; t++)
    MPI_Comm_free(&ths->comms_pm[t]);
  free(ths->comms_pm);

  /* take care of MPI variables for RMA algorithm */
  MPI_Group_free(&ths->grp);
  MPI_Win_free(&ths->win);
  MPI_Comm_free(&ths->comm_cart);

  PX(destroy_gctimer)(ths->timer_exg);
  PX(destroy_gctimer)(ths->timer_red);
  
  /* free memory */
  free(ths);
  /* ths=NULL; would be senseless, since we can not change the pointer itself */
}


static void create_neighbor_group(
    MPI_Comm comm_cart, MPI_Group *group_nbr
    )
{
  int rnk_pm, *dims, *periods, *coords;
  int mpi_self, mpi_left, mpi_right;
  int *ranks_nbr, n_nbr = 0, include_self = 0;
  MPI_Group group_cart;

  MPI_Comm_rank(comm_cart, &mpi_self);
  MPI_Cartdim_get(comm_cart, &rnk_pm);
  MPI_Comm_group(comm_cart, &group_cart);

  dims = PX(malloc_int)(rnk_pm);
  periods = PX(malloc_int)(rnk_pm);
  coords = PX(malloc_int)(rnk_pm);

  MPI_Cart_get(comm_cart, rnk_pm, dims, periods, coords);
  
  ranks_nbr = PX(malloc_int)(2*rnk_pm);

  /* avoid multiple inlcudes of same ranks */
  for(int t=0; t<rnk_pm; t++){
    MPI_Cart_shift(comm_cart, t, 1, &mpi_left, &mpi_right);
    switch(dims[t]){
      case 1:  
        /* this proc is its own right and left neighbor */
        include_self = 1; break;
      case 2:
        /* right and left neighbor are the same */
        ranks_nbr[n_nbr++] = mpi_left; break;
      default: 
        /* include right and left neighbors */
        ranks_nbr[n_nbr++] = mpi_left;
        ranks_nbr[n_nbr++] = mpi_right;
    }
  }

  /* include myrank only once, if neccessary */
  if( include_self )
    ranks_nbr[n_nbr++] = mpi_self;

  MPI_Group_incl(group_cart, n_nbr, ranks_nbr, group_nbr);
  
  free(ranks_nbr);
  free(dims); free(periods); free(coords);
  MPI_Group_free(&group_cart);
}

void PX(exchange_gc)(
    PX(gcplan) ths
    )
{
  if(ths == NULL)
    return;
  
  ths->timer_exg->pad_zeros -= MPI_Wtime();
  pad_array_with_zeros_outside_3d(
      ths->loc_n, ths->gc_below, ths->gc_above, ths->tuple,
      ths->data);
  ths->timer_exg->pad_zeros += MPI_Wtime();
  
  ths->timer_exg->exchange -= MPI_Wtime();
  if(ths->alg_flag & PFFT_GC_RMA)
    PX(exchange_gc_RMA)(ths);
  else
    PX(exchange_gc_sendrecv)(ths);
  ths->timer_exg->exchange += MPI_Wtime();
}


void PX(reduce_gc)(
    PX(gcplan) ths
    )
{
  if(ths == NULL)
    return;

  ths->timer_red->exchange -= MPI_Wtime();
  if(ths->alg_flag & PFFT_GC_RMA)
    PX(reduce_gc_RMA)(ths);
  else
    PX(reduce_gc_sendrecv)(ths);
  ths->timer_red->exchange += MPI_Wtime();

  ths->timer_red->pad_zeros -= MPI_Wtime();
  truncate_array_outside_3d(
      ths->ngc, ths->gc_below, ths->gc_above, ths->tuple,
      ths->data);
  ths->timer_red->pad_zeros += MPI_Wtime();
}



static void pad_array_with_zeros_outside_3d(
    const INT *n, const INT *zeros_below, const INT *zeros_above,
    INT tuple_size, R *real_data
    )
{
  INT m, mo;
  INT n_start[3], n_end[3], no[3];
  INT k0, k1, k2, h;
  int k0_is_outside, k1_is_outside, k2_is_outside, index_is_outside;

  PX(vcopy_INT)(3, zeros_below, n_start);
  PX(vadd_INT)(3, n_start, n, n_end);
  PX(vadd_INT)(3, n_end, zeros_above, no);

  if( PX(equal_INT)(3, n, no) )
    return;

  m  = tuple_size * PX(prod_INT)(3, n) - 1;
  mo = tuple_size * PX(prod_INT)(3, no) - 1;
  for(k0 = no[0] - 1 ; k0 >= 0; k0--){
    k0_is_outside = (k0 < n_start[0]) || (n_end[0] <= k0);
    for(k1 = no[1] - 1 ; k1 >= 0; k1--){
      k1_is_outside = (k1 < n_start[1]) || (n_end[1] <= k1);
      for(k2 = no[2] - 1 ; k2 >= 0; k2--){
        k2_is_outside = (k2 < n_start[2]) || (n_end[2] <= k2);
        index_is_outside = k0_is_outside || k1_is_outside || k2_is_outside;
        for(h = 0; h < tuple_size; h++)
          real_data[mo--] = (index_is_outside) ? 0 : real_data[m--];
      }
    }
  }
}


static void truncate_array_outside_3d(
    const INT *n, const INT *trunc_below, const INT *trunc_above,
    INT tuple_size, R *real_data
    )
{
  INT n_total, mo, k0, k1, k2, h;
  INT no_start[3], no_end[3], no[3];

  PX(vcopy_INT)(3, trunc_below, no_start);
  PX(vsub_INT)(3, n, trunc_above, no_end);
  PX(vsub_INT)(3, no_end, trunc_below, no);

  if( PX(equal_INT)(3, n, no) )
    return;

  mo = 0;
  for(k0 = no_start[0]; k0 < no_end[0]; k0++){
    for(k1 = no_start[1]; k1 < no_end[1]; k1++){
      for(k2 = no_start[2]; k2 < no_end[2]; k2++){
        for(h = 0; h < tuple_size; h++)
          real_data[mo++] = real_data[h+tuple_size*(k2+n[2]*(k1+n[1]*k0))];
      }
    }
  }

  n_total = tuple_size * PX(prod_INT)(3, n);
  for(k0 = mo; k0 < n_total; k0++)
    real_data[k0] = 0;
}

