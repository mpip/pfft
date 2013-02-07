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


/* PFFT util header file */
#ifndef PFFT_UTIL_H
#define PFFT_UTIL_H 1

#include "ipfft.h"

void PX(decompose)(
    const INT *pn, const INT* block,
    int rnk_pm, const int *coords_pm,
    INT *local_n, INT *local_start);
void PX(decompose_1d)(
    INT n, INT block_size, int which_block,
    INT *local_n, INT *local_n_start);
int PX(pos_mod)(int dividend, int divisor);
void PX(evaluate_user_block_size)(
    int rnk_pm, const INT *pn, const INT *block, const MPI_Comm *comms_pm,
    INT *block_intern);
void PX(evaluate_user_gcells)(
    int rnk_n, const INT *gc_below, const INT *gc_above,
    INT *gc_below_intern, INT *gc_above_intern);

int PX(flag_active)(
    unsigned flags, unsigned search_flag);
int PX(flag_not_active)(
    unsigned flags, unsigned search_flag);

void PX(physical_dft_size)(
    int rnk_n, const INT *n, unsigned trafo_flag,
    INT *pn);
INT PX(physical_dft_size_1d)(
    INT n, unsigned trafo_flag);

INT* PX(malloc_INT)(
    size_t size);
int* PX(malloc_int)(
    size_t size);
unsigned* PX(malloc_unsigned)(
    size_t size);

int PX(needs_3dto2d_remap)(
    int rnk_n, MPI_Comm comm_cart);


#endif /* !PFFT_UTIL_H */
