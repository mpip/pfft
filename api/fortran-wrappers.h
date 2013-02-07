/*
 * Copyright (c) 2003, 2007-8 Matteo Frigo
 * Copyright (c) 2003, 2007-8 Massachusetts Institute of Technology
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

/* Functions in the PFFT Fortran API, mangled according to the
   FORT(...) macro.  This file is designed to be #included by
   fortran_api.c, possibly multiple times in order to support multiple
   compiler manglings (via redefinition of FORT). */

/* Generate prototypes from fortran-wrappers.h using the vim search/replace
 * command :%s#)\n{[A-Za-z0-9_,\n ();*=&/\[\]-]*}\n#);#gc
 */

/* include prototypes of all fortran wrappers to avoid pedantic gcc warnings  */
#include <fortran-prototypes.h>


PFFT_VOIDFUNC FORT(init, INIT)(void)
{
     PX(init)();
}

PFFT_VOIDFUNC FORT(cleanup, CLEANUP)(void)
{
     PX(cleanup)();
}

PFFT_VOIDFUNC FORT(execute, EXECUTE)(
    PX(plan) * const p)
{
  PX(execute)(*p);
}

PFFT_VOIDFUNC FORT(destroy_plan, DESTROY_PLAN)(
    PX(plan) *p
    )
{
  PX(destroy_plan)(*p);
}

PFFT_VOIDFUNC FORT(init_input_c2c_3d, INIT_INPUT_C2C_3D)(
    const INT *n, const INT *local_n, const INT *local_start,
    C *data
    )
{
  INT n_rev[3], local_n_rev[3], local_start_rev[3];
  
  revert_INT(3, n, n_rev);
  revert_INT(3, local_n, local_n_rev);
  revert_and_sub_ones_INT(3, local_start, local_start_rev);

  PX(init_input_c2c_3d)(n_rev, local_n_rev, local_start_rev, data);
}

PFFT_VOIDFUNC FORT(init_input_c2c, INIT_INPUT_C2C)(
    int *rnk_n, const INT *n, const INT *local_n, const INT *local_start,
    C *data
    )
{
  INT *n_rev = malloc_and_revert_INT(*rnk_n, n);
  INT *local_n_rev = malloc_and_revert_INT(*rnk_n, local_n);
  INT *local_start_rev = PX(malloc_INT)(*rnk_n);
  revert_and_sub_ones_INT(*rnk_n, local_start, local_start_rev);

  PX(init_input_c2c)(*rnk_n, n_rev, local_n_rev, local_start_rev, data);

  free(n_rev); free(local_n_rev); free(local_start_rev);
}


PFFT_VOIDFUNC FORT(init_input_r2c_3d, INIT_INPUT_R2C_3D)(
    const INT *n, const INT *local_n, const INT *local_start,
    R *data
    )
{
  INT n_rev[3], local_n_rev[3], local_start_rev[3];
  
  revert_INT(3, n, n_rev);
  revert_INT(3, local_n, local_n_rev);
  revert_and_sub_ones_INT(3, local_start, local_start_rev);

  PX(init_input_r2c_3d)(n_rev, local_n_rev, local_start_rev, data);
}

PFFT_VOIDFUNC FORT(init_input_r2c, INIT_INPUT_R2C)(
    int *rnk_n, const INT *n, const INT *local_n, const INT *local_start,
    R *data
    )
{
  INT *n_rev = malloc_and_revert_INT(*rnk_n, n);
  INT *local_n_rev = malloc_and_revert_INT(*rnk_n, local_n);
  INT *local_start_rev = PX(malloc_INT)(*rnk_n);
  revert_and_sub_ones_INT(*rnk_n, local_start, local_start_rev);

  PX(init_input_r2c)(*rnk_n, n_rev, local_n_rev, local_start_rev, data);

  free(n_rev); free(local_n_rev); free(local_start_rev);
}


PFFT_VOIDFUNC FORT(init_input_r2r_3d, INIT_INPUT_R2R_3D)(
    const INT *n, const INT *local_n, const INT *local_start,
    R *data
    )
{
  INT n_rev[3], local_n_rev[3], local_start_rev[3];
  
  revert_INT(3, n, n_rev);
  revert_INT(3, local_n, local_n_rev);
  revert_and_sub_ones_INT(3, local_start, local_start_rev);

  PX(init_input_r2r_3d)(n_rev, local_n_rev, local_start_rev, data);
}

PFFT_VOIDFUNC FORT(init_input_r2r, INIT_INPUT_R2R)(
    int *rnk_n, const INT *n, const INT *local_n, const INT *local_start,
    R *data
    )
{
  INT *n_rev = malloc_and_revert_INT(*rnk_n, n);
  INT *local_n_rev = malloc_and_revert_INT(*rnk_n, local_n);
  INT *local_start_rev = PX(malloc_INT)(*rnk_n);
  revert_and_sub_ones_INT(*rnk_n, local_start, local_start_rev);

  PX(init_input_r2r)(*rnk_n, n_rev, local_n_rev, local_start_rev, data);

  free(n_rev); free(local_n_rev); free(local_start_rev);
}



PFFT_VOIDFUNC FORT(check_output_c2c_3d, CHECK_OUTPUT_C2C_3D)(
    R *err, const INT *n, const INT *local_n, const INT *local_start,
    const C *data, MPI_Fint *comm
    )
{
  INT n_rev[3], local_n_rev[3], local_start_rev[3];
  
  revert_INT(3, n, n_rev);
  revert_INT(3, local_n, local_n_rev);
  revert_and_sub_ones_INT(3, local_start, local_start_rev);

  *err = PX(check_output_c2c_3d)(n_rev, local_n_rev, local_start_rev,
      data, MPI_Comm_f2c(*comm));
}

PFFT_VOIDFUNC FORT(check_output_c2c, CHECK_OUTPUT_C2C)(
    R *err, int *rnk_n, const INT *n, const INT *local_n, const INT *local_start,
    const C *data, MPI_Fint *comm
    )
{
  INT *n_rev = malloc_and_revert_INT(*rnk_n, n);
  INT *local_n_rev = malloc_and_revert_INT(*rnk_n, local_n);
  INT *local_start_rev = PX(malloc_INT)(*rnk_n);
  revert_and_sub_ones_INT(*rnk_n, local_start, local_start_rev);

  *err = PX(check_output_c2c)(*rnk_n, n_rev, local_n_rev, local_start_rev,
      data, MPI_Comm_f2c(*comm));

  free(n_rev); free(local_n_rev); free(local_start_rev);
}


PFFT_VOIDFUNC FORT(check_output_c2r_3d, CHECK_OUTPUT_C2R_3D)(
    R *err, const INT *n, const INT *local_n, const INT *local_start,
    const R *data, MPI_Fint *comm
    )
{
  INT n_rev[3], local_n_rev[3], local_start_rev[3];
  
  revert_INT(3, n, n_rev);
  revert_INT(3, local_n, local_n_rev);
  revert_and_sub_ones_INT(3, local_start, local_start_rev);

  *err = PX(check_output_c2r_3d)(n_rev, local_n_rev, local_start_rev,
      data, MPI_Comm_f2c(*comm));
}

PFFT_VOIDFUNC FORT(check_output_c2r, CHECK_OUTPUT_C2R)(
    R *err, int *rnk_n, const INT *n, const INT *local_n, const INT *local_start,
    const R *data, MPI_Fint *comm
    )
{
  INT *n_rev = malloc_and_revert_INT(*rnk_n, n);
  INT *local_n_rev = malloc_and_revert_INT(*rnk_n, local_n);
  INT *local_start_rev = PX(malloc_INT)(*rnk_n);
  revert_and_sub_ones_INT(*rnk_n, local_start, local_start_rev);

  *err = PX(check_output_c2r)(*rnk_n, n_rev, local_n_rev, local_start_rev,
      data, MPI_Comm_f2c(*comm));

  free(n_rev); free(local_n_rev); free(local_start_rev);
}

PFFT_VOIDFUNC FORT(check_output_r2r_3d, CHECK_OUTPUT_R2R_3D)(
    R *err, const INT *n, const INT *local_n, const INT *local_start,
    const R *data, MPI_Fint *comm
    )
{
  INT n_rev[3], local_n_rev[3], local_start_rev[3];
  
  revert_INT(3, n, n_rev);
  revert_INT(3, local_n, local_n_rev);
  revert_and_sub_ones_INT(3, local_start, local_start_rev);

  *err = PX(check_output_r2r_3d)(n_rev, local_n_rev, local_start_rev,
      data, MPI_Comm_f2c(*comm));
}

PFFT_VOIDFUNC FORT(check_output_r2r, CHECK_OUTPUT_R2R)(
    R *err, int *rnk_n, const INT *n, const INT *local_n, const INT *local_start,
    const R *data, MPI_Fint *comm
    )
{
  INT *n_rev = malloc_and_revert_INT(*rnk_n, n);
  INT *local_n_rev = malloc_and_revert_INT(*rnk_n, local_n);
  INT *local_start_rev = PX(malloc_INT)(*rnk_n);
  revert_and_sub_ones_INT(*rnk_n, local_start, local_start_rev);

  *err = PX(check_output_r2r)(*rnk_n, n_rev, local_n_rev, local_start_rev,
      data, MPI_Comm_f2c(*comm));

  free(n_rev); free(local_n_rev); free(local_start_rev);
}

PFFT_VOIDFUNC FORT(local_size_dft_3d, LOCAL_SIZE_DFT_3D)(
    INT *alloc_local,
    const INT *n,
    MPI_Fint *comm_cart, int *pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  INT *n_rev = malloc_and_revert_INT(3, n);
  INT local_ni_rev[3], local_i_start_rev[3];
  INT local_no_rev[3], local_o_start_rev[3];
  
  *alloc_local = PX(local_size_dft_3d)(
      n_rev, MPI_Comm_f2c(*comm_cart), (unsigned) *pfft_flags,
      local_ni_rev, local_i_start_rev, local_no_rev, local_o_start_rev);
  
  revert_INT(3, local_ni_rev, local_ni);
  revert_and_add_ones_INT(3, local_i_start_rev, local_i_start);
  revert_INT(3, local_no_rev, local_no);
  revert_and_add_ones_INT(3, local_o_start_rev, local_o_start);
  
  PX(free)(n_rev);
}

PFFT_VOIDFUNC FORT(local_size_dft_r2c_3d, LOCAL_SIZE_DFT_R2C_3D)(
    INT *alloc_local,
    const INT *n,
    MPI_Fint *comm_cart, int *pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  INT *n_rev = malloc_and_revert_INT(3, n);
  INT local_ni_rev[3], local_i_start_rev[3];
  INT local_no_rev[3], local_o_start_rev[3];
  
  *alloc_local = PX(local_size_dft_r2c_3d)(
      n_rev, MPI_Comm_f2c(*comm_cart), (unsigned) *pfft_flags,
      local_ni_rev, local_i_start_rev, local_no_rev, local_o_start_rev);
  
  revert_INT(3, local_ni_rev, local_ni);
  revert_and_add_ones_INT(3, local_i_start_rev, local_i_start);
  revert_INT(3, local_no_rev, local_no);
  revert_and_add_ones_INT(3, local_o_start_rev, local_o_start);
  
  PX(free)(n_rev);
}

PFFT_VOIDFUNC FORT(local_size_dft_c2r_3d, LOCAL_SIZE_DFT_C2R_3D)(
    INT *alloc_local,
    const INT *n,
    MPI_Fint *comm_cart, int *pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  INT *n_rev = malloc_and_revert_INT(3, n);
  INT local_ni_rev[3], local_i_start_rev[3];
  INT local_no_rev[3], local_o_start_rev[3];
  
  *alloc_local = PX(local_size_dft_c2r_3d)(
      n_rev, MPI_Comm_f2c(*comm_cart), (unsigned) *pfft_flags,
      local_ni_rev, local_i_start_rev, local_no_rev, local_o_start_rev);
  
  revert_INT(3, local_ni_rev, local_ni);
  revert_and_add_ones_INT(3, local_i_start_rev, local_i_start);
  revert_INT(3, local_no_rev, local_no);
  revert_and_add_ones_INT(3, local_o_start_rev, local_o_start);
  
  PX(free)(n_rev);
}

PFFT_VOIDFUNC FORT(local_size_r2r_3d, LOCAL_SIZE_DFT_R2R_3D)(
    INT *alloc_local,
    const INT *n,
    MPI_Fint *comm_cart, int *pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  INT *n_rev = malloc_and_revert_INT(3, n);
  INT local_ni_rev[3], local_i_start_rev[3];
  INT local_no_rev[3], local_o_start_rev[3];
  
  *alloc_local = PX(local_size_r2r_3d)(
      n_rev, MPI_Comm_f2c(*comm_cart), (unsigned) *pfft_flags,
      local_ni_rev, local_i_start_rev, local_no_rev, local_o_start_rev);
  
  revert_INT(3, local_ni_rev, local_ni);
  revert_and_add_ones_INT(3, local_i_start_rev, local_i_start);
  revert_INT(3, local_no_rev, local_no);
  revert_and_add_ones_INT(3, local_o_start_rev, local_o_start);
  
  PX(free)(n_rev);
}


PFFT_VOIDFUNC FORT(local_size_dft, LOCAL_SIZE_DFT)(
    INT *alloc_local,
    int *rnk_n, const INT *n,
    MPI_Fint *comm_cart, int *pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  INT *n_rev = malloc_and_revert_INT(*rnk_n, n);
  INT *local_ni_rev = PX(malloc_INT)(*rnk_n);
  INT *local_i_start_rev = PX(malloc_INT)(*rnk_n);
  INT *local_no_rev = PX(malloc_INT)(*rnk_n);
  INT *local_o_start_rev = PX(malloc_INT)(*rnk_n);
  
  *alloc_local = PX(local_size_dft)(
      *rnk_n, n_rev, MPI_Comm_f2c(*comm_cart), (unsigned) *pfft_flags,
      local_ni_rev, local_i_start_rev, local_no_rev, local_o_start_rev);
  
  revert_INT(*rnk_n, local_ni_rev, local_ni);
  revert_and_add_ones_INT(*rnk_n, local_i_start_rev, local_i_start);
  revert_INT(*rnk_n, local_no_rev, local_no);
  revert_and_add_ones_INT(*rnk_n, local_o_start_rev, local_o_start);
  
  PX(free)(n_rev);
  PX(free)(local_ni_rev); PX(free)(local_i_start_rev);
  PX(free)(local_no_rev); PX(free)(local_o_start_rev);
}

PFFT_VOIDFUNC FORT(local_size_dft_r2c, LOCAL_SIZE_DFT_R2C)(
    INT *alloc_local,
    int *rnk_n, const INT *n,
    MPI_Fint *comm_cart, int *pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  INT *n_rev = malloc_and_revert_INT(*rnk_n, n);
  INT *local_ni_rev = PX(malloc_INT)(*rnk_n);
  INT *local_i_start_rev = PX(malloc_INT)(*rnk_n);
  INT *local_no_rev = PX(malloc_INT)(*rnk_n);
  INT *local_o_start_rev = PX(malloc_INT)(*rnk_n);
  
  *alloc_local = PX(local_size_dft_r2c)(
      *rnk_n, n_rev, MPI_Comm_f2c(*comm_cart), (unsigned) *pfft_flags,
      local_ni_rev, local_i_start_rev, local_no_rev, local_o_start_rev);
  
  revert_INT(*rnk_n, local_ni_rev, local_ni);
  revert_and_add_ones_INT(*rnk_n, local_i_start_rev, local_i_start);
  revert_INT(*rnk_n, local_no_rev, local_no);
  revert_and_add_ones_INT(*rnk_n, local_o_start_rev, local_o_start);
  
  PX(free)(n_rev);
  PX(free)(local_ni_rev); PX(free)(local_i_start_rev);
  PX(free)(local_no_rev); PX(free)(local_o_start_rev);
}

PFFT_VOIDFUNC FORT(local_size_dft_c2r, LOCAL_SIZE_DFT_C2R)(
    INT *alloc_local,
    int *rnk_n, const INT *n,
    MPI_Fint *comm_cart, int *pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  INT *n_rev = malloc_and_revert_INT(*rnk_n, n);
  INT *local_ni_rev = PX(malloc_INT)(*rnk_n);
  INT *local_i_start_rev = PX(malloc_INT)(*rnk_n);
  INT *local_no_rev = PX(malloc_INT)(*rnk_n);
  INT *local_o_start_rev = PX(malloc_INT)(*rnk_n);
  
  *alloc_local = PX(local_size_dft_c2r)(
      *rnk_n, n_rev, MPI_Comm_f2c(*comm_cart), (unsigned) *pfft_flags,
      local_ni_rev, local_i_start_rev, local_no_rev, local_o_start_rev);
  
  revert_INT(*rnk_n, local_ni_rev, local_ni);
  revert_and_add_ones_INT(*rnk_n, local_i_start_rev, local_i_start);
  revert_INT(*rnk_n, local_no_rev, local_no);
  revert_and_add_ones_INT(*rnk_n, local_o_start_rev, local_o_start);
  
  PX(free)(n_rev);
  PX(free)(local_ni_rev); PX(free)(local_i_start_rev);
  PX(free)(local_no_rev); PX(free)(local_o_start_rev);
}

PFFT_VOIDFUNC FORT(local_size_r2r, LOCAL_SIZE_R2R)(
    INT *alloc_local,
    int *rnk_n, const INT *n,
    MPI_Fint *comm_cart, int *pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  INT *n_rev = malloc_and_revert_INT(*rnk_n, n);
  INT *local_ni_rev = PX(malloc_INT)(*rnk_n);
  INT *local_i_start_rev = PX(malloc_INT)(*rnk_n);
  INT *local_no_rev = PX(malloc_INT)(*rnk_n);
  INT *local_o_start_rev = PX(malloc_INT)(*rnk_n);
  
  *alloc_local = PX(local_size_r2r)(
      *rnk_n, n_rev, MPI_Comm_f2c(*comm_cart), (unsigned) *pfft_flags,
      local_ni_rev, local_i_start_rev, local_no_rev, local_o_start_rev);
  
  revert_INT(*rnk_n, local_ni_rev, local_ni);
  revert_and_add_ones_INT(*rnk_n, local_i_start_rev, local_i_start);
  revert_INT(*rnk_n, local_no_rev, local_no);
  revert_and_add_ones_INT(*rnk_n, local_o_start_rev, local_o_start);
  
  PX(free)(n_rev);
  PX(free)(local_ni_rev); PX(free)(local_i_start_rev);
  PX(free)(local_no_rev); PX(free)(local_o_start_rev);
}


PFFT_VOIDFUNC FORT(local_size_many_dft, LOCAL_SIZE_MANY_DFT)(
    INT *alloc_local,
    int *rnk_n, const INT *n, const INT *ni, const INT *no,
    INT *howmany, const INT *iblock, const INT *oblock,
    MPI_Fint *comm_cart, int *pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  int rnk_pm;

  INT *n_rev = malloc_and_revert_INT(*rnk_n, n);
  INT *ni_rev = malloc_and_revert_INT(*rnk_n, ni);
  INT *no_rev = malloc_and_revert_INT(*rnk_n, no);
  INT *local_ni_rev = PX(malloc_INT)(*rnk_n);
  INT *local_i_start_rev = PX(malloc_INT)(*rnk_n);
  INT *local_no_rev = PX(malloc_INT)(*rnk_n);
  INT *local_o_start_rev = PX(malloc_INT)(*rnk_n);
  INT *iblock_rev, *oblock_rev;
  MPI_Cartdim_get(MPI_Comm_f2c(*comm_cart), &rnk_pm);
  iblock_rev = malloc_and_revert_blocks(rnk_pm, iblock);
  oblock_rev = malloc_and_revert_blocks(rnk_pm, oblock);

  *alloc_local = PX(local_size_many_dft)(
      *rnk_n, n_rev, ni_rev, no_rev, *howmany, iblock_rev, oblock_rev,
      MPI_Comm_f2c(*comm_cart), (unsigned) *pfft_flags,
      local_ni_rev, local_i_start_rev, local_no_rev, local_o_start_rev);
  
  revert_INT(*rnk_n, local_ni_rev, local_ni);
  revert_and_add_ones_INT(*rnk_n, local_i_start_rev, local_i_start);
  revert_INT(*rnk_n, local_no_rev, local_no);
  revert_and_add_ones_INT(*rnk_n, local_o_start_rev, local_o_start);

  PX(free)(n_rev); PX(free)(ni_rev); PX(free)(no_rev);
  PX(free)(local_ni_rev); PX(free)(local_i_start_rev);
  PX(free)(local_no_rev); PX(free)(local_o_start_rev);
  free_blocks(&iblock_rev); free_blocks(&oblock_rev);
}

PFFT_VOIDFUNC FORT(local_size_many_dft_r2c, LOCAL_SIZE_MANY_DFT_R2C)(
    INT *alloc_local,
    int *rnk_n, const INT *n, const INT *ni, const INT *no,
    INT *howmany, const INT *iblock, const INT *oblock,
    MPI_Fint *comm_cart, int *pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  int rnk_pm;

  INT *n_rev = malloc_and_revert_INT(*rnk_n, n);
  INT *ni_rev = malloc_and_revert_INT(*rnk_n, ni);
  INT *no_rev = malloc_and_revert_INT(*rnk_n, no);
  INT *local_ni_rev = PX(malloc_INT)(*rnk_n);
  INT *local_i_start_rev = PX(malloc_INT)(*rnk_n);
  INT *local_no_rev = PX(malloc_INT)(*rnk_n);
  INT *local_o_start_rev = PX(malloc_INT)(*rnk_n);
  INT *iblock_rev, *oblock_rev;
  MPI_Cartdim_get(MPI_Comm_f2c(*comm_cart), &rnk_pm);
  iblock_rev = malloc_and_revert_blocks(rnk_pm, iblock);
  oblock_rev = malloc_and_revert_blocks(rnk_pm, oblock);

  *alloc_local = PX(local_size_many_dft_r2c)(
      *rnk_n, n_rev, ni_rev, no_rev, *howmany, iblock_rev, oblock_rev,
      MPI_Comm_f2c(*comm_cart), (unsigned) *pfft_flags,
      local_ni_rev, local_i_start_rev, local_no_rev, local_o_start_rev);
  
  revert_INT(*rnk_n, local_ni_rev, local_ni);
  revert_and_add_ones_INT(*rnk_n, local_i_start_rev, local_i_start);
  revert_INT(*rnk_n, local_no_rev, local_no);
  revert_and_add_ones_INT(*rnk_n, local_o_start_rev, local_o_start);

  PX(free)(n_rev); PX(free)(ni_rev); PX(free)(no_rev);
  PX(free)(local_ni_rev); PX(free)(local_i_start_rev);
  PX(free)(local_no_rev); PX(free)(local_o_start_rev);
  free_blocks(&iblock_rev); free_blocks(&oblock_rev);
}

PFFT_VOIDFUNC FORT(local_size_many_dft_c2r, LOCAL_SIZE_MANY_DFT_C2R)(
    INT *alloc_local,
    int *rnk_n, const INT *n, const INT *ni, const INT *no,
    INT *howmany, const INT *iblock, const INT *oblock,
    MPI_Fint *comm_cart, int *pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  int rnk_pm;

  INT *n_rev = malloc_and_revert_INT(*rnk_n, n);
  INT *ni_rev = malloc_and_revert_INT(*rnk_n, ni);
  INT *no_rev = malloc_and_revert_INT(*rnk_n, no);
  INT *local_ni_rev = PX(malloc_INT)(*rnk_n);
  INT *local_i_start_rev = PX(malloc_INT)(*rnk_n);
  INT *local_no_rev = PX(malloc_INT)(*rnk_n);
  INT *local_o_start_rev = PX(malloc_INT)(*rnk_n);
  INT *iblock_rev, *oblock_rev;
  MPI_Cartdim_get(MPI_Comm_f2c(*comm_cart), &rnk_pm);
  iblock_rev = malloc_and_revert_blocks(rnk_pm, iblock);
  oblock_rev = malloc_and_revert_blocks(rnk_pm, oblock);

  *alloc_local = PX(local_size_many_dft_c2r)(
      *rnk_n, n_rev, ni_rev, no_rev, *howmany, iblock_rev, oblock_rev,
      MPI_Comm_f2c(*comm_cart), (unsigned) *pfft_flags,
      local_ni_rev, local_i_start_rev, local_no_rev, local_o_start_rev);
  
  revert_INT(*rnk_n, local_ni_rev, local_ni);
  revert_and_add_ones_INT(*rnk_n, local_i_start_rev, local_i_start);
  revert_INT(*rnk_n, local_no_rev, local_no);
  revert_and_add_ones_INT(*rnk_n, local_o_start_rev, local_o_start);

  PX(free)(n_rev); PX(free)(ni_rev); PX(free)(no_rev);
  PX(free)(local_ni_rev); PX(free)(local_i_start_rev);
  PX(free)(local_no_rev); PX(free)(local_o_start_rev);
  free_blocks(&iblock_rev); free_blocks(&oblock_rev);
}

PFFT_VOIDFUNC FORT(local_size_many_r2r, LOCAL_SIZE_MANY_R2R)(
    INT *alloc_local,
    int *rnk_n, const INT *n, const INT *ni, const INT *no,
    INT *howmany, const INT *iblock, const INT *oblock,
    MPI_Fint *comm_cart, int *pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    )
{
  int rnk_pm;

  INT *n_rev = malloc_and_revert_INT(*rnk_n, n);
  INT *ni_rev = malloc_and_revert_INT(*rnk_n, ni);
  INT *no_rev = malloc_and_revert_INT(*rnk_n, no);
  INT *local_ni_rev = PX(malloc_INT)(*rnk_n);
  INT *local_i_start_rev = PX(malloc_INT)(*rnk_n);
  INT *local_no_rev = PX(malloc_INT)(*rnk_n);
  INT *local_o_start_rev = PX(malloc_INT)(*rnk_n);
  INT *iblock_rev, *oblock_rev;
  MPI_Cartdim_get(MPI_Comm_f2c(*comm_cart), &rnk_pm);
  iblock_rev = malloc_and_revert_blocks(rnk_pm, iblock);
  oblock_rev = malloc_and_revert_blocks(rnk_pm, oblock);

  *alloc_local = PX(local_size_many_r2r)(
      *rnk_n, n_rev, ni_rev, no_rev, *howmany, iblock_rev, oblock_rev,
      MPI_Comm_f2c(*comm_cart), (unsigned) *pfft_flags,
      local_ni_rev, local_i_start_rev, local_no_rev, local_o_start_rev);
  
  revert_INT(*rnk_n, local_ni_rev, local_ni);
  revert_and_add_ones_INT(*rnk_n, local_i_start_rev, local_i_start);
  revert_INT(*rnk_n, local_no_rev, local_no);
  revert_and_add_ones_INT(*rnk_n, local_o_start_rev, local_o_start);

  PX(free)(n_rev); PX(free)(ni_rev); PX(free)(no_rev);
  PX(free)(local_ni_rev); PX(free)(local_i_start_rev);
  PX(free)(local_no_rev); PX(free)(local_o_start_rev);
  free_blocks(&iblock_rev); free_blocks(&oblock_rev);
}



PFFT_VOIDFUNC FORT(plan_dft_3d, PLAN_DFT_3D)(
    PX(plan) *p,
    const INT *n, C *in, C *out,
    MPI_Fint *comm_cart, int *sign,
    int *pfft_flags
    )
{
  INT *n_rev = malloc_and_revert_INT(3, n);
  
  *p = PX(plan_dft_3d)(n_rev, in, out, MPI_Comm_f2c(*comm_cart),
       *sign, (unsigned) *pfft_flags);
  
  PX(free)(n_rev);
}

PFFT_VOIDFUNC FORT(plan_dft_r2c_3d, PLAN_DFT_R2C_3D)(
    PX(plan) *p,
    const INT *n, R *in, C *out,
    MPI_Fint *comm_cart, int *sign,
    int *pfft_flags
    )
{
  INT *n_rev = malloc_and_revert_INT(3, n);
  
  *p = PX(plan_dft_r2c_3d)(n_rev, in, out, MPI_Comm_f2c(*comm_cart),
       *sign, (unsigned) *pfft_flags);
  
  PX(free)(n_rev);
}

PFFT_VOIDFUNC FORT(plan_dft_c2r_3d, PLAN_DFT_C2R_3D)(
    PX(plan) *p,
    const INT *n, C *in, R *out,
    MPI_Fint *comm_cart, int *sign,
    int *pfft_flags
    )
{
  INT *n_rev = malloc_and_revert_INT(3, n);
  
  *p = PX(plan_dft_c2r_3d)(n_rev, in, out, MPI_Comm_f2c(*comm_cart),
       *sign, (unsigned) *pfft_flags);
  
  PX(free)(n_rev);
}

PFFT_VOIDFUNC FORT(plan_r2r_3d, PLAN_R2R_3D)(
    PX(plan) *p,
    const INT *n, R *in, R *out,
    MPI_Fint *comm_cart, int *kinds,
    int *pfft_flags
    )
{
  INT *n_rev = malloc_and_revert_INT(3, n);
  PX(r2r_kind) *kinds_rev = malloc_and_revert_kinds(3, kinds);

  *p = PX(plan_r2r_3d)(n_rev, in, out, MPI_Comm_f2c(*comm_cart),
      kinds_rev, (unsigned) *pfft_flags);
  
  PX(free)(n_rev);
  PX(free)(kinds_rev);
}


PFFT_VOIDFUNC FORT(plan_dft, PLAN_DFT)(
    PX(plan) *p,
    int *rnk_n, const INT *n, C *in, C *out,
    MPI_Fint *comm_cart, int *sign,
    int *pfft_flags
    )
{
  INT *n_rev = malloc_and_revert_INT(*rnk_n, n);
  
  *p = PX(plan_dft)(*rnk_n, n_rev, in, out, MPI_Comm_f2c(*comm_cart),
       *sign, (unsigned) *pfft_flags);
  
  PX(free)(n_rev);
}

PFFT_VOIDFUNC FORT(plan_dft_r2c, PLAN_DFT_R2C)(
    PX(plan) *p,
    int *rnk_n, const INT *n, R *in, C *out,
    MPI_Fint *comm_cart, int *sign,
    int *pfft_flags
    )
{
  INT *n_rev = malloc_and_revert_INT(*rnk_n, n);
  
  *p = PX(plan_dft_r2c)(*rnk_n, n_rev, in, out, MPI_Comm_f2c(*comm_cart),
       *sign, (unsigned) *pfft_flags);
  
  PX(free)(n_rev);
}

PFFT_VOIDFUNC FORT(plan_dft_c2r, PLAN_DFT_C2R)(
    PX(plan) *p,
    int *rnk_n, const INT *n, C *in, R *out,
    MPI_Fint *comm_cart, int *sign,
    int *pfft_flags
    )
{
  INT *n_rev = malloc_and_revert_INT(*rnk_n, n);
  
  *p = PX(plan_dft_c2r)(*rnk_n, n_rev, in, out, MPI_Comm_f2c(*comm_cart),
       *sign, (unsigned) *pfft_flags);
  
  PX(free)(n_rev);
}

PFFT_VOIDFUNC FORT(plan_r2r, PLAN_R2R)(
    PX(plan) *p,
    int *rnk_n, const INT *n, R *in, R *out,
    MPI_Fint *comm_cart, int *kinds,
    int *pfft_flags
    )
{
  INT *n_rev = malloc_and_revert_INT(*rnk_n, n);
  PX(r2r_kind) *kinds_rev = malloc_and_revert_kinds(*rnk_n, kinds);
  
  *p = PX(plan_r2r)(*rnk_n, n_rev, in, out, MPI_Comm_f2c(*comm_cart),
      kinds_rev, (unsigned) *pfft_flags);
  
  PX(free)(n_rev);
  PX(free)(kinds_rev);
}


PFFT_VOIDFUNC FORT(plan_many_dft, PLAN_MANY_DFT)(
    PX(plan) *p,
    int *rnk_n, const INT *n, const INT *ni, const INT *no,
    INT *howmany, const INT *iblock, const INT *oblock,
    C *in, C *out, MPI_Fint *comm_cart,
    int *sign, int *pfft_flags
    )
{
  int rnk_pm;

  INT *n_rev = malloc_and_revert_INT(*rnk_n, n);
  INT *ni_rev = malloc_and_revert_INT(*rnk_n, ni);
  INT *no_rev = malloc_and_revert_INT(*rnk_n, no);
  INT *iblock_rev, *oblock_rev;
  MPI_Cartdim_get(MPI_Comm_f2c(*comm_cart), &rnk_pm);
  iblock_rev = malloc_and_revert_blocks(rnk_pm, iblock);
  oblock_rev = malloc_and_revert_blocks(rnk_pm, oblock);

  *p = PX(plan_many_dft)(*rnk_n, n_rev, ni_rev, no_rev,
      *howmany, iblock_rev, oblock_rev, in, out,
      MPI_Comm_f2c(*comm_cart), *sign,
      (unsigned) *pfft_flags);
  
  PX(free)(n_rev); PX(free)(ni_rev); PX(free)(no_rev);
  free_blocks(&iblock_rev); free_blocks(&oblock_rev);
}

PFFT_VOIDFUNC FORT(plan_many_dft_r2c, PLAN_MANY_DFT_R2C)(
    PX(plan) *p,
    int *rnk_n, const INT *n, const INT *ni, const INT *no,
    INT *howmany, const INT *iblock, const INT *oblock,
    R *in, C *out, MPI_Fint *comm_cart,
    int *sign, int *pfft_flags
    )
{
  int rnk_pm;

  INT *n_rev = malloc_and_revert_INT(*rnk_n, n);
  INT *ni_rev = malloc_and_revert_INT(*rnk_n, ni);
  INT *no_rev = malloc_and_revert_INT(*rnk_n, no);
  INT *iblock_rev, *oblock_rev;
  MPI_Cartdim_get(MPI_Comm_f2c(*comm_cart), &rnk_pm);
  iblock_rev = malloc_and_revert_blocks(rnk_pm, iblock);
  oblock_rev = malloc_and_revert_blocks(rnk_pm, oblock);

  *p = PX(plan_many_dft_r2c)(*rnk_n, n_rev, ni_rev, no_rev,
      *howmany, iblock_rev, oblock_rev, in, out,
      MPI_Comm_f2c(*comm_cart), *sign,
      (unsigned) *pfft_flags);
  
  PX(free)(n_rev); PX(free)(ni_rev); PX(free)(no_rev);
  free_blocks(&iblock_rev); free_blocks(&oblock_rev);
}

PFFT_VOIDFUNC FORT(plan_many_dft_c2r, PLAN_MANY_DFT_C2R)(
    PX(plan) *p,
    int *rnk_n, const INT *n, const INT *ni, const INT *no,
    INT *howmany, const INT *iblock, const INT *oblock,
    C *in, R *out, MPI_Fint *comm_cart,
    int *sign, int *pfft_flags
    )
{
  int rnk_pm;

  INT *n_rev = malloc_and_revert_INT(*rnk_n, n);
  INT *ni_rev = malloc_and_revert_INT(*rnk_n, ni);
  INT *no_rev = malloc_and_revert_INT(*rnk_n, no);
  INT *iblock_rev, *oblock_rev;
  MPI_Cartdim_get(MPI_Comm_f2c(*comm_cart), &rnk_pm);
  iblock_rev = malloc_and_revert_blocks(rnk_pm, iblock);
  oblock_rev = malloc_and_revert_blocks(rnk_pm, oblock);

  *p = PX(plan_many_dft_c2r)(*rnk_n, n_rev, ni_rev, no_rev,
      *howmany, iblock_rev, oblock_rev, in, out,
      MPI_Comm_f2c(*comm_cart), *sign,
      (unsigned) *pfft_flags);
  
  PX(free)(n_rev); PX(free)(ni_rev); PX(free)(no_rev);
  free_blocks(&iblock_rev); free_blocks(&oblock_rev);
}

PFFT_VOIDFUNC FORT(plan_many_r2r, PLAN_MANY_R2R)(
    PX(plan) *p,
    int *rnk_n, const INT *n, const INT *ni, const INT *no,
    INT *howmany, const INT *iblock, const INT *oblock,
    R *in, R *out, MPI_Fint *comm_cart,
    int *kinds, int *pfft_flags
    )
{
  int rnk_pm;

  INT *n_rev = malloc_and_revert_INT(*rnk_n, n);
  INT *ni_rev = malloc_and_revert_INT(*rnk_n, ni);
  INT *no_rev = malloc_and_revert_INT(*rnk_n, no);
  INT *iblock_rev, *oblock_rev;
  MPI_Cartdim_get(MPI_Comm_f2c(*comm_cart), &rnk_pm);
  iblock_rev = malloc_and_revert_blocks(rnk_pm, iblock);
  oblock_rev = malloc_and_revert_blocks(rnk_pm, oblock);
  PX(r2r_kind) *kinds_rev = malloc_and_revert_kinds(*rnk_n, kinds);
  
  *p = PX(plan_many_r2r)(*rnk_n, n_rev, ni_rev, no_rev,
      *howmany, iblock_rev, oblock_rev, in, out,
      MPI_Comm_f2c(*comm_cart), kinds_rev,
      (unsigned) *pfft_flags);
  
  PX(free)(n_rev); PX(free)(ni_rev); PX(free)(no_rev);
  free_blocks(&iblock_rev); free_blocks(&oblock_rev);
  PX(free)(kinds_rev);
}

PFFT_VOIDFUNC FORT(plan_many_dft_skipped, PLAN_MANY_DFT_SKIPPED)(
    PX(plan) *p,
    int *rnk_n, const INT *n, const INT *ni, const INT *no,
    INT *howmany, const INT *iblock, const INT *oblock,
    const int *skip_trafos, C *in, C *out, MPI_Fint *comm_cart,
    int *sign, int *pfft_flags
    )
{
  int rnk_pm;

  int *skip_trafos_rev;
  INT *n_rev = malloc_and_revert_INT(*rnk_n, n);
  INT *ni_rev = malloc_and_revert_INT(*rnk_n, ni);
  INT *no_rev = malloc_and_revert_INT(*rnk_n, no);
  INT *iblock_rev, *oblock_rev;
  MPI_Cartdim_get(MPI_Comm_f2c(*comm_cart), &rnk_pm);
  skip_trafos_rev = malloc_and_revert_int(rnk_pm, skip_trafos);
  iblock_rev = malloc_and_revert_blocks(rnk_pm, iblock);
  oblock_rev = malloc_and_revert_blocks(rnk_pm, oblock);

  *p = PX(plan_many_dft_skipped)(*rnk_n, n_rev, ni_rev, no_rev,
      *howmany, iblock_rev, oblock_rev, skip_trafos_rev, in, out,
      MPI_Comm_f2c(*comm_cart), *sign,
      (unsigned) *pfft_flags);
  
  PX(free)(n_rev); PX(free)(ni_rev); PX(free)(no_rev);
  PX(free)(skip_trafos_rev); 
  free_blocks(&iblock_rev); free_blocks(&oblock_rev);
}

PFFT_VOIDFUNC FORT(plan_rgc_3d, PLAN_RGC_3D)(
    PX(gcplan) *p,
    const INT *n, const INT *gc_below, const INT *gc_above,
    R *data, MPI_Fint *comm_cart, int *gc_flags
    )
{
  INT *n_rev = malloc_and_revert_INT(3, n);
  INT *gc_below_rev = malloc_and_revert_gcells(3, gc_below);
  INT *gc_above_rev = malloc_and_revert_gcells(3, gc_above);

  *p = PX(plan_rgc_3d)(n_rev, gc_below_rev, gc_above_rev,
      data, MPI_Comm_f2c(*comm_cart), (unsigned) *gc_flags);

  PX(free)(n_rev);
  free_gcells(&gc_below_rev); free_gcells(&gc_above_rev);
}

PFFT_VOIDFUNC FORT(local_size_gc_3d, LOCAL_SIZE_GC_3D)(
    INT *alloc_local_gc,
    const INT *local_n, const INT *local_n_start, INT *alloc_local,
    const INT *gc_below, const INT *gc_above,
    INT *local_ngc, INT *local_gc_start
    )
{
  INT *local_n_rev = malloc_and_revert_INT(3, local_n);
  INT *local_n_start_rev = malloc_and_revert_INT(3, local_n_start);
  INT *gc_below_rev = malloc_and_revert_gcells(3, gc_below);
  INT *gc_above_rev = malloc_and_revert_gcells(3, gc_above);
  INT *local_ngc_rev = PX(malloc_INT)(3);
  INT *local_gc_start_rev = PX(malloc_INT)(3);

  *alloc_local_gc = PX(local_size_gc_3d)(local_n_rev, local_n_start_rev, *alloc_local,
      gc_below_rev, gc_above_rev,
      local_ngc_rev, local_gc_start_rev);

  revert_INT(3, local_ngc_rev, local_ngc);
  revert_and_add_ones_INT(3, local_gc_start_rev, local_gc_start);

  PX(free)(local_n_rev); PX(free)(local_n_start_rev);
  PX(free)(local_ngc_rev); PX(free)(local_gc_start_rev);
  free_gcells(&gc_below_rev); free_gcells(&gc_above_rev);
}

PFFT_VOIDFUNC FORT(plan_cgc_3d, PLAN_CGC_3D)(
    PX(gcplan) *p,
    const INT *n, const INT *gc_below, const INT *gc_above,
    C *data, MPI_Fint *comm_cart, int *gc_flags
    )
{
  INT *n_rev = malloc_and_revert_INT(3, n);
  INT *gc_below_rev = malloc_and_revert_gcells(3, gc_below);
  INT *gc_above_rev = malloc_and_revert_gcells(3, gc_above);

  *p = PX(plan_cgc_3d)(n_rev, gc_below_rev, gc_above_rev,
      data, MPI_Comm_f2c(*comm_cart), (unsigned) *gc_flags);

  PX(free)(n_rev);
  free_gcells(&gc_below_rev); free_gcells(&gc_above_rev);
}

PFFT_VOIDFUNC FORT(local_size_gc, LOCAL_SIZE_GC)(
    INT *alloc_local_gc,
    int *rnk_n, const INT *local_n, const INT *local_n_start, INT *alloc_local,
    const INT *gc_below, const INT *gc_above,
    INT *local_ngc, INT *local_gc_start
    )
{
  INT *local_n_rev = malloc_and_revert_INT(*rnk_n, local_n);
  INT *local_n_start_rev = malloc_and_revert_INT(*rnk_n, local_n_start);
  INT *gc_below_rev = malloc_and_revert_gcells(*rnk_n, gc_below);
  INT *gc_above_rev = malloc_and_revert_gcells(*rnk_n, gc_above);
  INT *local_ngc_rev = PX(malloc_INT)(*rnk_n);
  INT *local_gc_start_rev = PX(malloc_INT)(*rnk_n);

  *alloc_local_gc = PX(local_size_gc)(*rnk_n, local_n_rev, local_n_start_rev, *alloc_local,
      gc_below_rev, gc_above_rev,
      local_ngc_rev, local_gc_start_rev);

  revert_INT(*rnk_n, local_ngc_rev, local_ngc);
  revert_and_add_ones_INT(*rnk_n, local_gc_start_rev, local_gc_start);

  PX(free)(local_n_rev); PX(free)(local_n_start_rev);
  PX(free)(local_ngc_rev); PX(free)(local_gc_start_rev);
  free_gcells(&gc_below_rev); free_gcells(&gc_above_rev);
}

PFFT_VOIDFUNC FORT(plan_rgc, PLAN_RGC)(
    PX(gcplan) *p,
    int *rnk_n, const INT *n, const INT *gc_below, const INT *gc_above,
    R *data, MPI_Fint *comm_cart, int *gc_flags
    )
{
  INT *n_rev = malloc_and_revert_INT(*rnk_n, n);
  INT *gc_below_rev = malloc_and_revert_gcells(*rnk_n, gc_below);
  INT *gc_above_rev = malloc_and_revert_gcells(*rnk_n, gc_above);

  *p = PX(plan_rgc)(*rnk_n, n_rev, gc_below_rev, gc_above_rev,
      data, MPI_Comm_f2c(*comm_cart), (unsigned) *gc_flags);

  PX(free)(n_rev);
  free_gcells(&gc_below_rev); free_gcells(&gc_above_rev);
}

PFFT_VOIDFUNC FORT(plan_cgc, PLAN_CGC)(
    PX(gcplan) *p,
    int *rnk_n, const INT *n, const INT *gc_below, const INT *gc_above,
    C *data, MPI_Fint *comm_cart, int *gc_flags
    )
{
  INT *n_rev = malloc_and_revert_INT(*rnk_n, n);
  INT *gc_below_rev = malloc_and_revert_gcells(*rnk_n, gc_below);
  INT *gc_above_rev = malloc_and_revert_gcells(*rnk_n, gc_above);

  *p = PX(plan_cgc)(*rnk_n, n_rev, gc_below_rev, gc_above_rev,
      data, MPI_Comm_f2c(*comm_cart), (unsigned) *gc_flags);

  PX(free)(n_rev);
  free_gcells(&gc_below_rev); free_gcells(&gc_above_rev);
}

PFFT_VOIDFUNC FORT(exchange, EXCHANGE)(
    PX(gcplan) *p
    )
{
  PX(exchange)(*p);
}

PFFT_VOIDFUNC FORT(reduce, REDUCE)(
    PX(gcplan) *p
    )
{
  PX(reduce)(*p);
}

PFFT_VOIDFUNC FORT(destroy_gcplan, DESTROY_GCPLAN)(
    PX(gcplan) *p
    )
{
  PX(destroy_gcplan)(*p);
}

PFFT_VOIDFUNC FORT(local_size_many_gc, LOCAL_SIZE_MANY_GC)(
    INT *alloc_local_gc,
    int *rnk_n, const INT *local_n, const INT *local_n_start, INT *alloc_local,
    INT *howmany, const INT *gc_below, const INT *gc_above,
    INT *local_ngc, INT *local_gc_start
    )
{
  INT *local_n_rev = malloc_and_revert_INT(*rnk_n, local_n);
  INT *local_n_start_rev = malloc_and_revert_INT(*rnk_n, local_n_start);
  INT *gc_below_rev = malloc_and_revert_gcells(*rnk_n, gc_below);
  INT *gc_above_rev = malloc_and_revert_gcells(*rnk_n, gc_above);
  INT *local_ngc_rev = PX(malloc_INT)(*rnk_n);
  INT *local_gc_start_rev = PX(malloc_INT)(*rnk_n);

  *alloc_local_gc = PX(local_size_many_gc)(*rnk_n, local_n_rev, local_n_start_rev, *alloc_local,
      *howmany, gc_below_rev, gc_above_rev,
      local_ngc_rev, local_gc_start_rev);

  revert_INT(*rnk_n, local_ngc_rev, local_ngc);
  revert_and_add_ones_INT(*rnk_n, local_gc_start_rev, local_gc_start);

  PX(free)(local_n_rev); PX(free)(local_n_start_rev);
  PX(free)(local_ngc_rev); PX(free)(local_gc_start_rev);
  free_gcells(&gc_below_rev); free_gcells(&gc_above_rev);
}

PFFT_VOIDFUNC FORT(plan_many_rgc, PLAN_MANY_RGC)(
    PX(gcplan) *p,
    int *rnk_n, const INT *n,
    INT *howmany, const INT *block,
    const INT *gc_below, const INT *gc_above,
    R *data, MPI_Fint *comm_cart, int *gc_flags
    )
{
  int rnk_pm;

  INT *n_rev = malloc_and_revert_INT(*rnk_n, n);
  INT *gc_below_rev = malloc_and_revert_gcells(*rnk_n, gc_below);
  INT *gc_above_rev = malloc_and_revert_gcells(*rnk_n, gc_above);
  INT *block_rev;
  MPI_Cartdim_get(MPI_Comm_f2c(*comm_cart), &rnk_pm);
  block_rev = malloc_and_revert_blocks(rnk_pm, block);

  *p = PX(plan_many_rgc)(*rnk_n, n_rev, *howmany, block_rev, gc_below_rev, gc_above_rev,
      data, MPI_Comm_f2c(*comm_cart), (unsigned) *gc_flags);

  PX(free)(n_rev);
  free_blocks(&block_rev);
  free_gcells(&gc_below_rev); free_gcells(&gc_above_rev);
}

PFFT_VOIDFUNC FORT(plan_many_cgc, PLAN_MANY_CGC)(
    PX(gcplan) *p,
    int *rnk_n, const INT *n,
    INT *howmany, const INT *block,
    const INT *gc_below, const INT *gc_above,
    C *data, MPI_Fint *comm_cart, int *gc_flags
    )
{
  int rnk_pm;

  INT *n_rev = malloc_and_revert_INT(*rnk_n, n);
  INT *gc_below_rev = malloc_and_revert_gcells(*rnk_n, gc_below);
  INT *gc_above_rev = malloc_and_revert_gcells(*rnk_n, gc_above);
  INT *block_rev;
  MPI_Cartdim_get(MPI_Comm_f2c(*comm_cart), &rnk_pm);
  block_rev = malloc_and_revert_blocks(rnk_pm, block);

  *p = PX(plan_many_cgc)(*rnk_n, n_rev, *howmany, block_rev, gc_below_rev, gc_above_rev,
      data, MPI_Comm_f2c(*comm_cart), (unsigned) *gc_flags);

  PX(free)(n_rev);
  free_blocks(&block_rev);
  free_gcells(&gc_below_rev); free_gcells(&gc_above_rev);
}





PFFT_VOIDFUNC FORT(prod_INT, PROD_INT)(
    INT *prod, int *d, const INT *vec
    )
{
  *prod = PX(prod_INT)(*d, vec);
}

PFFT_VOIDFUNC FORT(sum_INT, SUM_INT)(
    INT *sum, int *d, const INT *vec
    )
{
  *sum = PX(sum_INT)(*d, vec);
}

PFFT_VOIDFUNC FORT(equal_INT, EQUAL_INT)(
    INT *equal, int *d, const INT *vec1, const INT *vec2
    )
{
  *equal = PX(equal_INT)(*d, vec1, vec2);
}

PFFT_VOIDFUNC FORT(vcopy_INT, VCOPY_INT)(
    int *d, const INT *vec1,
    INT *vec2
    )
{
  PX(vcopy_INT)(*d, vec1,
      vec2);
}

PFFT_VOIDFUNC FORT(vadd_INT, VADD_INT)(
    int *d, const INT *vec1, const INT *vec2,
    INT *sum
    )
{
  PX(vadd_INT)(*d, vec1, vec2,
      sum);
}

PFFT_VOIDFUNC FORT(vsub_INT, VSUB_INT)(
    int *d, const INT *vec1, const INT *vec2,
    INT *sub
    )
{
  PX(vsub_INT)(*d, vec1, vec2,
      sub);
}

PFFT_VOIDFUNC  FORT(apr_complex_3d, APR_COMPLEX_3D)(
    const C *data, const INT *local_n, const INT *local_start,
    const char *name, MPI_Fint *comm
    )
{
  PX(apr_complex_permuted_3d)(data, local_n, local_start,
      2, 1, 0, name, MPI_Comm_f2c(*comm));
}


PFFT_VOIDFUNC  FORT(apr_complex_permuted_3d, APR_COMPLEX_PERMUTED_3D)(
    const C *data, const INT *local_n, const INT *local_start,
    int *perm1, int *perm2, int *perm3, const char *name, MPI_Fint *comm
    )
{
  PX(apr_complex_permuted_3d)(data, local_n, local_start,
      *perm3 - 1, *perm2 - 1, *perm1 - 1, name, MPI_Comm_f2c(*comm));
}


PFFT_VOIDFUNC FORT(message, MESSAGE)(
    MPI_Fint *comm, const char *msg
    )
{
  PX(printf)(MPI_Comm_f2c(*comm), msg);
}


PFFT_VOIDFUNC FORT(error_message, ERROR_MESSAGE)(
    MPI_Fint *comm, const char *msg
    )
{
  PX(fprintf)(MPI_Comm_f2c(*comm), stderr, msg);
}




PFFT_VOIDFUNC FORT(create_procmesh, CREATE_PROCMESH)(
    int *ierror,
    int *rnk, MPI_Fint *comm, int *np,
    MPI_Fint *comm_cart
    )
{
  MPI_Comm comm_cart_clike;
  int *np_rev = malloc_and_revert_int(*rnk, np);

  *ierror = PX(create_procmesh)(*rnk, MPI_Comm_f2c(*comm), np_rev, &comm_cart_clike);
  
  /* C to Fortran type conversion */
  *comm_cart = MPI_Comm_c2f(comm_cart_clike);

  free(np_rev);
}

PFFT_VOIDFUNC FORT(create_procmesh_2d, CREATE_PROCMESH_2D)(
    int *ierror,
    MPI_Fint *comm, int *np2, int *np3,
    MPI_Fint *comm_cart_2d
    )
{
  MPI_Comm comm_cart_2d_clike;

  *ierror = PX(create_procmesh_2d)(MPI_Comm_f2c(*comm), *np3, *np2, &comm_cart_2d_clike);
  
  /* C to Fortran type conversion */
  *comm_cart_2d = MPI_Comm_c2f(comm_cart_2d_clike);
}


