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


PFFT_VOIDFUNC FORT(init, INIT)(void);
PFFT_VOIDFUNC FORT(cleanup, CLEANUP)(void);
PFFT_VOIDFUNC FORT(execute, EXECUTE)(
    PX(plan) * const p);
PFFT_VOIDFUNC FORT(destroy_plan, DESTROY_PLAN)(
    PX(plan) *p
    );
PFFT_VOIDFUNC FORT(init_input_c2c_3d, INIT_INPUT_C2C_3D)(
    const INT *n, const INT *local_n, const INT *local_start,
    C *data
    );
PFFT_VOIDFUNC FORT(init_input_c2c, INIT_INPUT_C2C)(
    int *rnk_n, const INT *n, const INT *local_n, const INT *local_start,
    C *data
    );

PFFT_VOIDFUNC FORT(init_input_r2c_3d, INIT_INPUT_R2C_3D)(
    const INT *n, const INT *local_n, const INT *local_start,
    R *data
    );
PFFT_VOIDFUNC FORT(init_input_r2c, INIT_INPUT_R2C)(
    int *rnk_n, const INT *n, const INT *local_n, const INT *local_start,
    R *data
    );

PFFT_VOIDFUNC FORT(init_input_r2r_3d, INIT_INPUT_R2R_3D)(
    const INT *n, const INT *local_n, const INT *local_start,
    R *data
    );
PFFT_VOIDFUNC FORT(init_input_r2r, INIT_INPUT_R2R)(
    int *rnk_n, const INT *n, const INT *local_n, const INT *local_start,
    R *data
    );


PFFT_VOIDFUNC FORT(check_output_c2c_3d, CHECK_OUTPUT_C2C_3D)(
    R *err, const INT *n, const INT *local_n, const INT *local_start,
    const C *data, MPI_Fint *comm
    );
PFFT_VOIDFUNC FORT(check_output_c2c, CHECK_OUTPUT_C2C)(
    R *err, int *rnk_n, const INT *n, const INT *local_n, const INT *local_start,
    const C *data, MPI_Fint *comm
    );

PFFT_VOIDFUNC FORT(check_output_c2r_3d, CHECK_OUTPUT_C2R_3D)(
    R *err, const INT *n, const INT *local_n, const INT *local_start,
    const R *data, MPI_Fint *comm
    );
PFFT_VOIDFUNC FORT(check_output_c2r, CHECK_OUTPUT_C2R)(
    R *err, int *rnk_n, const INT *n, const INT *local_n, const INT *local_start,
    const R *data, MPI_Fint *comm
    );
PFFT_VOIDFUNC FORT(check_output_r2r_3d, CHECK_OUTPUT_R2R_3D)(
    R *err, const INT *n, const INT *local_n, const INT *local_start,
    const R *data, MPI_Fint *comm
    );
PFFT_VOIDFUNC FORT(check_output_r2r, CHECK_OUTPUT_R2R)(
    R *err, int *rnk_n, const INT *n, const INT *local_n, const INT *local_start,
    const R *data, MPI_Fint *comm
    );
PFFT_VOIDFUNC FORT(local_size_dft_3d, LOCAL_SIZE_DFT_3D)(
    INT *alloc_local,
    const INT *n,
    MPI_Fint *comm_cart, int *pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    );
PFFT_VOIDFUNC FORT(local_size_dft_r2c_3d, LOCAL_SIZE_DFT_R2C_3D)(
    INT *alloc_local,
    const INT *n,
    MPI_Fint *comm_cart, int *pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    );
PFFT_VOIDFUNC FORT(local_size_dft_c2r_3d, LOCAL_SIZE_DFT_C2R_3D)(
    INT *alloc_local,
    const INT *n,
    MPI_Fint *comm_cart, int *pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    );
PFFT_VOIDFUNC FORT(local_size_r2r_3d, LOCAL_SIZE_DFT_R2R_3D)(
    INT *alloc_local,
    const INT *n,
    MPI_Fint *comm_cart, int *pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    );

PFFT_VOIDFUNC FORT(local_size_dft, LOCAL_SIZE_DFT)(
    INT *alloc_local,
    int *rnk_n, const INT *n,
    MPI_Fint *comm_cart, int *pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    );
PFFT_VOIDFUNC FORT(local_size_dft_r2c, LOCAL_SIZE_DFT_R2C)(
    INT *alloc_local,
    int *rnk_n, const INT *n,
    MPI_Fint *comm_cart, int *pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    );
PFFT_VOIDFUNC FORT(local_size_dft_c2r, LOCAL_SIZE_DFT_C2R)(
    INT *alloc_local,
    int *rnk_n, const INT *n,
    MPI_Fint *comm_cart, int *pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    );
PFFT_VOIDFUNC FORT(local_size_r2r, LOCAL_SIZE_R2R)(
    INT *alloc_local,
    int *rnk_n, const INT *n,
    MPI_Fint *comm_cart, int *pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    );

PFFT_VOIDFUNC FORT(local_size_many_dft, LOCAL_SIZE_MANY_DFT)(
    INT *alloc_local,
    int *rnk_n, const INT *n, const INT *ni, const INT *no,
    INT *howmany, const INT *iblock, const INT *oblock,
    MPI_Fint *comm_cart, int *pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    );
PFFT_VOIDFUNC FORT(local_size_many_dft_r2c, LOCAL_SIZE_MANY_DFT_R2C)(
    INT *alloc_local,
    int *rnk_n, const INT *n, const INT *ni, const INT *no,
    INT *howmany, const INT *iblock, const INT *oblock,
    MPI_Fint *comm_cart, int *pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    );
PFFT_VOIDFUNC FORT(local_size_many_dft_c2r, LOCAL_SIZE_MANY_DFT_C2R)(
    INT *alloc_local,
    int *rnk_n, const INT *n, const INT *ni, const INT *no,
    INT *howmany, const INT *iblock, const INT *oblock,
    MPI_Fint *comm_cart, int *pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    );
PFFT_VOIDFUNC FORT(local_size_many_r2r, LOCAL_SIZE_MANY_R2R)(
    INT *alloc_local,
    int *rnk_n, const INT *n, const INT *ni, const INT *no,
    INT *howmany, const INT *iblock, const INT *oblock,
    MPI_Fint *comm_cart, int *pfft_flags,
    INT *local_ni, INT *local_i_start,
    INT *local_no, INT *local_o_start
    );


PFFT_VOIDFUNC FORT(plan_dft_3d, PLAN_DFT_3D)(
    PX(plan) *p,
    const INT *n, C *in, C *out,
    MPI_Fint *comm_cart, int *sign,
    int *pfft_flags
    );
PFFT_VOIDFUNC FORT(plan_dft_r2c_3d, PLAN_DFT_R2C_3D)(
    PX(plan) *p,
    const INT *n, R *in, C *out,
    MPI_Fint *comm_cart, int *sign,
    int *pfft_flags
    );
PFFT_VOIDFUNC FORT(plan_dft_c2r_3d, PLAN_DFT_C2R_3D)(
    PX(plan) *p,
    const INT *n, C *in, R *out,
    MPI_Fint *comm_cart, int *sign,
    int *pfft_flags
    );
PFFT_VOIDFUNC FORT(plan_r2r_3d, PLAN_R2R_3D)(
    PX(plan) *p,
    const INT *n, R *in, R *out,
    MPI_Fint *comm_cart, int *kinds,
    int *pfft_flags
    );

PFFT_VOIDFUNC FORT(plan_dft, PLAN_DFT)(
    PX(plan) *p,
    int *rnk_n, const INT *n, C *in, C *out,
    MPI_Fint *comm_cart, int *sign,
    int *pfft_flags
    );
PFFT_VOIDFUNC FORT(plan_dft_r2c, PLAN_DFT_R2C)(
    PX(plan) *p,
    int *rnk_n, const INT *n, R *in, C *out,
    MPI_Fint *comm_cart, int *sign,
    int *pfft_flags
    );
PFFT_VOIDFUNC FORT(plan_dft_c2r, PLAN_DFT_C2R)(
    PX(plan) *p,
    int *rnk_n, const INT *n, C *in, R *out,
    MPI_Fint *comm_cart, int *sign,
    int *pfft_flags
    );
PFFT_VOIDFUNC FORT(plan_r2r, PLAN_R2R)(
    PX(plan) *p,
    int *rnk_n, const INT *n, R *in, R *out,
    MPI_Fint *comm_cart, int *kinds,
    int *pfft_flags
    );

PFFT_VOIDFUNC FORT(plan_many_dft, PLAN_MANY_DFT)(
    PX(plan) *p,
    int *rnk, const INT *n, const INT *ni, const INT *no,
    INT *howmany, const INT *iblock, const INT *oblock,
    C *in, C *out, MPI_Fint *comm_cart,
    int *sign, int *pfft_flags
    );
PFFT_VOIDFUNC FORT(plan_many_dft_r2c, PLAN_MANY_DFT_R2C)(
    PX(plan) *p,
    int *rnk, const INT *n, const INT *ni, const INT *no,
    INT *howmany, const INT *iblock, const INT *oblock,
    R *in, C *out, MPI_Fint *comm_cart,
    int *sign, int *pfft_flags
    );
PFFT_VOIDFUNC FORT(plan_many_dft_c2r, PLAN_MANY_DFT_C2R)(
    PX(plan) *p,
    int *rnk_n, const INT *n, const INT *ni, const INT *no,
    INT *howmany, const INT *iblock, const INT *oblock,
    C *in, R *out, MPI_Fint *comm_cart,
    int *sign, int *pfft_flags
    );
PFFT_VOIDFUNC FORT(plan_many_r2r, PLAN_MANY_R2R)(
    PX(plan) *p,
    int *rnk_n, const INT *n, const INT *ni, const INT *no,
    INT *howmany, const INT *iblock, const INT *oblock,
    R *in, R *out, MPI_Fint *comm_cart,
    int *kinds, int *pfft_flags
    );
PFFT_VOIDFUNC FORT(plan_many_dft_skipped, PLAN_MANY_DFT_SKIPPED)(
    PX(plan) *p,
    int *rnk, const INT *n, const INT *ni, const INT *no,
    INT *howmany, const INT *iblock, const INT *oblock,
    const int *skip_trafos, C *in, C *out, MPI_Fint *comm_cart,
    int *sign, int *pfft_flags
    );
PFFT_VOIDFUNC FORT(plan_rgc_3d, PLAN_RGC_3D)(
    PX(gcplan) *p,
    const INT *n, const INT *gc_below, const INT *gc_above,
    R *data, MPI_Fint *comm_cart, int *gc_flags
    );
PFFT_VOIDFUNC FORT(local_size_gc_3d, LOCAL_SIZE_GC_3D)(
    INT *alloc_local_gc,
    const INT *local_n, const INT *local_n_start, INT *alloc_local,
    const INT *gc_below, const INT *gc_above,
    INT *local_ngc, INT *local_gc_start
    );
PFFT_VOIDFUNC FORT(plan_cgc_3d, PLAN_CGC_3D)(
    PX(gcplan) *p,
    const INT *n, const INT *gc_below, const INT *gc_above,
    C *data, MPI_Fint *comm_cart, int *gc_flags
    );
PFFT_VOIDFUNC FORT(local_size_gc, LOCAL_SIZE_GC)(
    INT *alloc_local_gc,
    int *rnk, const INT *local_n, const INT *local_n_start, INT *alloc_local,
    const INT *gc_below, const INT *gc_above,
    INT *local_ngc, INT *local_gc_start
    );
PFFT_VOIDFUNC FORT(plan_rgc, PLAN_RGC)(
    PX(gcplan) *p,
    int *rnk, const INT *n, const INT *gc_below, const INT *gc_above,
    R *data, MPI_Fint *comm_cart, int *gc_flags
    );
PFFT_VOIDFUNC FORT(plan_cgc, PLAN_CGC)(
    PX(gcplan) *p,
    int *rnk, const INT *n, const INT *gc_below, const INT *gc_above,
    C *data, MPI_Fint *comm_cart, int *gc_flags
    );
PFFT_VOIDFUNC FORT(exchange, EXCHANGE)(
    PX(gcplan) *p
    );
PFFT_VOIDFUNC FORT(reduce, REDUCE)(
    PX(gcplan) *p
    );
PFFT_VOIDFUNC FORT(destroy_gcplan, DESTROY_GCPLAN)(
    PX(gcplan) *p
    );
PFFT_VOIDFUNC FORT(local_size_many_gc, LOCAL_SIZE_MANY_GC)(
    INT *alloc_local_gc,
    int *rnk, const INT *local_n, const INT *local_n_start, INT *alloc_local,
    INT *howmany, const INT *gc_below, const INT *gc_above,
    INT *local_ngc, INT *local_gc_start
    );
PFFT_VOIDFUNC FORT(plan_many_rgc, PLAN_MANY_RGC)(
    PX(gcplan) *p,
    int *rnk_n, const INT *n,
    INT *howmany, const INT *block,
    const INT *gc_below, const INT *gc_above,
    R *data, MPI_Fint *comm_cart, int *gc_flags
    );
PFFT_VOIDFUNC FORT(plan_many_cgc, PLAN_MANY_CGC)(
    PX(gcplan) *p,
    int *rnk_n, const INT *n,
    INT *howmany, const INT *block,
    const INT *gc_below, const INT *gc_above,
    C *data, MPI_Fint *comm_cart, int *gc_flags
    );




PFFT_VOIDFUNC FORT(prod_INT, PROD_INT)(
    INT *prod, int *d, const INT *vec
    );
PFFT_VOIDFUNC FORT(sum_INT, SUM_INT)(
    INT *sum, int *d, const INT *vec
    );
PFFT_VOIDFUNC FORT(equal_INT, EQUAL_INT)(
    INT *equal, int *d, const INT *vec1, const INT *vec2
    );
PFFT_VOIDFUNC FORT(vcopy_INT, VCOPY_INT)(
    int *d, const INT *vec1,
    INT *vec2
    );
PFFT_VOIDFUNC FORT(vadd_INT, VADD_INT)(
    int *d, const INT *vec1, const INT *vec2,
    INT *sum
    );
PFFT_VOIDFUNC FORT(vsub_INT, VSUB_INT)(
    int *d, const INT *vec1, const INT *vec2,
    INT *sub
    );
PFFT_VOIDFUNC  FORT(apr_complex_3d, APR_COMPLEX_3D)(
    const C *data, const INT *local_n, const INT *local_start,
    const char *name, MPI_Fint *comm
    );

PFFT_VOIDFUNC  FORT(apr_complex_permuted_3d, APR_COMPLEX_PERMUTED_3D)(
    const C *data, const INT *local_n, const INT *local_start,
    int *perm1, int *perm2, int *perm3, const char *name, MPI_Fint *comm
    );

PFFT_VOIDFUNC FORT(message, MESSAGE)(
    MPI_Fint *comm, const char *msg
    );

PFFT_VOIDFUNC FORT(error_message, ERROR_MESSAGE)(
    MPI_Fint *comm, const char *msg
    );



PFFT_VOIDFUNC FORT(create_procmesh, CREATE_PROCMESH)(
    int *ierror,
    int *rnk, MPI_Fint *comm, int *np,
    MPI_Fint *comm_cart
    );
PFFT_VOIDFUNC FORT(create_procmesh_2d, CREATE_PROCMESH_2D)(
    int *ierror,
    MPI_Fint *comm, int *np2, int *np3,
    MPI_Fint *comm_cart_2d
    );

