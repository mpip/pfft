/* Generated automatically.  DO NOT EDIT! */

#include "pfft.h"
#include "ipfft.h"

PFFT_EXTERN R PX(check_output_c2c_3d_f03)(const INT * Nos, const INT * local_n, const INT * local_n_start, const C * data, MPI_Fint f_comm);
PFFT_EXTERN R PX(check_output_c2c_f03)(int rnk_n, const INT * Nos, const INT * local_n, const INT * local_start, const C * data, MPI_Fint f_comm);
PFFT_EXTERN R PX(check_output_c2r_3d_f03)(const INT * Nos, const INT * local_n, const INT * local_n_start, const R * data, MPI_Fint f_comm);
PFFT_EXTERN R PX(check_output_c2r_f03)(int rnk_n, const INT * Nos, const INT * local_n, const INT * local_start, const R * data, MPI_Fint f_comm);
PFFT_EXTERN R PX(check_output_r2r_3d_f03)(const INT * Nos, const INT * local_n, const INT * local_n_start, const R * data, MPI_Fint f_comm);
PFFT_EXTERN R PX(check_output_r2r_f03)(int rnk_n, const INT * Nos, const INT * local_n, const INT * local_start, const R * data, MPI_Fint f_comm);
PFFT_EXTERN INT PX(local_size_dft_3d_f03)(const INT * Nos, MPI_Fint f_comm_cart, unsigned pfft_flags, INT * local_ni, INT * local_i_start, INT * local_no, INT * local_o_start);
PFFT_EXTERN INT PX(local_size_dft_r2c_3d_f03)(const INT * Nos, MPI_Fint f_comm_cart, unsigned pfft_flags, INT * local_ni, INT * local_i_start, INT * local_no, INT * local_o_start);
PFFT_EXTERN INT PX(local_size_dft_c2r_3d_f03)(const INT * Nos, MPI_Fint f_comm_cart, unsigned pfft_flags, INT * local_ni, INT * local_i_start, INT * local_no, INT * local_o_start);
PFFT_EXTERN INT PX(local_size_r2r_3d_f03)(const INT * Nos, MPI_Fint f_comm_cart, unsigned pfft_flags, INT * local_ni, INT * local_i_start, INT * local_no, INT * local_o_start);
PFFT_EXTERN INT PX(local_size_dft_f03)(int rnk, const INT * Nos, MPI_Fint f_comm_cart, unsigned pfft_flags, INT * local_ni, INT * local_i_start, INT * local_no, INT * local_o_start);
PFFT_EXTERN INT PX(local_size_dft_r2c_f03)(int rnk, const INT * Nos, MPI_Fint f_comm_cart, unsigned pfft_flags, INT * local_ni, INT * local_i_start, INT * local_no, INT * local_o_start);
PFFT_EXTERN INT PX(local_size_dft_c2r_f03)(int rnk, const INT * Nos, MPI_Fint f_comm_cart, unsigned pfft_flags, INT * local_ni, INT * local_i_start, INT * local_no, INT * local_o_start);
PFFT_EXTERN INT PX(local_size_r2r_f03)(int rnk_n, const INT * Nos, MPI_Fint f_comm_cart, unsigned pfft_flags, INT * local_ni, INT * local_i_start, INT * local_no, INT * local_o_start);
PFFT_EXTERN INT PX(local_size_many_dft_f03)(int rnk, const INT * Nos, const INT * ni, const INT * no, INT howmany, const INT * iblock, const INT * oblock, MPI_Fint f_comm_cart, unsigned pfft_flags, INT * local_ni, INT * local_i_start, INT * local_no, INT * local_o_start);
PFFT_EXTERN INT PX(local_size_many_dft_r2c_f03)(int rnk_n, const INT * Nos, const INT * ni, const INT * no, INT howmany, const INT * iblock, const INT * oblock, MPI_Fint f_comm_cart, unsigned pfft_flags, INT * local_ni, INT * local_i_start, INT * local_no, INT * local_o_start);
PFFT_EXTERN INT PX(local_size_many_dft_c2r_f03)(int rnk_n, const INT * Nos, const INT * ni, const INT * no, INT howmany, const INT * iblock, const INT * oblock, MPI_Fint f_comm_cart, unsigned pfft_flags, INT * local_ni, INT * local_i_start, INT * local_no, INT * local_o_start);
PFFT_EXTERN INT PX(local_size_many_r2r_f03)(int rnk_n, const INT * Nos, const INT * ni, const INT * no, INT howmany, const INT * iblock, const INT * oblock, MPI_Fint f_comm_cart, unsigned pfft_flags, INT * local_ni, INT * local_i_start, INT * local_no, INT * local_o_start);
PFFT_EXTERN PX(plan) PX(plan_dft_3d_f03)(const INT * Nos, C * in, C * out, MPI_Fint f_comm_cart, int sign, unsigned pfft_flags);
PFFT_EXTERN PX(plan) PX(plan_dft_r2c_3d_f03)(const INT * Nos, R * in, C * out, MPI_Fint f_comm_cart, int sign, unsigned pfft_flags);
PFFT_EXTERN PX(plan) PX(plan_dft_c2r_3d_f03)(const INT * Nos, C * in, R * out, MPI_Fint f_comm_cart, int sign, unsigned pfft_flags);
PFFT_EXTERN PX(plan) PX(plan_r2r_3d_f03)(const INT * Nos, R * in, R * out, MPI_Fint f_comm_cart, const PX(r2r_kind) * kinds, unsigned pfft_flags);
PFFT_EXTERN PX(plan) PX(plan_dft_f03)(int rnk, const INT * Nos, C * in, C * out, MPI_Fint f_comm_cart, int sign, unsigned pfft_flags);
PFFT_EXTERN PX(plan) PX(plan_dft_r2c_f03)(int rnk_n, const INT * Nos, R * in, C * out, MPI_Fint f_comm_cart, int sign, unsigned pfft_flags);
PFFT_EXTERN PX(plan) PX(plan_dft_c2r_f03)(int rnk_n, const INT * Nos, C * in, R * out, MPI_Fint f_comm_cart, int sign, unsigned pfft_flags);
PFFT_EXTERN PX(plan) PX(plan_many_dft_f03)(int rnk, const INT * Nos, const INT * ni, const INT * no, INT howmany, const INT * iblock, const INT * oblock, C * in, C * out, MPI_Fint f_comm_cart, int sign, unsigned pfft_flags);
PFFT_EXTERN PX(plan) PX(plan_many_dft_r2c_f03)(int rnk_n, const INT * Nos, const INT * ni, const INT * no, INT howmany, const INT * iblock, const INT * oblock, R * in, C * out, MPI_Fint f_comm_cart, int sign, unsigned pfft_flags);
PFFT_EXTERN PX(plan) PX(plan_many_dft_c2r_f03)(int rnk_n, const INT * Nos, const INT * ni, const INT * no, INT howmany, const INT * iblock, const INT * oblock, C * in, R * out, MPI_Fint f_comm_cart, int sign, unsigned pfft_flags);
PFFT_EXTERN PX(plan) PX(plan_many_r2r_f03)(int rnk_n, const INT * Nos, const INT * ni, const INT * no, INT howmany, const INT * iblock, const INT * oblock, R * in, R * out, MPI_Fint f_comm_cart, const PX(r2r_kind) * kinds, unsigned pfft_flags);
PFFT_EXTERN PX(plan) PX(plan_many_dft_skipped_f03)(int rnk, const INT * Nos, const INT * ni, const INT * no, INT howmany, const INT * iblock, const INT * oblock, const int * skip_trafos, C * in, C * out, MPI_Fint f_comm_cart, int sign, unsigned pfft_flags);
PFFT_EXTERN void PX(print_average_timer_f03)(const PX(plan) ths, MPI_Fint f_comm);
PFFT_EXTERN void PX(print_average_timer_adv_f03)(const PX(plan) ths, MPI_Fint f_comm);
PFFT_EXTERN void PX(write_average_timer_f03)(const PX(plan) ths, const char * name, MPI_Fint f_comm);
PFFT_EXTERN void PX(write_average_timer_adv_f03)(const PX(plan) ths, const char * name, MPI_Fint f_comm);
PFFT_EXTERN PX(timer) PX(reduce_max_timer_f03)(const PX(timer) ths, MPI_Fint f_comm);
PFFT_EXTERN int PX(create_procmesh_f03)(int rnk, MPI_Fint f_comm, const int * np, MPI_Fint * f_comm_cart);
PFFT_EXTERN int PX(create_procmesh_2d_f03)(MPI_Fint f_comm, int np0, int np1, MPI_Fint * f_comm_cart_2d);
PFFT_EXTERN PX(gcplan) PX(plan_rgc_3d_f03)(const INT * Nos, const INT * gc_below, const INT * gc_above, R * data, MPI_Fint f_comm_cart, unsigned gc_flags);
PFFT_EXTERN PX(gcplan) PX(plan_cgc_3d_f03)(const INT * Nos, const INT * gc_below, const INT * gc_above, C * data, MPI_Fint f_comm_cart, unsigned gc_flags);
PFFT_EXTERN PX(gcplan) PX(plan_rgc_f03)(int rnk_n, const INT * Nos, const INT * gc_below, const INT * gc_above, R * data, MPI_Fint f_comm_cart, unsigned gc_flags);
PFFT_EXTERN PX(gcplan) PX(plan_cgc_f03)(int rnk_n, const INT * Nos, const INT * gc_below, const INT * gc_above, C * data, MPI_Fint f_comm_cart, unsigned gc_flags);
PFFT_EXTERN PX(gcplan) PX(plan_many_rgc_f03)(int rnk_n, const INT * Nos, INT howmany, const INT * block, const INT * gc_below, const INT * gc_above, R * data, MPI_Fint f_comm_cart, unsigned gc_flags);
PFFT_EXTERN PX(gcplan) PX(plan_many_cgc_f03)(int rnk_n, const INT * Nos, INT howmany, const INT * block, const INT * gc_below, const INT * gc_above, C * data, MPI_Fint f_comm_cart, unsigned gc_flags);
PFFT_EXTERN void PX(print_average_gctimer_f03)(const PX(gcplan) ths, MPI_Fint f_comm);
PFFT_EXTERN void PX(print_average_gctimer_adv_f03)(const PX(gcplan) ths, MPI_Fint f_comm);
PFFT_EXTERN void PX(write_average_gctimer_f03)(const PX(gcplan) ths, const char * name, MPI_Fint f_comm);
PFFT_EXTERN void PX(write_average_gctimer_adv_f03)(const PX(gcplan) ths, const char * name, MPI_Fint f_comm);
PFFT_EXTERN PX(gctimer) PX(reduce_max_gctimer_f03)(const PX(gctimer) ths, MPI_Fint f_comm);

R PX(check_output_c2c_3d_f03)(const INT * Nos, const INT * local_n, const INT * local_n_start, const C * data, MPI_Fint f_comm)
{
  MPI_Comm comm;

  comm = MPI_Comm_f2c(f_comm);
  R ret = PX(check_output_c2c_3d)(Nos, local_n, local_n_start, data, comm);
  return ret;
}

R PX(check_output_c2c_f03)(int rnk_n, const INT * Nos, const INT * local_n, const INT * local_start, const C * data, MPI_Fint f_comm)
{
  MPI_Comm comm;

  comm = MPI_Comm_f2c(f_comm);
  R ret = PX(check_output_c2c)(rnk_n, Nos, local_n, local_start, data, comm);
  return ret;
}

R PX(check_output_c2r_3d_f03)(const INT * Nos, const INT * local_n, const INT * local_n_start, const R * data, MPI_Fint f_comm)
{
  MPI_Comm comm;

  comm = MPI_Comm_f2c(f_comm);
  R ret = PX(check_output_c2r_3d)(Nos, local_n, local_n_start, data, comm);
  return ret;
}

R PX(check_output_c2r_f03)(int rnk_n, const INT * Nos, const INT * local_n, const INT * local_start, const R * data, MPI_Fint f_comm)
{
  MPI_Comm comm;

  comm = MPI_Comm_f2c(f_comm);
  R ret = PX(check_output_c2r)(rnk_n, Nos, local_n, local_start, data, comm);
  return ret;
}

R PX(check_output_r2r_3d_f03)(const INT * Nos, const INT * local_n, const INT * local_n_start, const R * data, MPI_Fint f_comm)
{
  MPI_Comm comm;

  comm = MPI_Comm_f2c(f_comm);
  R ret = PX(check_output_r2r_3d)(Nos, local_n, local_n_start, data, comm);
  return ret;
}

R PX(check_output_r2r_f03)(int rnk_n, const INT * Nos, const INT * local_n, const INT * local_start, const R * data, MPI_Fint f_comm)
{
  MPI_Comm comm;

  comm = MPI_Comm_f2c(f_comm);
  R ret = PX(check_output_r2r)(rnk_n, Nos, local_n, local_start, data, comm);
  return ret;
}

INT PX(local_size_dft_3d_f03)(const INT * Nos, MPI_Fint f_comm_cart, unsigned pfft_flags, INT * local_ni, INT * local_i_start, INT * local_no, INT * local_o_start)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  INT ret = PX(local_size_dft_3d)(Nos, comm_cart, pfft_flags, local_ni, local_i_start, local_no, local_o_start);
  return ret;
}

INT PX(local_size_dft_r2c_3d_f03)(const INT * Nos, MPI_Fint f_comm_cart, unsigned pfft_flags, INT * local_ni, INT * local_i_start, INT * local_no, INT * local_o_start)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  INT ret = PX(local_size_dft_r2c_3d)(Nos, comm_cart, pfft_flags, local_ni, local_i_start, local_no, local_o_start);
  return ret;
}

INT PX(local_size_dft_c2r_3d_f03)(const INT * Nos, MPI_Fint f_comm_cart, unsigned pfft_flags, INT * local_ni, INT * local_i_start, INT * local_no, INT * local_o_start)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  INT ret = PX(local_size_dft_c2r_3d)(Nos, comm_cart, pfft_flags, local_ni, local_i_start, local_no, local_o_start);
  return ret;
}

INT PX(local_size_r2r_3d_f03)(const INT * Nos, MPI_Fint f_comm_cart, unsigned pfft_flags, INT * local_ni, INT * local_i_start, INT * local_no, INT * local_o_start)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  INT ret = PX(local_size_r2r_3d)(Nos, comm_cart, pfft_flags, local_ni, local_i_start, local_no, local_o_start);
  return ret;
}

INT PX(local_size_dft_f03)(int rnk, const INT * Nos, MPI_Fint f_comm_cart, unsigned pfft_flags, INT * local_ni, INT * local_i_start, INT * local_no, INT * local_o_start)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  INT ret = PX(local_size_dft)(rnk, Nos, comm_cart, pfft_flags, local_ni, local_i_start, local_no, local_o_start);
  return ret;
}

INT PX(local_size_dft_r2c_f03)(int rnk, const INT * Nos, MPI_Fint f_comm_cart, unsigned pfft_flags, INT * local_ni, INT * local_i_start, INT * local_no, INT * local_o_start)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  INT ret = PX(local_size_dft_r2c)(rnk, Nos, comm_cart, pfft_flags, local_ni, local_i_start, local_no, local_o_start);
  return ret;
}

INT PX(local_size_dft_c2r_f03)(int rnk, const INT * Nos, MPI_Fint f_comm_cart, unsigned pfft_flags, INT * local_ni, INT * local_i_start, INT * local_no, INT * local_o_start)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  INT ret = PX(local_size_dft_c2r)(rnk, Nos, comm_cart, pfft_flags, local_ni, local_i_start, local_no, local_o_start);
  return ret;
}

INT PX(local_size_r2r_f03)(int rnk_n, const INT * Nos, MPI_Fint f_comm_cart, unsigned pfft_flags, INT * local_ni, INT * local_i_start, INT * local_no, INT * local_o_start)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  INT ret = PX(local_size_r2r)(rnk_n, Nos, comm_cart, pfft_flags, local_ni, local_i_start, local_no, local_o_start);
  return ret;
}

INT PX(local_size_many_dft_f03)(int rnk, const INT * Nos, const INT * ni, const INT * no, INT howmany, const INT * iblock, const INT * oblock, MPI_Fint f_comm_cart, unsigned pfft_flags, INT * local_ni, INT * local_i_start, INT * local_no, INT * local_o_start)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  INT ret = PX(local_size_many_dft)(rnk, Nos, ni, no, howmany, iblock, oblock, comm_cart, pfft_flags, local_ni, local_i_start, local_no, local_o_start);
  return ret;
}

INT PX(local_size_many_dft_r2c_f03)(int rnk_n, const INT * Nos, const INT * ni, const INT * no, INT howmany, const INT * iblock, const INT * oblock, MPI_Fint f_comm_cart, unsigned pfft_flags, INT * local_ni, INT * local_i_start, INT * local_no, INT * local_o_start)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  INT ret = PX(local_size_many_dft_r2c)(rnk_n, Nos, ni, no, howmany, iblock, oblock, comm_cart, pfft_flags, local_ni, local_i_start, local_no, local_o_start);
  return ret;
}

INT PX(local_size_many_dft_c2r_f03)(int rnk_n, const INT * Nos, const INT * ni, const INT * no, INT howmany, const INT * iblock, const INT * oblock, MPI_Fint f_comm_cart, unsigned pfft_flags, INT * local_ni, INT * local_i_start, INT * local_no, INT * local_o_start)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  INT ret = PX(local_size_many_dft_c2r)(rnk_n, Nos, ni, no, howmany, iblock, oblock, comm_cart, pfft_flags, local_ni, local_i_start, local_no, local_o_start);
  return ret;
}

INT PX(local_size_many_r2r_f03)(int rnk_n, const INT * Nos, const INT * ni, const INT * no, INT howmany, const INT * iblock, const INT * oblock, MPI_Fint f_comm_cart, unsigned pfft_flags, INT * local_ni, INT * local_i_start, INT * local_no, INT * local_o_start)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  INT ret = PX(local_size_many_r2r)(rnk_n, Nos, ni, no, howmany, iblock, oblock, comm_cart, pfft_flags, local_ni, local_i_start, local_no, local_o_start);
  return ret;
}

PX(plan) PX(plan_dft_3d_f03)(const INT * Nos, C * in, C * out, MPI_Fint f_comm_cart, int sign, unsigned pfft_flags)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  PX(plan) ret = PX(plan_dft_3d)(Nos, in, out, comm_cart, sign, pfft_flags);
  return ret;
}

PX(plan) PX(plan_dft_r2c_3d_f03)(const INT * Nos, R * in, C * out, MPI_Fint f_comm_cart, int sign, unsigned pfft_flags)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  PX(plan) ret = PX(plan_dft_r2c_3d)(Nos, in, out, comm_cart, sign, pfft_flags);
  return ret;
}

PX(plan) PX(plan_dft_c2r_3d_f03)(const INT * Nos, C * in, R * out, MPI_Fint f_comm_cart, int sign, unsigned pfft_flags)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  PX(plan) ret = PX(plan_dft_c2r_3d)(Nos, in, out, comm_cart, sign, pfft_flags);
  return ret;
}

PX(plan) PX(plan_r2r_3d_f03)(const INT * Nos, R * in, R * out, MPI_Fint f_comm_cart, const PX(r2r_kind) * kinds, unsigned pfft_flags)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  PX(plan) ret = PX(plan_r2r_3d)(Nos, in, out, comm_cart, kinds, pfft_flags);
  return ret;
}

PX(plan) PX(plan_dft_f03)(int rnk, const INT * Nos, C * in, C * out, MPI_Fint f_comm_cart, int sign, unsigned pfft_flags)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  PX(plan) ret = PX(plan_dft)(rnk, Nos, in, out, comm_cart, sign, pfft_flags);
  return ret;
}

PX(plan) PX(plan_dft_r2c_f03)(int rnk_n, const INT * Nos, R * in, C * out, MPI_Fint f_comm_cart, int sign, unsigned pfft_flags)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  PX(plan) ret = PX(plan_dft_r2c)(rnk_n, Nos, in, out, comm_cart, sign, pfft_flags);
  return ret;
}

PX(plan) PX(plan_dft_c2r_f03)(int rnk_n, const INT * Nos, C * in, R * out, MPI_Fint f_comm_cart, int sign, unsigned pfft_flags)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  PX(plan) ret = PX(plan_dft_c2r)(rnk_n, Nos, in, out, comm_cart, sign, pfft_flags);
  return ret;
}

PX(plan) PX(plan_many_dft_f03)(int rnk, const INT * Nos, const INT * ni, const INT * no, INT howmany, const INT * iblock, const INT * oblock, C * in, C * out, MPI_Fint f_comm_cart, int sign, unsigned pfft_flags)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  PX(plan) ret = PX(plan_many_dft)(rnk, Nos, ni, no, howmany, iblock, oblock, in, out, comm_cart, sign, pfft_flags);
  return ret;
}

PX(plan) PX(plan_many_dft_r2c_f03)(int rnk_n, const INT * Nos, const INT * ni, const INT * no, INT howmany, const INT * iblock, const INT * oblock, R * in, C * out, MPI_Fint f_comm_cart, int sign, unsigned pfft_flags)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  PX(plan) ret = PX(plan_many_dft_r2c)(rnk_n, Nos, ni, no, howmany, iblock, oblock, in, out, comm_cart, sign, pfft_flags);
  return ret;
}

PX(plan) PX(plan_many_dft_c2r_f03)(int rnk_n, const INT * Nos, const INT * ni, const INT * no, INT howmany, const INT * iblock, const INT * oblock, C * in, R * out, MPI_Fint f_comm_cart, int sign, unsigned pfft_flags)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  PX(plan) ret = PX(plan_many_dft_c2r)(rnk_n, Nos, ni, no, howmany, iblock, oblock, in, out, comm_cart, sign, pfft_flags);
  return ret;
}

PX(plan) PX(plan_many_r2r_f03)(int rnk_n, const INT * Nos, const INT * ni, const INT * no, INT howmany, const INT * iblock, const INT * oblock, R * in, R * out, MPI_Fint f_comm_cart, const PX(r2r_kind) * kinds, unsigned pfft_flags)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  PX(plan) ret = PX(plan_many_r2r)(rnk_n, Nos, ni, no, howmany, iblock, oblock, in, out, comm_cart, kinds, pfft_flags);
  return ret;
}

PX(plan) PX(plan_many_dft_skipped_f03)(int rnk, const INT * Nos, const INT * ni, const INT * no, INT howmany, const INT * iblock, const INT * oblock, const int * skip_trafos, C * in, C * out, MPI_Fint f_comm_cart, int sign, unsigned pfft_flags)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  PX(plan) ret = PX(plan_many_dft_skipped)(rnk, Nos, ni, no, howmany, iblock, oblock, skip_trafos, in, out, comm_cart, sign, pfft_flags);
  return ret;
}

void PX(print_average_timer_f03)(const PX(plan) ths, MPI_Fint f_comm)
{
  MPI_Comm comm;

  comm = MPI_Comm_f2c(f_comm);
  PX(print_average_timer)(ths, comm);
}

void PX(print_average_timer_adv_f03)(const PX(plan) ths, MPI_Fint f_comm)
{
  MPI_Comm comm;

  comm = MPI_Comm_f2c(f_comm);
  PX(print_average_timer_adv)(ths, comm);
}

void PX(write_average_timer_f03)(const PX(plan) ths, const char * name, MPI_Fint f_comm)
{
  MPI_Comm comm;

  comm = MPI_Comm_f2c(f_comm);
  PX(write_average_timer)(ths, name, comm);
}

void PX(write_average_timer_adv_f03)(const PX(plan) ths, const char * name, MPI_Fint f_comm)
{
  MPI_Comm comm;

  comm = MPI_Comm_f2c(f_comm);
  PX(write_average_timer_adv)(ths, name, comm);
}

PX(timer) PX(reduce_max_timer_f03)(const PX(timer) ths, MPI_Fint f_comm)
{
  MPI_Comm comm;

  comm = MPI_Comm_f2c(f_comm);
  PX(timer) ret = PX(reduce_max_timer)(ths, comm);
  return ret;
}

int PX(create_procmesh_f03)(int rnk, MPI_Fint f_comm, const int * np, MPI_Fint * f_comm_cart)
{
  MPI_Comm comm, comm_cart;

  comm = MPI_Comm_f2c(f_comm);
  int ret = PX(create_procmesh)(rnk, comm, np, &comm_cart);
  *f_comm_cart = MPI_Comm_c2f(comm_cart);
  return ret;
}

int PX(create_procmesh_2d_f03)(MPI_Fint f_comm, int np0, int np1, MPI_Fint * f_comm_cart_2d)
{
  MPI_Comm comm, comm_cart_2d;

  comm = MPI_Comm_f2c(f_comm);
  int ret = PX(create_procmesh_2d)(comm, np0, np1, &comm_cart_2d);
  *f_comm_cart_2d = MPI_Comm_c2f(comm_cart_2d);
  return ret;
}

PX(gcplan) PX(plan_rgc_3d_f03)(const INT * Nos, const INT * gc_below, const INT * gc_above, R * data, MPI_Fint f_comm_cart, unsigned gc_flags)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  PX(gcplan) ret = PX(plan_rgc_3d)(Nos, gc_below, gc_above, data, comm_cart, gc_flags);
  return ret;
}

PX(gcplan) PX(plan_cgc_3d_f03)(const INT * Nos, const INT * gc_below, const INT * gc_above, C * data, MPI_Fint f_comm_cart, unsigned gc_flags)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  PX(gcplan) ret = PX(plan_cgc_3d)(Nos, gc_below, gc_above, data, comm_cart, gc_flags);
  return ret;
}

PX(gcplan) PX(plan_rgc_f03)(int rnk_n, const INT * Nos, const INT * gc_below, const INT * gc_above, R * data, MPI_Fint f_comm_cart, unsigned gc_flags)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  PX(gcplan) ret = PX(plan_rgc)(rnk_n, Nos, gc_below, gc_above, data, comm_cart, gc_flags);
  return ret;
}

PX(gcplan) PX(plan_cgc_f03)(int rnk_n, const INT * Nos, const INT * gc_below, const INT * gc_above, C * data, MPI_Fint f_comm_cart, unsigned gc_flags)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  PX(gcplan) ret = PX(plan_cgc)(rnk_n, Nos, gc_below, gc_above, data, comm_cart, gc_flags);
  return ret;
}

PX(gcplan) PX(plan_many_rgc_f03)(int rnk_n, const INT * Nos, INT howmany, const INT * block, const INT * gc_below, const INT * gc_above, R * data, MPI_Fint f_comm_cart, unsigned gc_flags)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  PX(gcplan) ret = PX(plan_many_rgc)(rnk_n, Nos, howmany, block, gc_below, gc_above, data, comm_cart, gc_flags);
  return ret;
}

PX(gcplan) PX(plan_many_cgc_f03)(int rnk_n, const INT * Nos, INT howmany, const INT * block, const INT * gc_below, const INT * gc_above, C * data, MPI_Fint f_comm_cart, unsigned gc_flags)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  PX(gcplan) ret = PX(plan_many_cgc)(rnk_n, Nos, howmany, block, gc_below, gc_above, data, comm_cart, gc_flags);
  return ret;
}

void PX(print_average_gctimer_f03)(const PX(gcplan) ths, MPI_Fint f_comm)
{
  MPI_Comm comm;

  comm = MPI_Comm_f2c(f_comm);
  PX(print_average_gctimer)(ths, comm);
}

void PX(print_average_gctimer_adv_f03)(const PX(gcplan) ths, MPI_Fint f_comm)
{
  MPI_Comm comm;

  comm = MPI_Comm_f2c(f_comm);
  PX(print_average_gctimer_adv)(ths, comm);
}

void PX(write_average_gctimer_f03)(const PX(gcplan) ths, const char * name, MPI_Fint f_comm)
{
  MPI_Comm comm;

  comm = MPI_Comm_f2c(f_comm);
  PX(write_average_gctimer)(ths, name, comm);
}

void PX(write_average_gctimer_adv_f03)(const PX(gcplan) ths, const char * name, MPI_Fint f_comm)
{
  MPI_Comm comm;

  comm = MPI_Comm_f2c(f_comm);
  PX(write_average_gctimer_adv)(ths, name, comm);
}

PX(gctimer) PX(reduce_max_gctimer_f03)(const PX(gctimer) ths, MPI_Fint f_comm)
{
  MPI_Comm comm;

  comm = MPI_Comm_f2c(f_comm);
  PX(gctimer) ret = PX(reduce_max_gctimer)(ths, comm);
  return ret;
}
