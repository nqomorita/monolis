#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "monolis_def_solver_util_c.h"
#include "monolis_def_solver_c.h"

void monolis_set_method(
  MONOLIS* mat,
  int      param)
{
  mat->prm.Iarray[monolis_prm_I_method] = param;
}

void monolis_set_precond(
  MONOLIS* mat,
  int      param)
{
  mat->prm.Iarray[monolis_prm_I_precond] = param;
}

void monolis_set_maxiter(
  MONOLIS* mat,
  int      param)
{
  mat->prm.Iarray[monolis_prm_I_max_iter] = param;
}

void monolis_get_converge_iter(
  MONOLIS* mat,
  int*     param)
{
  param = &mat->prm.Iarray[monolis_prm_I_cur_iter];
}

void monolis_get_error_tag(
  MONOLIS* mat,
  int*     param)
{
  param = &mat->prm.Iarray[monolis_prm_I_ierr];
}

void monolis_set_init_x(
  MONOLIS* mat,
  bool     param)
{
  mat->prm.Iarray[monolis_prm_I_is_init_x] = (int)param;
}

void monolis_set_sym_matrix(
  MONOLIS* mat,
  bool     param)
{
  mat->prm.Iarray[monolis_prm_I_is_sym_matrix] = param;
}

void monolis_set_debug(
  MONOLIS* mat,
  bool     param)
{
  mat->prm.Iarray[monolis_prm_I_is_debug] = param;
}

void monolis_set_performance_measurement(
  MONOLIS* mat,
  bool     param)
{
  mat->prm.Iarray[monolis_prm_I_is_measurement] = param;
}

void monolis_set_check_diag(
  MONOLIS* mat,
  bool     param)
{
  mat->prm.Iarray[monolis_prm_I_is_check_diag] = param;
}

void monolis_set_prec_stored(
  MONOLIS* mat,
  bool     param)
{
  mat->prm.Iarray[monolis_prm_I_is_prec_stored] = param;
}

void monolis_show_iterlog(
  MONOLIS* mat,
  bool     param)
{
  mat->prm.Iarray[monolis_prm_I_show_iterlog] = param;
}

void monolis_show_timelog(
  MONOLIS* mat,
  bool     param)
{
  mat->prm.Iarray[monolis_prm_I_show_time] = param;
}

void monolis_show_summary(
  MONOLIS* mat,
  bool     param)
{
  mat->prm.Iarray[monolis_prm_I_show_summary] = param;
}

void monolis_show_timelog_statistics(
  MONOLIS* mat,
  bool     param)
{
  mat->prm.Iarray[monolis_prm_I_show_time_statistics] = param;
}

void monolis_set_tolerance(
  MONOLIS* mat,
  double   val)
{
  mat->prm.Rarray[monolis_prm_R_tol] = val;
}

void monolis_get_converge_residual(
  MONOLIS* mat,
  double*  val)
{
  val = &mat->prm.Rarray[monolis_prm_I_method];
}

void monolis_get_time_solver(
  MONOLIS* mat,
  double*  val)
{
  val = &mat->prm.Rarray[monolis_R_time_sol];
}

void monolis_get_time_preparing(
  MONOLIS* mat,
  double*  val)
{
  val = &mat->prm.Rarray[monolis_R_time_prep];
}

void monolis_get_time_spmv(
  MONOLIS* mat,
  double*  val)
{
  val = &mat->prm.Rarray[monolis_R_time_spmv];
}

void monolis_get_time_inner_product(
  MONOLIS* mat,
  double*  val)
{
  val = &mat->prm.Rarray[monolis_R_time_dotp];
}

void monolis_get_time_precondition(
  MONOLIS* mat,
  double*  val)
{
  val = &mat->prm.Rarray[monolis_R_time_prec];
}

void monolis_get_time_comm_inner_product(
  MONOLIS* mat,
  double*  val)
{
  val = &mat->prm.Rarray[monolis_R_time_comm_dotp];
}

void monolis_get_time_comm_spmv(
  MONOLIS* mat,
  double*  val)
{
  val = &mat->prm.Rarray[monolis_R_time_comm_spmv];
}
