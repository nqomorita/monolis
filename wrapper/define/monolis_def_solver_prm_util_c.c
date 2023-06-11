#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "monolis_def_solver_prm_util_c.h"
#include "monolis_def_solver_prm_c.h"

void monolis_set_method(
  MONOLIS* mat,
  int      param)
{
  mat->prm.Iarray[MONOLIS_PRM_I_METHOD] = param;
}

void monolis_set_precond(
  MONOLIS* mat,
  int      param)
{
  mat->prm.Iarray[MONOLIS_PRM_I_PRECOND] = param;
}

void monolis_set_maxiter(
  MONOLIS* mat,
  int      param)
{
  mat->prm.Iarray[MONOLIS_PRM_I_MAX_ITER] = param;
}

void monolis_get_converge_iter(
  MONOLIS* mat,
  int*     param)
{
  *param = mat->prm.Iarray[MONOLIS_PRM_I_CUR_ITER];
}

void monolis_get_error_tag(
  MONOLIS* mat,
  int*     param)
{
  *param = mat->prm.Iarray[MONOLIS_PRM_I_IERR];
}

void monolis_set_init_x(
  MONOLIS* mat,
  bool     param)
{
  mat->prm.Iarray[MONOLIS_PRM_I_IS_INIT_X] = (int)param;
}

void monolis_set_sym_matrix(
  MONOLIS* mat,
  bool     param)
{
  mat->prm.Iarray[MONOLIS_PRM_I_IS_SYM_MATRIX] = param;
}

void monolis_set_debug(
  MONOLIS* mat,
  bool     param)
{
  mat->prm.Iarray[MONOLIS_PRM_I_IS_DEBUG] = param;
}

void monolis_set_performance_measurement(
  MONOLIS* mat,
  bool     param)
{
  mat->prm.Iarray[MONOLIS_PRM_I_IS_MEASUREMENT] = param;
}

void monolis_set_check_diag(
  MONOLIS* mat,
  bool     param)
{
  mat->prm.Iarray[MONOLIS_PRM_I_IS_CHECK_DIAG] = param;
}

void monolis_set_prec_stored(
  MONOLIS* mat,
  bool     param)
{
  mat->prm.Iarray[MONOLIS_PRM_I_IS_PREC_STORED] = param;
}

void monolis_show_iterlog(
  MONOLIS* mat,
  bool     param)
{
  mat->prm.Iarray[MONOLIS_PRM_I_SHOW_ITERLOG] = param;
}

void monolis_show_timelog(
  MONOLIS* mat,
  bool     param)
{
  mat->prm.Iarray[MONOLIS_PRM_I_SHOW_TIME] = param;
}

void monolis_show_summary(
  MONOLIS* mat,
  bool     param)
{
  mat->prm.Iarray[MONOLIS_PRM_I_SHOW_SUMMARY] = param;
}

void monolis_show_timelog_statistics(
  MONOLIS* mat,
  bool     param)
{
  mat->prm.Iarray[MONOLIS_PRM_I_SHOW_TIME_STATISTICS] = param;
}

void monolis_set_tolerance(
  MONOLIS* mat,
  double   val)
{
  mat->prm.Rarray[MONOLIS_PRM_R_TOL] = val;
}

void monolis_get_converge_residual(
  MONOLIS* mat,
  double*  val)
{
  *val = mat->prm.Rarray[MONOLIS_PRM_R_CUR_RESID];
}

void monolis_get_time_solver(
  MONOLIS* mat,
  double*  val)
{
  *val = mat->prm.Rarray[MONOLIS_R_TIME_SOL];
}

void monolis_get_time_preparing(
  MONOLIS* mat,
  double*  val)
{
  *val = mat->prm.Rarray[MONOLIS_R_TIME_PREP];
}

void monolis_get_time_spmv(
  MONOLIS* mat,
  double*  val)
{
  *val = mat->prm.Rarray[MONOLIS_R_TIME_SPMV];
}

void monolis_get_time_inner_product(
  MONOLIS* mat,
  double*  val)
{
  *val = mat->prm.Rarray[MONOLIS_R_TIME_DOTP];
}

void monolis_get_time_precondition(
  MONOLIS* mat,
  double*  val)
{
  *val = mat->prm.Rarray[MONOLIS_R_TIME_PREC];
}

void monolis_get_time_comm_inner_product(
  MONOLIS* mat,
  double*  val)
{
  *val = mat->prm.Rarray[MONOLIS_R_TIME_COMM_DOTP];
}

void monolis_get_time_comm_spmv(
  MONOLIS* mat,
  double*  val)
{
  *val = mat->prm.Rarray[MONOLIS_R_TIME_COMM_SPMV];
}
