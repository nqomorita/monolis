#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "monolis.h"
#include "monolis_utils.h"

void monolis_def_solver_util_c_test(){
  MONOLIS mat;
  int i_param;
  double r_param;
  bool b_param;

  monolis_std_log_string("monolis_def_solver_util_c_test");

  monolis_set_method(&mat, i_param);

  monolis_set_precond(&mat, i_param);

  monolis_set_maxiter(&mat, i_param);

  monolis_get_converge_iter(&mat, &i_param);

  monolis_get_error_tag(&mat, &i_param);

  monolis_set_init_x(&mat, b_param);

  monolis_set_sym_matrix(&mat, b_param);

  monolis_set_debug(&mat, b_param);

  monolis_set_performance_measurement(&mat, b_param);

  monolis_set_check_diag(&mat, b_param);

  monolis_set_prec_stored(&mat, b_param);

  monolis_show_iterlog(&mat, b_param);

  monolis_show_timelog(&mat, b_param);

  monolis_show_summary(&mat, b_param);

  monolis_show_timelog_statistics(&mat, b_param);

  monolis_set_tolerance(&mat, r_param);

  monolis_get_converge_residual(&mat, &r_param);

  monolis_get_time_solver(&mat, &r_param);

  monolis_get_time_preparing(&mat, &r_param);

  monolis_get_time_spmv(&mat, &r_param);

  monolis_get_time_inner_product(&mat, &r_param);

  monolis_get_time_precondition(&mat, &r_param);

  monolis_get_time_comm_inner_product(&mat, &r_param);

  monolis_get_time_comm_spmv(&mat, &r_param);
}