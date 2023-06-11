#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "monolis.h"
#include "monolis_utils.h"

void clear_Iarray(
  MONOLIS *mat){
  int i;
  for (i = 0; i < 100; ++i){
    mat->prm.Iarray[i] = 0;
  }
}

void monolis_def_solver_util_c_test(){
  MONOLIS mat;
  int i_param;
  double r_param;
  bool b_param;

  monolis_std_log_string("monolis_def_solver_util_c_test");

  i_param = 1;
  monolis_set_method(&mat, i_param);

  i_param = 2;
  monolis_set_precond(&mat, i_param);

  i_param = 3;
  monolis_set_maxiter(&mat, i_param);

  monolis_test_check_eq_I1("monolis_def_solver_util_c_test 1", mat.prm.Iarray[MONOLIS_PRM_I_METHOD], 1);
  monolis_test_check_eq_I1("monolis_def_solver_util_c_test 2", mat.prm.Iarray[MONOLIS_PRM_I_PRECOND], 2);
  monolis_test_check_eq_I1("monolis_def_solver_util_c_test 3", mat.prm.Iarray[MONOLIS_PRM_I_MAX_ITER], 3);

  mat.prm.Iarray[MONOLIS_PRM_I_CUR_ITER] = 4;
  monolis_get_converge_iter(&mat, &i_param);
  monolis_test_check_eq_I1("monolis_def_solver_util_c_test 4", i_param, 4);

  mat.prm.Iarray[MONOLIS_PRM_I_IERR] = 5;
  monolis_get_error_tag(&mat, &i_param);
  monolis_test_check_eq_I1("monolis_def_solver_util_c_test 5", i_param, 5);

  clear_Iarray(&mat);
  b_param = true;
  monolis_set_init_x(&mat, b_param);
  monolis_test_check_eq_I1("monolis_def_solver_util_c_test 6", mat.prm.Iarray[MONOLIS_PRM_I_IS_INIT_X], 1);

  clear_Iarray(&mat);
  b_param = true;
  monolis_set_sym_matrix(&mat, b_param);
  monolis_test_check_eq_I1("monolis_def_solver_util_c_test 7", mat.prm.Iarray[MONOLIS_PRM_I_IS_SYM_MATRIX], 1);

  clear_Iarray(&mat);
  b_param = true;
  monolis_set_debug(&mat, b_param);
  monolis_test_check_eq_I1("monolis_def_solver_util_c_test 8", mat.prm.Iarray[MONOLIS_PRM_I_IS_DEBUG], 1);

  clear_Iarray(&mat);
  b_param = true;
  monolis_set_performance_measurement(&mat, b_param);
  monolis_test_check_eq_I1("monolis_def_solver_util_c_test 9", mat.prm.Iarray[MONOLIS_PRM_I_IS_MEASUREMENT], 1);

  clear_Iarray(&mat);
  b_param = true;
  monolis_set_check_diag(&mat, b_param);
  monolis_test_check_eq_I1("monolis_def_solver_util_c_test 10", mat.prm.Iarray[MONOLIS_PRM_I_IS_CHECK_DIAG], 1);

  clear_Iarray(&mat);
  b_param = true;
  monolis_set_prec_stored(&mat, b_param);
  monolis_test_check_eq_I1("monolis_def_solver_util_c_test 11", mat.prm.Iarray[MONOLIS_PRM_I_IS_PREC_STORED], 1);

  clear_Iarray(&mat);
  b_param = true;
  monolis_show_iterlog(&mat, b_param);
  monolis_test_check_eq_I1("monolis_def_solver_util_c_test 12", mat.prm.Iarray[MONOLIS_PRM_I_SHOW_ITERLOG], 1);

  clear_Iarray(&mat);
  b_param = true;
  monolis_show_timelog(&mat, b_param);
  monolis_test_check_eq_I1("monolis_def_solver_util_c_test 13", mat.prm.Iarray[MONOLIS_PRM_I_SHOW_TIME], 1);

  clear_Iarray(&mat);
  b_param = true;
  monolis_show_summary(&mat, b_param);
  monolis_test_check_eq_I1("monolis_def_solver_util_c_test 14", mat.prm.Iarray[MONOLIS_PRM_I_SHOW_SUMMARY], 1);

  clear_Iarray(&mat);
  b_param = true;
  monolis_show_timelog_statistics(&mat, b_param);
  monolis_test_check_eq_I1("monolis_def_solver_util_c_test 15", mat.prm.Iarray[MONOLIS_PRM_I_SHOW_TIME_STATISTICS], 1);

  r_param = 1.0;
  monolis_set_tolerance(&mat, r_param);
  monolis_test_check_eq_R1("monolis_def_solver_util_c_test 16", mat.prm.Rarray[MONOLIS_PRM_R_TOL], 1.0);

  r_param = 0.0;
  mat.prm.Rarray[MONOLIS_PRM_R_CUR_RESID] = 2.0;
  monolis_get_converge_residual(&mat, &r_param);
  monolis_test_check_eq_R1("monolis_def_solver_util_c_test 17", r_param, 2.0);

  r_param = 0.0;
  mat.prm.Rarray[MONOLIS_R_TIME_SOL] = 3.0;
  monolis_get_time_solver(&mat, &r_param);
  monolis_test_check_eq_R1("monolis_def_solver_util_c_test 18", r_param, 3.0);

  r_param = 0.0;
  mat.prm.Rarray[MONOLIS_R_TIME_PREP] = 4.0;
  monolis_get_time_preparing(&mat, &r_param);
  monolis_test_check_eq_R1("monolis_def_solver_util_c_test 19", r_param, 4.0);

  r_param = 0.0;
  mat.prm.Rarray[MONOLIS_R_TIME_SPMV] = 5.0;
  monolis_get_time_spmv(&mat, &r_param);
  monolis_test_check_eq_R1("monolis_def_solver_util_c_test 20", r_param, 5.0);

  r_param = 0.0;
  mat.prm.Rarray[MONOLIS_R_TIME_DOTP] = 6.0;
  monolis_get_time_inner_product(&mat, &r_param);
  monolis_test_check_eq_R1("monolis_def_solver_util_c_test 21", r_param, 6.0);

  r_param = 0.0;
  mat.prm.Rarray[MONOLIS_R_TIME_PREC] = 7.0;
  monolis_get_time_precondition(&mat, &r_param);
  monolis_test_check_eq_R1("monolis_def_solver_util_c_test 22", r_param, 7.0);

  r_param = 0.0;
  mat.prm.Rarray[MONOLIS_R_TIME_COMM_DOTP] = 8.0;
  monolis_get_time_comm_inner_product(&mat, &r_param);
  monolis_test_check_eq_R1("monolis_def_solver_util_c_test 23", r_param, 8.0);

  r_param = 0.0;
  mat.prm.Rarray[MONOLIS_R_TIME_COMM_SPMV] = 9.0;
  monolis_get_time_comm_spmv(&mat, &r_param);
  monolis_test_check_eq_R1("monolis_def_solver_util_c_test 24", r_param, 9.0);
}
