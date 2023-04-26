#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "monolis.h"

void monolis_def_solver_c_test(){
  MONOLIS mat;
  int i;

  monolis_std_log_string("monolis_def_solver_c_test");

  monolis_prm_initialize(&mat.prm);

  monolis_test_check_eq_I1("monolis_def_solver_c_test 1", mat.prm.Iarray[MONOLIS_PRM_I_METHOD], 1);
  monolis_test_check_eq_I1("monolis_def_solver_c_test 2", mat.prm.Iarray[MONOLIS_PRM_I_PRECOND], 1);
  monolis_test_check_eq_I1("monolis_def_solver_c_test 3", mat.prm.Iarray[MONOLIS_PRM_I_MAX_ITER], 1000);
  monolis_test_check_eq_I1("monolis_def_solver_c_test 4", mat.prm.Iarray[MONOLIS_PRM_I_IS_INIT_X], 1);
  monolis_test_check_eq_I1("monolis_def_solver_c_test 5", mat.prm.Iarray[MONOLIS_PRM_I_IS_SYM_MATRIX], 0);
  monolis_test_check_eq_I1("monolis_def_solver_c_test 6", mat.prm.Iarray[MONOLIS_PRM_I_IS_DEBUG], 0);
  monolis_test_check_eq_I1("monolis_def_solver_c_test 7", mat.prm.Iarray[MONOLIS_PRM_I_IS_MEASUREMENT], 0);
  monolis_test_check_eq_I1("monolis_def_solver_c_test 8", mat.prm.Iarray[MONOLIS_PRM_I_IS_CHECK_DIAG], 0);
  monolis_test_check_eq_I1("monolis_def_solver_c_test 9", mat.prm.Iarray[MONOLIS_PRM_I_IS_PREC_STORED], 0);
  monolis_test_check_eq_I1("monolis_def_solver_c_test 10", mat.prm.Iarray[MONOLIS_PRM_I_SHOW_ITERLOG], 1);
  monolis_test_check_eq_I1("monolis_def_solver_c_test 11", mat.prm.Iarray[MONOLIS_PRM_I_SHOW_TIME], 1);
  monolis_test_check_eq_I1("monolis_def_solver_c_test 12", mat.prm.Iarray[MONOLIS_PRM_I_SHOW_SUMMARY], 1);
  monolis_test_check_eq_I1("monolis_def_solver_c_test 13", mat.prm.Iarray[MONOLIS_PRM_I_SHOW_TIME_STATISTICS], 0);

  monolis_prm_finalize(&mat.prm);

  for (i = 0; i < 100; ++i){
    monolis_test_check_eq_I1("monolis_def_solver_c_test I", mat.prm.Iarray[i], 0);
  }

  for (i = 0; i < 100; ++i){
    monolis_test_check_eq_R1("monolis_def_solver_c_test R", mat.prm.Rarray[i], 0.0);
  }
}