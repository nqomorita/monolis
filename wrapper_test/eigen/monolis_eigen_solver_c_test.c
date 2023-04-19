#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "monolis.h"
#include "monolis_utils.h"
#include "monolis_eigen_solver_c.h"

void monolis_eigen_solve_c_test(){
  MONOLIS mat;
  MONOLIS_COM com;
  int n_node;
  int n_base;
  int n_dof;
  int n_elem;
  int n_get_eigen;
  int i, j, k;
  int** elem;
  double val;
  double eig_val[6];
  double** eig_mode;
  bool is_bc[6];

  monolis_std_log_string("monolis_eigen_solve_c_test");

  n_node = 5;
  n_base = 2;
  n_dof  = 1;
  n_elem = 4;

  n_get_eigen = 5;

  elem = monolis_alloc_I_2d(elem, n_node, n_base);

  elem[0][0] = 0; elem[0][1] = 1;
  elem[1][0] = 1; elem[1][1] = 2;
  elem[2][0] = 2; elem[2][1] = 3;
  elem[3][0] = 3; elem[3][1] = 4;

  monolis_initialize(&mat);
  monolis_com_initialize_by_self(&com);

  monolis_get_nonzero_pattern_by_simple_mesh_R(&mat, n_node, n_base, n_dof, n_elem, elem);

  for(i = 0; i < n_node; ++i){
    monolis_set_scalar_to_sparse_matrix_R(&mat, i, i, 0, 0, 2.0);
  }

  for(i = 0; i < 4; ++i){
    monolis_set_scalar_to_sparse_matrix_R(&mat, i, i + 1, 0, 0, 1.0);
    monolis_set_scalar_to_sparse_matrix_R(&mat, i + 1, i, 0, 0, 1.0);
  }

  for(i = 0; i < n_node; ++i){
    is_bc[i] = false;
  }

  eig_mode = monolis_alloc_R_2d(eig_mode, n_get_eigen, n_node);

  monolis_eigen_inverted_standard_lanczos_R(&mat, &com, &n_get_eigen, 1.0e-6, 100, eig_val, eig_mode, is_bc);

  monolis_test_check_eq_R1("monolis_eigen_inverted_standard_lanczos_R 1 a", eig_val[0], 0.267949192431122);
  monolis_test_check_eq_R1("monolis_eigen_inverted_standard_lanczos_R 2 a", eig_val[1], 1.0);
  monolis_test_check_eq_R1("monolis_eigen_inverted_standard_lanczos_R 3 a", eig_val[2], 2.0);
  monolis_test_check_eq_R1("monolis_eigen_inverted_standard_lanczos_R 4 a", eig_val[3], 3.0);
  monolis_test_check_eq_R1("monolis_eigen_inverted_standard_lanczos_R 5 a", eig_val[4], 3.732050807568877);

  monolis_test_check_eq_R1("monolis_eigen_inverted_standard_lanczos_R 1 b", fabs(eig_mode[0][0]), 0.28867513459481281);
  monolis_test_check_eq_R1("monolis_eigen_inverted_standard_lanczos_R 2 b", fabs(eig_mode[0][1]),fabs(-0.5));
  monolis_test_check_eq_R1("monolis_eigen_inverted_standard_lanczos_R 3 b", fabs(eig_mode[0][2]), 0.57735026918962640);
  monolis_test_check_eq_R1("monolis_eigen_inverted_standard_lanczos_R 4 b", fabs(eig_mode[0][3]),fabs(-0.5));
  monolis_test_check_eq_R1("monolis_eigen_inverted_standard_lanczos_R 5 b", fabs(eig_mode[0][4]), 0.28867513459481270);

  monolis_test_check_eq_R1("monolis_eigen_inverted_standard_lanczos_R 1 c", fabs(eig_mode[1][0]), 0.5);
  monolis_test_check_eq_R1("monolis_eigen_inverted_standard_lanczos_R 2 c", fabs(eig_mode[1][1]),fabs(-0.5));
  monolis_test_check_eq_R1("monolis_eigen_inverted_standard_lanczos_R 3 c", 0.0, fabs(eig_mode[1][2]));
  monolis_test_check_eq_R1("monolis_eigen_inverted_standard_lanczos_R 4 c", fabs(eig_mode[1][3]), 0.5);
  monolis_test_check_eq_R1("monolis_eigen_inverted_standard_lanczos_R 5 c", fabs(eig_mode[1][4]),fabs(-0.5));

  monolis_test_check_eq_R1("monolis_eigen_inverted_standard_lanczos_R 1 d", fabs(eig_mode[2][0]), 0.5773502691896257);
  monolis_test_check_eq_R1("monolis_eigen_inverted_standard_lanczos_R 2 d", 0.0, fabs(eig_mode[2][1]));
  monolis_test_check_eq_R1("monolis_eigen_inverted_standard_lanczos_R 3 d", fabs(eig_mode[2][2]),fabs(-0.5773502691896257));
  monolis_test_check_eq_R1("monolis_eigen_inverted_standard_lanczos_R 4 d", 0.0, fabs(eig_mode[2][3]));
  monolis_test_check_eq_R1("monolis_eigen_inverted_standard_lanczos_R 5 d", fabs(eig_mode[2][4]), 0.5773502691896257);

  monolis_test_check_eq_R1("monolis_eigen_inverted_standard_lanczos_R 1 f", fabs(eig_mode[3][0]), 0.5);
  monolis_test_check_eq_R1("monolis_eigen_inverted_standard_lanczos_R 2 f", fabs(eig_mode[3][1]), 0.5);
  monolis_test_check_eq_R1("monolis_eigen_inverted_standard_lanczos_R 3 f", 0.0, fabs(eig_mode[3][2]));
  monolis_test_check_eq_R1("monolis_eigen_inverted_standard_lanczos_R 4 f", fabs(eig_mode[3][3]),fabs(-0.5));
  monolis_test_check_eq_R1("monolis_eigen_inverted_standard_lanczos_R 5 f", fabs(eig_mode[3][4]),fabs(-0.5));

  monolis_test_check_eq_R1("monolis_eigen_inverted_standard_lanczos_R 1 g", fabs(eig_mode[4][0]),fabs(-0.28867513459481281));
  monolis_test_check_eq_R1("monolis_eigen_inverted_standard_lanczos_R 2 g", fabs(eig_mode[4][1]),fabs(-0.5));
  monolis_test_check_eq_R1("monolis_eigen_inverted_standard_lanczos_R 3 g", fabs(eig_mode[4][2]),fabs(-0.57735026918962640));
  monolis_test_check_eq_R1("monolis_eigen_inverted_standard_lanczos_R 4 g", fabs(eig_mode[4][3]),fabs(-0.5));
  monolis_test_check_eq_R1("monolis_eigen_inverted_standard_lanczos_R 5 g", fabs(eig_mode[4][4]),fabs(-0.28867513459481270));

  for(i = 0; i < n_get_eigen; ++i){
    eig_val[i] = 0.0;
  }

  for(i = 0; i < n_get_eigen; ++i){
    for(j = 0; j < n_node; ++j){
      eig_mode[i][j] = 0.0;
    }
  }

  monolis_std_log_string("monolis_eigen_standard_lanczos_R");

  monolis_eigen_standard_lanczos_R(&mat, &com, &n_get_eigen, 1.0e-6, 100, eig_val, eig_mode, is_bc);

  monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 1 a", eig_val[4], 0.267949192431122);
  monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 2 a", eig_val[3], 1.0);
  monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 3 a", eig_val[2], 2.0);
  monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 4 a", eig_val[1], 3.0);
  monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 5 a", eig_val[0], 3.732050807568877);

  monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 1 b", fabs(eig_mode[4][0]), 0.28867513459481281);
  monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 2 b", fabs(eig_mode[4][1]),fabs(-0.5));
  monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 3 b", fabs(eig_mode[4][2]), 0.57735026918962640);
  monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 4 b", fabs(eig_mode[4][3]),fabs(-0.5));
  monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 5 b", fabs(eig_mode[4][4]), 0.28867513459481270);

  monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 1 c", fabs(eig_mode[3][0]), 0.5);
  monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 2 c", fabs(eig_mode[3][1]),fabs(-0.5));
  monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 3 c", 0.0, fabs(eig_mode[3][2]));
  monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 4 c", fabs(eig_mode[3][3]), 0.5);
  monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 5 c", fabs(eig_mode[3][4]),fabs(-0.5));

  monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 1 d", fabs(eig_mode[2][0]), 0.5773502691896257);
  monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 2 d", 0.0, fabs(eig_mode[2][1]));
  monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 3 d", fabs(eig_mode[2][2]),fabs(-0.5773502691896257));
  monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 4 d", 0.0, fabs(eig_mode[2][3]));
  monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 5 d", fabs(eig_mode[2][4]), 0.5773502691896257);

  monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 1 f", fabs(eig_mode[1][0]), 0.5);
  monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 2 f", fabs(eig_mode[1][1]), 0.5);
  monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 3 f", 0.0, fabs(eig_mode[1][2]));
  monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 4 f", fabs(eig_mode[1][3]),fabs(-0.5));
  monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 5 f", fabs(eig_mode[1][4]),fabs(-0.5));

  monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 1 f", fabs(eig_mode[0][0]),fabs(-0.28867513459481281));
  monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 2 f", fabs(eig_mode[0][1]),fabs(-0.5));
  monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 3 f", fabs(eig_mode[0][2]),fabs(-0.57735026918962640));
  monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 4 f", fabs(eig_mode[0][3]),fabs(-0.5));
  monolis_test_check_eq_R1("monolis_eigen_standard_lanczos_R 5 f", fabs(eig_mode[0][4]),fabs(-0.28867513459481270));

  monolis_finalize(&mat);
}
