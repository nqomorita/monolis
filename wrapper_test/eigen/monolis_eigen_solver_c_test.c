#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "monolis.h"
#include "monolis_utils.h"
#include "monolis_eigen_solver_c.h"

void monolis_eigen_solve_c_test(){
  MONOLIS mat;
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

  n_node = 6;
  n_base = 2;
  n_dof  = 1;
  n_elem = 5;

  n_get_eigen = 5;

  elem = monolis_alloc_I_2d(elem, n_node, n_base);

  elem[0][0] = 0; elem[0][1] = 1;
  elem[1][0] = 1; elem[1][1] = 2;
  elem[2][0] = 2; elem[2][1] = 3;
  elem[3][0] = 3; elem[3][1] = 4;
  elem[4][0] = 4; elem[4][1] = 5;

  monolis_initialize(&mat, "./");

  monolis_get_nonzero_pattern_by_simple_mesh_R(&mat, n_node, n_base, n_dof, n_elem, elem);

  for(int i = 0; i < 6; ++i){
    monolis_set_scalar_to_sparse_matrix_R(&mat, i, i, 0, 0, 2.0);
  }

  for(int i = 0; i < 4; ++i){
    monolis_set_scalar_to_sparse_matrix_R(&mat, i, i + 1, 0, 0, 1.0);
    monolis_set_scalar_to_sparse_matrix_R(&mat, i + 1, i, 0, 0, 1.0);
  }

  for(int i = 0; i < 6; ++i){
    is_bc[i] = false;
  }

  eig_mode = monolis_alloc_R_2d(eig_mode, n_node, n_get_eigen);

  monolis_eigen_inverted_standard_lanczos_R(&mat, n_get_eigen, 1.0e-6, 100, eig_val, eig_mode, is_bc);

  monolis_finalize(&mat);
}
