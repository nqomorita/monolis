#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "monolis.h"
#include "monolis_utils.h"

void monolis_spmat_handler_c_test_R(){
  MONOLIS mat;
  int n_node;
  int n_base;
  int n_dof;
  int n_elem;
  int i;
  int e[1];
  int f[1];
  int** elem;
  double val;
  double b[8];
  double** bval;
  bool is_find;

  monolis_std_log_string("monolis_spmat_handler_c_test_R");

  n_node = 4;
  n_base = 2;
  n_dof  = 2;
  n_elem = 3;

  elem = monolis_alloc_I_2d(elem, n_node, n_base);

  elem[0][0] = 0; elem[0][1] = 1;
  elem[1][0] = 1; elem[1][1] = 2;
  elem[2][0] = 2; elem[2][1] = 3;

  monolis_get_nonzero_pattern_by_simple_mesh_R(&mat, n_node, n_base, n_dof, n_elem, elem);

  val = 1.0;

  monolis_set_scalar_to_sparse_matrix_R(&mat, 1, 1, 0, 0, val);

  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_R 1", mat.mat.R.A[12], 1.0);

  val = 0.0;

  monolis_get_scalar_from_sparse_matrix_R(&mat, 1, 1, 0, 0, &val, &is_find);

  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_R 2", val, 1.0);

  val = 1.0;

  monolis_add_scalar_to_sparse_matrix_R(&mat, 1, 1, 0, 0, val);

  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_R 3", mat.mat.R.A[12], 2.0);

  bval = monolis_alloc_R_2d(bval, 2, 2);

  bval[0][0] = 1.0; bval[0][1] = 2.0;
  bval[1][0] = 3.0; bval[1][1] = 4.0;

  e[0] = 1;

  monolis_add_matrix_to_sparse_matrix_R(&mat, 1, e, bval);

  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_R 4", mat.mat.R.A[12], 3.0);
  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_R 4", mat.mat.R.A[13], 2.0);
  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_R 4", mat.mat.R.A[14], 3.0);
  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_R 4", mat.mat.R.A[15], 4.0);

  e[0] = 1;
  f[0] = 2;

  monolis_add_matrix_to_sparse_matrix_offdiag_R(&mat, 1, 1, e, f, bval);

  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_R 5", mat.mat.R.A[16], 1.0);
  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_R 5", mat.mat.R.A[17], 2.0);
  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_R 5", mat.mat.R.A[18], 3.0);
  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_R 5", mat.mat.R.A[19], 4.0);

  for(i = 0; i < 8; ++i){
    b[i] = 0.0;
  }

  for(i = 0; i < 40; ++i){
    mat.mat.R.A[i] = 2.0;
  }

  monolis_set_Dirichlet_bc_R(&mat, b, 1, 0, 1.0);

  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_R 6", b[0],-2.0);
  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_R 6", b[1],-2.0);
  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_R 6", b[2], 1.0);
  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_R 6", b[3],-2.0);
  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_R 6", b[4],-2.0);
  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_R 6", b[5],-2.0);
  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_R 6", b[6], 0.0);
  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_R 6", b[7], 0.0);
}

void monolis_spmat_handler_c_test_C(){
  MONOLIS mat;
  int n_node;
  int n_base;
  int n_dof;
  int n_elem;
  int i;
  int e[1];
  int f[1];
  int** elem;
  double _Complex val;
  double _Complex b[8];
  double _Complex** bval;
  bool is_find;

  monolis_std_log_string("monolis_spmat_handler_c_test_C");

  n_node = 4;
  n_base = 2;
  n_dof  = 2;
  n_elem = 3;

  elem = monolis_alloc_I_2d(elem, n_node, n_base);

  elem[0][0] = 0; elem[0][1] = 1;
  elem[1][0] = 1; elem[1][1] = 2;
  elem[2][0] = 2; elem[2][1] = 3;

  monolis_get_nonzero_pattern_by_simple_mesh_C(&mat, n_node, n_base, n_dof, n_elem, elem);

  val = 1.0;

  monolis_set_scalar_to_sparse_matrix_C(&mat, 1, 1, 0, 0, val);

  monolis_test_check_eq_C1("monolis_spmat_handler_c_test_C 1", mat.mat.C.A[12], 1.0);

  val = 0.0;

  monolis_get_scalar_from_sparse_matrix_C(&mat, 1, 1, 0, 0, &val, &is_find);

  monolis_test_check_eq_C1("monolis_spmat_handler_c_test_C 2", val, 1.0);

  val = 1.0;

  monolis_add_scalar_to_sparse_matrix_C(&mat, 1, 1, 0, 0, val);

  monolis_test_check_eq_C1("monolis_spmat_handler_c_test_C 3", mat.mat.C.A[12], 2.0);

  bval = monolis_alloc_C_2d(bval, 2, 2);

  bval[0][0] = 1.0; bval[0][1] = 2.0;
  bval[1][0] = 3.0; bval[1][1] = 4.0;

  e[0] = 1;

  monolis_add_matrix_to_sparse_matrix_C(&mat, 1, e, bval);

  monolis_test_check_eq_C1("monolis_spmat_handler_c_test_C 4", mat.mat.C.A[12], 3.0);
  monolis_test_check_eq_C1("monolis_spmat_handler_c_test_C 4", mat.mat.C.A[13], 2.0);
  monolis_test_check_eq_C1("monolis_spmat_handler_c_test_C 4", mat.mat.C.A[14], 3.0);
  monolis_test_check_eq_C1("monolis_spmat_handler_c_test_C 4", mat.mat.C.A[15], 4.0);

  e[0] = 1;
  f[0] = 2;

  monolis_add_matrix_to_sparse_matrix_offdiag_C(&mat, 1, 1, e, f, bval);

  monolis_test_check_eq_C1("monolis_spmat_handler_c_test_C 5", mat.mat.C.A[16], 1.0);
  monolis_test_check_eq_C1("monolis_spmat_handler_c_test_C 5", mat.mat.C.A[17], 2.0);
  monolis_test_check_eq_C1("monolis_spmat_handler_c_test_C 5", mat.mat.C.A[18], 3.0);
  monolis_test_check_eq_C1("monolis_spmat_handler_c_test_C 5", mat.mat.C.A[19], 4.0);

  for(i = 0; i < 8; ++i){
    b[i] = 0.0;
  }

  for(i = 0; i < 40; ++i){
    mat.mat.C.A[i] = 2.0;
  }

  monolis_set_Dirichlet_bc_C(&mat, b, 1, 0, 1.0);

  monolis_test_check_eq_C1("monolis_spmat_handler_c_test_C 6", b[0],-2.0);
  monolis_test_check_eq_C1("monolis_spmat_handler_c_test_C 6", b[1],-2.0);
  monolis_test_check_eq_C1("monolis_spmat_handler_c_test_C 6", b[2], 1.0);
  monolis_test_check_eq_C1("monolis_spmat_handler_c_test_C 6", b[3],-2.0);
  monolis_test_check_eq_C1("monolis_spmat_handler_c_test_C 6", b[4],-2.0);
  monolis_test_check_eq_C1("monolis_spmat_handler_c_test_C 6", b[5],-2.0);
  monolis_test_check_eq_C1("monolis_spmat_handler_c_test_C 6", b[6], 0.0);
  monolis_test_check_eq_C1("monolis_spmat_handler_c_test_C 6", b[7], 0.0);
}

void monolis_spmat_handler_c_test(){
  monolis_spmat_handler_c_test_R();
  monolis_spmat_handler_c_test_C();
}