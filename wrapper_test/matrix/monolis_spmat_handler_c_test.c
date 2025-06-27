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

void monolis_spmat_handler_c_test_V_R(){
  MONOLIS mat;
  int n_node;
  int n_base;
  int* n_dof_list;
  int n_elem;
  int i;
  int e[1];
  int f[1];
  int** elem;
  double val;
  double* b;
  double** bval;
  bool is_find;
  int total_dof;

  monolis_std_log_string("monolis_spmat_handler_c_test_V_R");

  // 4節点、2要素のメッシュ設定
  n_node = 4;
  n_base = 3;
  n_elem = 2;

  // 各節点の自由度を異なる値に設定
  // 節点1: 1自由度, 節点2: 2自由度, 節点3: 3自由度, 節点4: 1自由度
  n_dof_list = monolis_alloc_I_1d(n_dof_list, n_node);
  n_dof_list[0] = 1;
  n_dof_list[1] = 2;
  n_dof_list[2] = 3;
  n_dof_list[3] = 1;

  // 総自由度数を計算
  total_dof = 0;
  for (i = 0; i < n_node; ++i) {
    total_dof += n_dof_list[i];
  }

  // 三角形要素の設定
  elem = monolis_alloc_I_2d(elem, n_elem, n_base);
  elem[0][0] = 0; elem[0][1] = 1; elem[0][2] = 2;
  elem[1][0] = 1; elem[1][1] = 2; elem[1][2] = 3;

  // monolis_get_nonzero_pattern_by_simple_mesh_V_Rで非ゼロ構造を決定
  monolis_get_nonzero_pattern_by_simple_mesh_V_R(&mat, n_node, n_base, n_dof_list, n_elem, elem);

  // スカラー値の設定・取得・加算テスト
  val = 1.0;
  monolis_set_scalar_to_sparse_matrix_R(&mat, 0, 0, 0, 0, val);
  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_V_R 1", mat.mat.R.A[0], 1.0);

  val = 0.0;
  monolis_get_scalar_from_sparse_matrix_R(&mat, 0, 0, 0, 0, &val, &is_find);
  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_V_R 2", val, 1.0);

  val = 1.0;
  monolis_add_scalar_to_sparse_matrix_R(&mat, 0, 0, 0, 0, val);
  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_V_R 3", mat.mat.R.A[0], 2.0);

  val = 10.0;
  monolis_add_scalar_to_sparse_matrix_R(&mat, 0, 1, 0, 0, val);
  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_V_R 3a", mat.mat.R.A[1], 10.0);

  val = 11.0;
  monolis_add_scalar_to_sparse_matrix_R(&mat, 0, 1, 0, 1, val);
  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_V_R 3b", mat.mat.R.A[2], 11.0);

  // 節点2（2自由度）に対するテスト
  val = 3.0;
  monolis_set_scalar_to_sparse_matrix_R(&mat, 1, 1, 0, 0, val);
  val = 4.0;
  monolis_set_scalar_to_sparse_matrix_R(&mat, 1, 1, 1, 1, val);
  
  val = 0.0;
  monolis_get_scalar_from_sparse_matrix_R(&mat, 1, 1, 0, 0, &val, &is_find);
  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_V_R 4", val, 3.0);
  
  val = 0.0;
  monolis_get_scalar_from_sparse_matrix_R(&mat, 1, 1, 1, 1, &val, &is_find);
  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_V_R 5", val, 4.0);

  val = 4.0;
  monolis_add_scalar_to_sparse_matrix_R(&mat, 1, 1, 1, 1, val);

  val = 0.0;
  monolis_get_scalar_from_sparse_matrix_R(&mat, 1, 1, 1, 1, &val, &is_find);
  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_V_R 5", val, 8.0);

  // 節点3（3自由度）に対するテスト
  val = 5.0;
  monolis_set_scalar_to_sparse_matrix_R(&mat, 2, 2, 2, 2, val);
  
  val = 0.0;
  monolis_get_scalar_from_sparse_matrix_R(&mat, 2, 2, 2, 2, &val, &is_find);
  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_V_R 6", val, 5.0);

  // 非対角成分のテスト（節点1-節点2間）
  val = 6.0;
  monolis_set_scalar_to_sparse_matrix_R(&mat, 0, 1, 0, 0, val);
  
  val = 0.0;
  monolis_get_scalar_from_sparse_matrix_R(&mat, 0, 1, 0, 0, &val, &is_find);
  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_V_R 7", val, 6.0);

  // 小行列の設定テスト（節点1の1x1行列）
  bval = monolis_alloc_R_2d(bval, 1, 1);
  bval[0][0] = 7.0;
  e[0] = 0;
  monolis_add_matrix_to_sparse_matrix_R(&mat, 1, e, bval);
  
  val = 0.0;
  monolis_get_scalar_from_sparse_matrix_R(&mat, 0, 0, 0, 0, &val, &is_find);
  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_V_R 8", val, 9.0); // 2.0 + 7.0

  // 2x2行列のテスト（節点2）
  bval = monolis_alloc_R_2d(bval, 2, 2);
  bval[0][0] = 1.0; bval[0][1] = 2.0;
  bval[1][0] = 3.0; bval[1][1] = 4.0;
  e[0] = 1;
  monolis_add_matrix_to_sparse_matrix_R(&mat, 1, e, bval);

  val = 0.0;
  monolis_get_scalar_from_sparse_matrix_R(&mat, 1, 1, 0, 0, &val, &is_find);
  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_V_R 9", val, 4.0); // 3.0 + 1.0

  val = 0.0;
  monolis_get_scalar_from_sparse_matrix_R(&mat, 1, 1, 1, 1, &val, &is_find);
  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_V_R 10", val, 12.0); // 4.0 + 4.0 + 4.0

  // Dirichlet境界条件のテスト
  b = monolis_alloc_R_1d(b, total_dof);
  for(i = 0; i < total_dof; ++i){
    b[i] = 0.0;
  }

  // 行列を一様な値で初期化
  for(i = 0; i < mat.mat.CSR.index[n_node]; ++i){
    mat.mat.R.A[i] = 2.0;
  }

  // 節点1の第1自由度に境界条件を設定
  monolis_set_Dirichlet_bc_R(&mat, b, 0, 0, 10.0);

  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_V_R 11", b[0], 10.0);
  
  // 節点2の第2自由度に境界条件を設定
  monolis_set_Dirichlet_bc_R(&mat, b, 1, 1, 20.0);
  
  monolis_test_check_eq_R1("monolis_spmat_handler_c_test_V_R 12", b[2], 20.0);
}

void monolis_spmat_handler_c_test(){
  monolis_spmat_handler_c_test_R();
  monolis_spmat_handler_c_test_C();
  monolis_spmat_handler_c_test_V_R();
}