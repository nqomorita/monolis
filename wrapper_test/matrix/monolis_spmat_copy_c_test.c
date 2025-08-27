#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "monolis.h"
#include "monolis_utils.h"

void monolis_spmat_copy_c_test_R()
{
  MONOLIS mat_in;
  MONOLIS mat_out;
  int n_node;
  int n_base;
  int n_dof;
  int n_elem;
  int i, nz;
  int** elem;

  monolis_initialize(&mat_in);
  monolis_initialize(&mat_out);

  n_node = 4;
  n_base = 2;
  n_dof  = 2;
  n_elem = 3;

  elem = monolis_alloc_I_2d(elem, n_node, n_base);

  elem[0][0] = 0; elem[0][1] = 1;
  elem[1][0] = 1; elem[1][1] = 2;
  elem[2][0] = 2; elem[2][1] = 3;

  monolis_get_nonzero_pattern_by_simple_mesh_R(&mat_in, n_node, n_base, n_dof, n_elem, elem);

  nz = mat_in.mat.CSR.index[n_node];
  for (i = 0; i < nz; ++i) {
    mat_in.mat.R.A[i] = 1.0;
  }

  for (i = 0; i < n_node*n_dof; ++i) {
    mat_in.mat.R.B[i] = 2.0;
  }

  for (i = 0; i < n_node*n_dof; ++i) {
    mat_in.mat.R.X[i] = 3.0;
  }

  monolis_copy_mat_R(&mat_in, &mat_out);

  for (i = 0; i < nz; ++i) {
    monolis_test_check_eq_R1("monolis_copy_mat_R a", mat_out.mat.R.A[i], 1.0);
  }

  for (i = 0; i < n_node*n_dof; ++i) {
    monolis_test_check_eq_R1("monolis_copy_mat_R b", mat_out.mat.R.B[i], 2.0);
  }

  for (i = 0; i < n_node*n_dof; ++i) {
    monolis_test_check_eq_R1("monolis_copy_mat_R c", mat_out.mat.R.X[i], 3.0);
  }

  for (i = 0; i < n_node + 1; ++i) {
    monolis_test_check_eq_I1("monolis_copy_mat_R d", mat_in.mat.CSR.index[i], mat_out.mat.CSR.index[i]);
  }

  for (i = 0; i < nz; ++i) {
    monolis_test_check_eq_I1("monolis_copy_mat_R e", mat_in.mat.CSR.item[i], mat_out.mat.CSR.item[i]);
  }

  for (i = 0; i < n_node + 1; ++i) {
    monolis_test_check_eq_I1("monolis_copy_mat_R f", mat_in.mat.CSC.index[i], mat_out.mat.CSC.index[i]);
  }

  for (i = 0; i < nz; ++i) {
    monolis_test_check_eq_I1("monolis_copy_mat_R g", mat_in.mat.CSC.item[i], mat_out.mat.CSC.item[i]);
  }

  for (i = 0; i < nz; ++i) {
    monolis_test_check_eq_I1("monolis_copy_mat_R h", mat_in.mat.CSC.perm[i], mat_out.mat.CSC.perm[i]);
  }

  monolis_clear_mat_value_R(&mat_out);

  for (i = 0; i < nz; ++i) {
    monolis_test_check_eq_R1("monolis_copy_mat_R i", mat_out.mat.R.A[i], 0.0);
  }

  for (i = 0; i < n_node*n_dof; ++i) {
    monolis_test_check_eq_R1("monolis_copy_mat_R j", mat_out.mat.R.B[i], 0.0);
  }

  for (i = 0; i < n_node*n_dof; ++i) {
    monolis_test_check_eq_R1("monolis_copy_mat_R k", mat_out.mat.R.X[i], 0.0);
  }
}

void monolis_spmat_copy_c_test_C()
{
  MONOLIS mat_in;
  MONOLIS mat_out;
  int n_node;
  int n_base;
  int n_dof;
  int n_elem;
  int i, nz;
  int** elem;

  monolis_initialize(&mat_in);
  monolis_initialize(&mat_out);

  n_node = 4;
  n_base = 2;
  n_dof  = 2;
  n_elem = 3;

  elem = monolis_alloc_I_2d(elem, n_node, n_base);

  elem[0][0] = 0; elem[0][1] = 1;
  elem[1][0] = 1; elem[1][1] = 2;
  elem[2][0] = 2; elem[2][1] = 3;

  monolis_get_nonzero_pattern_by_simple_mesh_C(&mat_in, n_node, n_base, n_dof, n_elem, elem);

  nz = mat_in.mat.CSR.index[n_node];
  for (i = 0; i < nz; ++i) {
    mat_in.mat.C.A[i] = 1.0;
  }

  for (i = 0; i < n_node*n_dof; ++i) {
    mat_in.mat.C.B[i] = 2.0;
  }

  for (i = 0; i < n_node*n_dof; ++i) {
    mat_in.mat.C.X[i] = 3.0;
  }

  monolis_copy_mat_C(&mat_in, &mat_out);

  for (i = 0; i < nz; ++i) {
    monolis_test_check_eq_C1("monolis_copy_mat_C a", mat_out.mat.C.A[i], 1.0);
  }

  for (i = 0; i < n_node*n_dof; ++i) {
    monolis_test_check_eq_C1("monolis_copy_mat_C b", mat_out.mat.C.B[i], 2.0);
  }

  for (i = 0; i < n_node*n_dof; ++i) {
    monolis_test_check_eq_C1("monolis_copy_mat_C c", mat_out.mat.C.X[i], 3.0);
  }

  for (i = 0; i < n_node + 1; ++i) {
    monolis_test_check_eq_I1("monolis_copy_mat_C d", mat_in.mat.CSR.index[i], mat_out.mat.CSR.index[i]);
  }

  for (i = 0; i < nz; ++i) {
    monolis_test_check_eq_I1("monolis_copy_mat_C e", mat_in.mat.CSR.item[i], mat_out.mat.CSR.item[i]);
  }

  for (i = 0; i < n_node + 1; ++i) {
    monolis_test_check_eq_I1("monolis_copy_mat_C f", mat_in.mat.CSC.index[i], mat_out.mat.CSC.index[i]);
  }

  for (i = 0; i < nz; ++i) {
    monolis_test_check_eq_I1("monolis_copy_mat_C g", mat_in.mat.CSC.item[i], mat_out.mat.CSC.item[i]);
  }

  for (i = 0; i < nz; ++i) {
    monolis_test_check_eq_I1("monolis_copy_mat_C h", mat_in.mat.CSC.perm[i], mat_out.mat.CSC.perm[i]);
  }

  monolis_clear_mat_value_C(&mat_out);

  for (i = 0; i < nz; ++i) {
    monolis_test_check_eq_C1("monolis_copy_mat_C i", mat_out.mat.C.A[i], 0.0);
  }

  for (i = 0; i < n_node*n_dof; ++i) {
    monolis_test_check_eq_C1("monolis_copy_mat_C j", mat_out.mat.C.B[i], 0.0);
  }

  for (i = 0; i < n_node*n_dof; ++i) {
    monolis_test_check_eq_C1("monolis_copy_mat_C k", mat_out.mat.C.X[i], 0.0);
  }
}

void monolis_spmat_copy_c_test_V_R()
{
  MONOLIS mat_in;
  MONOLIS mat_out;
  int n_node;
  int n_base;
  int* n_dof_list;
  int n_elem;
  int i, nz, total_dof;
  int** elem;

  monolis_initialize(&mat_in);
  monolis_initialize(&mat_out);

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
  monolis_get_nonzero_pattern_by_simple_mesh_V_R(&mat_in, n_node, n_base, n_dof_list, n_elem, elem);

  // テスト用の値を設定
  nz = mat_in.mat.CSR.index[n_node];
  for (i = 0; i < nz; ++i) {
    mat_in.mat.R.A[i] = 1.0;
  }

  for (i = 0; i < total_dof; ++i) {
    mat_in.mat.R.B[i] = 2.0;
  }

  for (i = 0; i < total_dof; ++i) {
    mat_in.mat.R.X[i] = 3.0;
  }

  // 行列をコピー
  monolis_copy_mat_R(&mat_in, &mat_out);

  // 基本パラメータのテスト
  monolis_test_check_eq_I1("monolis_copy_mat_V_R N", mat_in.mat.N, mat_out.mat.N);
  monolis_test_check_eq_I1("monolis_copy_mat_V_R NP", mat_in.mat.NP, mat_out.mat.NP);
  monolis_test_check_eq_I1("monolis_copy_mat_V_R NDOF", mat_in.mat.NDOF, mat_out.mat.NDOF);

  // 行列値のテスト
  for (i = 0; i < nz; ++i) {
    monolis_test_check_eq_R1("monolis_copy_mat_V_R A", mat_out.mat.R.A[i], 1.0);
  }

  for (i = 0; i < total_dof; ++i) {
    monolis_test_check_eq_R1("monolis_copy_mat_V_R B", mat_out.mat.R.B[i], 2.0);
  }

  for (i = 0; i < total_dof; ++i) {
    monolis_test_check_eq_R1("monolis_copy_mat_V_R X", mat_out.mat.R.X[i], 3.0);
  }

  // CSR構造のテスト
  for (i = 0; i < n_node + 1; ++i) {
    monolis_test_check_eq_I1("monolis_copy_mat_V_R CSR index", mat_in.mat.CSR.index[i], mat_out.mat.CSR.index[i]);
  }

  for (i = 0; i < nz; ++i) {
    monolis_test_check_eq_I1("monolis_copy_mat_V_R CSR item", mat_in.mat.CSR.item[i], mat_out.mat.CSR.item[i]);
  }

  // CSC構造のテスト
  for (i = 0; i < n_node + 1; ++i) {
    monolis_test_check_eq_I1("monolis_copy_mat_V_R CSC index", mat_in.mat.CSC.index[i], mat_out.mat.CSC.index[i]);
  }

  for (i = 0; i < nz; ++i) {
    monolis_test_check_eq_I1("monolis_copy_mat_V_R CSC item", mat_in.mat.CSC.item[i], mat_out.mat.CSC.item[i]);
  }

  for (i = 0; i < nz; ++i) {
    monolis_test_check_eq_I1("monolis_copy_mat_V_R CSC perm", mat_in.mat.CSC.perm[i], mat_out.mat.CSC.perm[i]);
  }

  // 自由度インデックスのテスト
  for (i = 0; i < n_node + 1; ++i) {
    monolis_test_check_eq_I1("monolis_copy_mat_V_R n_dof_index", mat_in.mat.n_dof_index[i], mat_out.mat.n_dof_index[i]);
  }

  for (i = 0; i < n_node + 1; ++i) {
    monolis_test_check_eq_I1("monolis_copy_mat_V_R n_dof_index2", mat_in.mat.n_dof_index2[i], mat_out.mat.n_dof_index2[i]);
  }

  for (i = 0; i < n_node; ++i) {
    monolis_test_check_eq_I1("monolis_copy_mat_V_R n_dof_list", mat_in.mat.n_dof_list[i], mat_out.mat.n_dof_list[i]);
  }

  // 値をクリアしてテスト
  monolis_clear_mat_value_R(&mat_out);

  for (i = 0; i < nz; ++i) {
    monolis_test_check_eq_R1("monolis_copy_mat_V_R clear A", mat_out.mat.R.A[i], 0.0);
  }

  for (i = 0; i < total_dof; ++i) {
    monolis_test_check_eq_R1("monolis_copy_mat_V_R clear B", mat_out.mat.R.B[i], 0.0);
  }

  for (i = 0; i < total_dof; ++i) {
    monolis_test_check_eq_R1("monolis_copy_mat_V_R clear X", mat_out.mat.R.X[i], 0.0);
  }
}

void monolis_spmat_copy_c_test()
{
  monolis_std_log_string("monolis_spmat_copy_c_test");
  monolis_spmat_copy_c_test_R();
  monolis_spmat_copy_c_test_C();
  monolis_spmat_copy_c_test_V_R();
}
