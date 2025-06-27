#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "monolis.h"
#include "monolis_utils.h"

void monolis_spmat_nzpattern_c_test_R(){
  MONOLIS mat;
  int n_node;
  int n_base;
  int n_dof;
  int n_elem;
  int** elem;

  monolis_std_log_string("monolis_spmat_nzpattern_c_test_R");

  n_node = 4;
  n_base = 2;
  n_dof  = 2;
  n_elem = 3;

  elem = monolis_alloc_I_2d(elem, n_node, n_base);

  elem[0][0] = 0; elem[0][1] = 1;
  elem[1][0] = 1; elem[1][1] = 2;
  elem[2][0] = 2; elem[2][1] = 3;

  monolis_get_nonzero_pattern_by_simple_mesh_R(&mat, n_node, n_base, n_dof, n_elem, elem);

  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 1",  mat.mat.N, 4);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 2",  mat.mat.NP, 4);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 3",  mat.mat.NDOF, 2);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 4",  mat.mat.CSR.index[0], 0);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 5",  mat.mat.CSR.index[1], 2);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 6",  mat.mat.CSR.index[2], 5);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 7",  mat.mat.CSR.index[3], 8);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 8",  mat.mat.CSR.index[4], 10);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 9",  mat.mat.CSR.item[0], 1);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 10", mat.mat.CSR.item[1], 2);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 11", mat.mat.CSR.item[2], 1);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 12", mat.mat.CSR.item[3], 2);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 13", mat.mat.CSR.item[4], 3);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 14", mat.mat.CSR.item[5], 2);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 15", mat.mat.CSR.item[6], 3);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 16", mat.mat.CSR.item[7], 4);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 17", mat.mat.CSR.item[8], 3);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 18", mat.mat.CSR.item[9], 4);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 19", mat.mat.CSC.index[0], 0);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 20", mat.mat.CSC.index[1], 2);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 21", mat.mat.CSC.index[2], 5);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 22", mat.mat.CSC.index[3], 8);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 23", mat.mat.CSC.index[4], 10);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 24", mat.mat.CSC.item[0], 1);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 25", mat.mat.CSC.item[1], 2);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 26", mat.mat.CSC.item[2], 1);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 27", mat.mat.CSC.item[3], 2);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 28", mat.mat.CSC.item[4], 3);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 29", mat.mat.CSC.item[5], 2);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 30", mat.mat.CSC.item[6], 3);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 31", mat.mat.CSC.item[7], 4);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 32", mat.mat.CSC.item[8], 3);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 33", mat.mat.CSC.item[9], 4);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 34", mat.mat.CSC.perm[0], 1);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 35", mat.mat.CSC.perm[1], 3);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 36", mat.mat.CSC.perm[2], 2);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 37", mat.mat.CSC.perm[3], 4);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 38", mat.mat.CSC.perm[4], 6);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 39", mat.mat.CSC.perm[5], 5);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 40", mat.mat.CSC.perm[6], 7);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 41", mat.mat.CSC.perm[7], 9);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 42", mat.mat.CSC.perm[8], 8);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_R_test 43", mat.mat.CSC.perm[9], 10);
}

void monolis_spmat_nzpattern_c_test_C(){
  MONOLIS mat;
  int n_node;
  int n_base;
  int n_dof;
  int n_elem;
  int** elem;

  monolis_std_log_string("monolis_spmat_nzpattern_c_test_C");

  n_node = 4;
  n_base = 2;
  n_dof  = 2;
  n_elem = 3;

  elem = monolis_alloc_I_2d(elem, n_node, n_base);

  elem[0][0] = 0; elem[0][1] = 1;
  elem[1][0] = 1; elem[1][1] = 2;
  elem[2][0] = 2; elem[2][1] = 3;

  monolis_get_nonzero_pattern_by_simple_mesh_C(&mat, n_node, n_base, n_dof, n_elem, elem);

  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 1",  mat.mat.N, 4);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 2",  mat.mat.NP, 4);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 3",  mat.mat.NDOF, 2);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 4",  mat.mat.CSR.index[0], 0);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 5",  mat.mat.CSR.index[1], 2);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 6",  mat.mat.CSR.index[2], 5);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 7",  mat.mat.CSR.index[3], 8);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 8",  mat.mat.CSR.index[4], 10);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 9",  mat.mat.CSR.item[0], 1);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 10", mat.mat.CSR.item[1], 2);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 11", mat.mat.CSR.item[2], 1);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 12", mat.mat.CSR.item[3], 2);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 13", mat.mat.CSR.item[4], 3);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 14", mat.mat.CSR.item[5], 2);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 15", mat.mat.CSR.item[6], 3);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 16", mat.mat.CSR.item[7], 4);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 17", mat.mat.CSR.item[8], 3);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 18", mat.mat.CSR.item[9], 4);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 19", mat.mat.CSC.index[0], 0);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 20", mat.mat.CSC.index[1], 2);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 21", mat.mat.CSC.index[2], 5);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 22", mat.mat.CSC.index[3], 8);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 23", mat.mat.CSC.index[4], 10);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 24", mat.mat.CSC.item[0], 1);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 25", mat.mat.CSC.item[1], 2);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 26", mat.mat.CSC.item[2], 1);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 27", mat.mat.CSC.item[3], 2);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 28", mat.mat.CSC.item[4], 3);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 29", mat.mat.CSC.item[5], 2);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 30", mat.mat.CSC.item[6], 3);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 31", mat.mat.CSC.item[7], 4);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 32", mat.mat.CSC.item[8], 3);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 33", mat.mat.CSC.item[9], 4);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 34", mat.mat.CSC.perm[0], 1);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 35", mat.mat.CSC.perm[1], 3);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 36", mat.mat.CSC.perm[2], 2);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 37", mat.mat.CSC.perm[3], 4);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 38", mat.mat.CSC.perm[4], 6);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 39", mat.mat.CSC.perm[5], 5);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 40", mat.mat.CSC.perm[6], 7);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 41", mat.mat.CSC.perm[7], 9);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 42", mat.mat.CSC.perm[8], 8);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_C_test 43", mat.mat.CSC.perm[9], 10);
}

void monolis_spmat_nzpattern_c_test_V_R(){
  MONOLIS mat;
  int n_node;
  int n_base;
  int* n_dof_list;
  int n_elem;
  int** elem;
  int total_dof;
  int i;

  monolis_std_log_string("monolis_spmat_nzpattern_c_test_V_R");

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

  // 基本パラメータのテスト
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 1", mat.mat.N, 4);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 2", mat.mat.NP, 4);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 3", mat.mat.NDOF, -1); // 任意自由度の場合は-1

  // 自由度インデックスのテスト
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 4", mat.mat.n_dof_index[0], 0);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 5", mat.mat.n_dof_index[1], 1);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 6", mat.mat.n_dof_index[2], 3);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 7", mat.mat.n_dof_index[3], 6);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 8", mat.mat.n_dof_index[4], 7);

  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 9", mat.mat.n_dof_index2[0], 0);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 10", mat.mat.n_dof_index2[1], 1);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 11", mat.mat.n_dof_index2[2], 3);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 12", mat.mat.n_dof_index2[3], 6);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 13", mat.mat.n_dof_index2[4], 8);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 13", mat.mat.n_dof_index2[5], 12);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 13", mat.mat.n_dof_index2[6], 18);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 13", mat.mat.n_dof_index2[7], 20);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 13", mat.mat.n_dof_index2[8], 23);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 13", mat.mat.n_dof_index2[9], 29);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 13", mat.mat.n_dof_index2[10], 38);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 13", mat.mat.n_dof_index2[11], 41);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 13", mat.mat.n_dof_index2[12], 43);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 13", mat.mat.n_dof_index2[13], 46);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 13", mat.mat.n_dof_index2[14], 47);

  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 14", mat.mat.n_dof_list[0], 1);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 15", mat.mat.n_dof_list[1], 2);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 16", mat.mat.n_dof_list[2], 3);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 17", mat.mat.n_dof_list[3], 1);


  // CSR構造のテスト - 節点グラフ構造を確認
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 9",  mat.mat.CSR.index[0], 0);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 10", mat.mat.CSR.index[1], 3);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 11", mat.mat.CSR.index[2], 7);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 12", mat.mat.CSR.index[3], 11);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 13", mat.mat.CSR.index[4], 14);

  // CSR item配列のテスト
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 14", mat.mat.CSR.item[0], 1);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 15", mat.mat.CSR.item[1], 2);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 16", mat.mat.CSR.item[2], 3);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 17", mat.mat.CSR.item[3], 1);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 18", mat.mat.CSR.item[4], 2);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 19", mat.mat.CSR.item[5], 3);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 20", mat.mat.CSR.item[6], 4);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 21", mat.mat.CSR.item[7], 1);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 22", mat.mat.CSR.item[8], 2);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 23", mat.mat.CSR.item[9], 3);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 24", mat.mat.CSR.item[10], 4);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 25", mat.mat.CSR.item[11], 2);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 26", mat.mat.CSR.item[12], 3);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 27", mat.mat.CSR.item[13], 4);

  // CSC構造のテスト
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 24", mat.mat.CSC.index[0], 0);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 25", mat.mat.CSC.index[1], 3);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 26", mat.mat.CSC.index[2], 7);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 27", mat.mat.CSC.index[3], 11);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 28", mat.mat.CSC.index[4], 14);

  // CSC item配列のテスト
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 29", mat.mat.CSC.item[0], 1);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 30", mat.mat.CSC.item[1], 2);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 31", mat.mat.CSC.item[2], 3);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 32", mat.mat.CSC.item[3], 1);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 33", mat.mat.CSC.item[4], 2);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 34", mat.mat.CSC.item[5], 3);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 35", mat.mat.CSC.item[6], 4);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 36", mat.mat.CSC.item[7], 1);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 37", mat.mat.CSC.item[8], 2);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 38", mat.mat.CSC.item[9], 3);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 38", mat.mat.CSC.item[10], 4);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 38", mat.mat.CSC.item[11], 2);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 38", mat.mat.CSC.item[12], 3);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 38", mat.mat.CSC.item[13], 4);

  // CSC perm配列のテスト
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 39", mat.mat.CSC.perm[0], 1);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 40", mat.mat.CSC.perm[1], 4);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 41", mat.mat.CSC.perm[2], 8);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 42", mat.mat.CSC.perm[3], 2);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 43", mat.mat.CSC.perm[4], 5);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 44", mat.mat.CSC.perm[5], 9);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 45", mat.mat.CSC.perm[6], 12);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 46", mat.mat.CSC.perm[7], 3);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 47", mat.mat.CSC.perm[8], 6);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 48", mat.mat.CSC.perm[9], 10);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 48", mat.mat.CSC.perm[10], 13);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 48", mat.mat.CSC.perm[11], 7);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 48", mat.mat.CSC.perm[12], 11);
  monolis_test_check_eq_I1("monolis_get_nzpattern_by_simple_mesh_V_R_test 48", mat.mat.CSC.perm[13], 14);
}

void monolis_spmat_nzpattern_c_test(){
  monolis_spmat_nzpattern_c_test_R();
  monolis_spmat_nzpattern_c_test_C();
  monolis_spmat_nzpattern_c_test_V_R();
}
