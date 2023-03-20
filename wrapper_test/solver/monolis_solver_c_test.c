#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "monolis.h"
#include "monolis_utils.h"

void monolis_solve_c_test(){
  MONOLIS mat;
  int n_node;
  int n_base;
  int n_dof;
  int n_elem;
  int i, j, k;
  int** elem;
  double val;
  double a[20];
  double b[20];

  monolis_std_log_string("monolis_solve_c_test");

  n_node = 10;
  n_base = 2;
  n_dof  = 2;
  n_elem = 9;

  elem = monolis_alloc_I_2d(elem, n_node, n_base);

  elem[0][0] = 0; elem[0][1] = 1;
  elem[1][0] = 1; elem[1][1] = 2;
  elem[2][0] = 2; elem[2][1] = 3;
  elem[3][0] = 3; elem[3][1] = 4;
  elem[4][0] = 4; elem[4][1] = 5;
  elem[5][0] = 5; elem[5][1] = 6;
  elem[6][0] = 6; elem[6][1] = 7;
  elem[7][0] = 7; elem[7][1] = 8;
  elem[8][0] = 8; elem[8][1] = 9;

  monolis_get_nonzero_pattern_by_simple_mesh_R(&mat, n_node, n_base, n_dof, n_elem, elem);

  for(int i = 0; i < n_node; ++i){
    for(int j = 0; j < n_dof; ++j){
      val = rand()%1001;
      val = val + 2.0;
      monolis_set_scalar_to_sparse_matrix_R(&mat, i, i, j, j, val);
    }
  }

  for(int i = 0; i < n_elem; ++i){
    for(int j = 0; j < n_dof; ++j){
    for(int k = 0; k < n_dof; ++k){
      val = rand()%1001;
      monolis_set_scalar_to_sparse_matrix_R(&mat, elem[i][0], elem[i][1], j, k, val);
      monolis_set_scalar_to_sparse_matrix_R(&mat, elem[i][1], elem[i][0], j, k, val);
    }
    }
  }

  for(int i = 0; i < 20; ++i){
    a[i] = 1.0;
  }

  //monolis_matvec_product_R(&mat, a, b);

  for(int i = 0; i < 20; ++i){
    a[i] = 0.0;
  }

  //monolis_solve_R(&mat, b, a);

  monolis_test_check_eq_R1("monolis_solve_c_test", a[0], 1.0);
}