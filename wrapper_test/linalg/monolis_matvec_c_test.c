#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "monolis.h"
#include "monolis_utils.h"

void monolis_matvec_c_test_R(){
  MONOLIS mat;
  int n_dof;
  int n_node;
  int n_elem;
  int** elem;
  double rx[5];
  double ry[5];

  monolis_std_log_string("monolis_matvec_c_test_R");

  elem = monolis_alloc_I_2d(elem, 4, 2);

  n_dof = 1;

  n_node = 5;

  n_elem = 4;

  elem[0][0] = 0; elem[0][1] = 1;
  elem[1][0] = 1; elem[1][1] = 2;
  elem[2][0] = 2; elem[2][1] = 3;
  elem[3][0] = 3; elem[3][1] = 4;

  monolis_initialize(&mat, "./");

  monolis_get_nonzero_pattern_by_simple_mesh_R(&mat, n_node, 2, n_dof, n_elem, elem);

  monolis_add_scalar_to_sparse_matrix_R(&mat, 0, 0, 0, 0, 2.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 0, 1, 0, 0, 1.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 1, 0, 0, 0, 1.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 1, 1, 0, 0, 2.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 1, 2, 0, 0, 3.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 2, 1, 0, 0, 1.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 2, 2, 0, 0, 2.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 2, 3, 0, 0, 4.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 3, 2, 0, 0, 1.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 3, 3, 0, 0, 2.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 3, 4, 0, 0, 5.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 4, 3, 0, 0, 1.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 4, 4, 0, 0, 2.0);

  for(int i = 0; i < 5; ++i){
    rx[i] = 1.0;
  }

  for(int i = 0; i < 5; ++i){
    ry[i] = 0.0;
  }

  monolis_matvec_product_R(&mat, rx, ry);

  //monolis_finalize(&mat);
}

void monolis_matvec_c_test(){
  monolis_matvec_c_test_R();
}