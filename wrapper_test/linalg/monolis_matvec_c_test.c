#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "monolis.h"
#include "monolis_utils.h"

void monolis_matvec_c_test(){
  MONOLIS mat;
  int n_dof;
  int n_node;
  int n_elem;
  int** elem;
  double rx[4];
  double ry[4];
  double rsum;
  double complex cx[4];
  double complex cy[4];
  double complex csum;

  monolis_std_log_string("monolis_matvec_c_test");

  elem = monolis_alloc_I_2d(elem, 4, 2);

  elem[0][0] = 0; elem[0][1] = 1;
  elem[1][0] = 1; elem[1][1] = 2;
  elem[2][0] = 2; elem[2][1] = 3;
  elem[3][0] = 3; elem[3][1] = 4;

  //monolis_get_nonzero_pattern_by_simple_mesh_R(&mat, n_node, 2, 1, n_elem, elem);

  //monolis_add_scalar_to_sparse_matrix_R(&mat, 0, 0, 0, 0, 2.0);
  //monolis_add_scalar_to_sparse_matrix_R(&mat, 0, 1, 0, 0, 1.0);
  //monolis_add_scalar_to_sparse_matrix_R(&mat, 1, 0, 0, 0, 1.0);
  //monolis_add_scalar_to_sparse_matrix_R(&mat, 1, 1, 0, 0, 2.0);
  //monolis_add_scalar_to_sparse_matrix_R(&mat, 1, 2, 0, 0, 3.0);
  //monolis_add_scalar_to_sparse_matrix_R(&mat, 2, 1, 0, 0, 1.0);
  //monolis_add_scalar_to_sparse_matrix_R(&mat, 2, 2, 0, 0, 2.0);
  //monolis_add_scalar_to_sparse_matrix_R(&mat, 2, 3, 0, 0, 4.0);
  //monolis_add_scalar_to_sparse_matrix_R(&mat, 3, 2, 0, 0, 1.0);
  //monolis_add_scalar_to_sparse_matrix_R(&mat, 3, 3, 0, 0, 2.0);
  //monolis_add_scalar_to_sparse_matrix_R(&mat, 3, 4, 0, 0, 5.0);
  //monolis_add_scalar_to_sparse_matrix_R(&mat, 4, 3, 0, 0, 1.0);
  //monolis_add_scalar_to_sparse_matrix_R(&mat, 4, 4, 0, 0, 2.0);

  //monolis_matvec_product_R(mat, rx, ry);

  //monolis_matvec_product_C(mat, cx, cy);

}