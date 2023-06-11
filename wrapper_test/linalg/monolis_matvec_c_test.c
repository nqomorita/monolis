#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "monolis.h"
#include "monolis_utils.h"

void monolis_matvec_c_test_R11(){
  MONOLIS mat;
  MONOLIS_COM com;
  int i;
  int n_dof;
  int n_node;
  int n_elem;
  int** elem;
  double rx[5];
  double ry[5];

  monolis_std_log_string("monolis_matvec_c_test_R11");

  elem = monolis_alloc_I_2d(elem, 4, 2);

  n_dof = 1;

  n_node = 5;

  n_elem = 4;

  elem[0][0] = 0; elem[0][1] = 1;
  elem[1][0] = 1; elem[1][1] = 2;
  elem[2][0] = 2; elem[2][1] = 3;
  elem[3][0] = 3; elem[3][1] = 4;

  monolis_initialize(&mat);
  monolis_com_initialize_by_self(&com);

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

  for(i = 0; i < 5; ++i){
    rx[i] = 1.0;
  }

  for(i = 0; i < 5; ++i){
    ry[i] = 0.0;
  }

  monolis_matvec_product_R(&mat, &com, rx, ry);

  monolis_test_check_eq_R1("monolis_matvec_c_test_R11", ry[0], 3.0);
  monolis_test_check_eq_R1("monolis_matvec_c_test_R11", ry[1], 6.0);
  monolis_test_check_eq_R1("monolis_matvec_c_test_R11", ry[2], 7.0);
  monolis_test_check_eq_R1("monolis_matvec_c_test_R11", ry[3], 8.0);
  monolis_test_check_eq_R1("monolis_matvec_c_test_R11", ry[4], 3.0);

  monolis_finalize(&mat);
}

void monolis_matvec_c_test_R22(){
  MONOLIS mat;
  MONOLIS_COM com;
  int i;
  int n_dof;
  int n_node;
  int n_elem;
  int** elem;
  double rx[10];
  double ry[10];

  monolis_std_log_string("monolis_matvec_c_test_R22");

  elem = monolis_alloc_I_2d(elem, 4, 2);

  n_dof = 2;

  n_node = 5;

  n_elem = 4;

  elem[0][0] = 0; elem[0][1] = 1;
  elem[1][0] = 1; elem[1][1] = 2;
  elem[2][0] = 2; elem[2][1] = 3;
  elem[3][0] = 3; elem[3][1] = 4;

  monolis_initialize(&mat);
  monolis_com_initialize_by_self(&com);

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

  monolis_add_scalar_to_sparse_matrix_R(&mat, 0, 0, 1, 1, 4.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 0, 1, 1, 1, 2.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 1, 0, 1, 1, 2.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 1, 1, 1, 1, 4.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 1, 2, 1, 1, 6.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 2, 1, 1, 1, 2.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 2, 2, 1, 1, 4.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 2, 3, 1, 1, 8.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 3, 2, 1, 1, 2.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 3, 3, 1, 1, 4.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 3, 4, 1, 1,10.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 4, 3, 1, 1, 2.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 4, 4, 1, 1, 4.0);

  for(i = 0; i < 10; ++i){
    rx[i] = 1.0;
  }

  for(i = 0; i < 10; ++i){
    ry[i] = 0.0;
  }

  monolis_matvec_product_R(&mat, &com, rx, ry);

  monolis_test_check_eq_R1("monolis_matvec_c_test_R22", ry[0], 3.0);
  monolis_test_check_eq_R1("monolis_matvec_c_test_R22", ry[1], 6.0);
  monolis_test_check_eq_R1("monolis_matvec_c_test_R22", ry[2], 6.0);
  monolis_test_check_eq_R1("monolis_matvec_c_test_R22", ry[3],12.0);
  monolis_test_check_eq_R1("monolis_matvec_c_test_R22", ry[4], 7.0);
  monolis_test_check_eq_R1("monolis_matvec_c_test_R22", ry[5],14.0);
  monolis_test_check_eq_R1("monolis_matvec_c_test_R22", ry[6], 8.0);
  monolis_test_check_eq_R1("monolis_matvec_c_test_R22", ry[7],16.0);
  monolis_test_check_eq_R1("monolis_matvec_c_test_R22", ry[8], 3.0);
  monolis_test_check_eq_R1("monolis_matvec_c_test_R22", ry[9], 6.0);

  monolis_finalize(&mat);
}

void monolis_matvec_c_test_R33(){
  MONOLIS mat;
  MONOLIS_COM com;
  int i;
  int n_dof;
  int n_node;
  int n_elem;
  int** elem;
  double rx[15];
  double ry[15];

  monolis_std_log_string("monolis_matvec_c_test_R22");

  elem = monolis_alloc_I_2d(elem, 4, 2);

  n_dof = 3;

  n_node = 5;

  n_elem = 4;

  elem[0][0] = 0; elem[0][1] = 1;
  elem[1][0] = 1; elem[1][1] = 2;
  elem[2][0] = 2; elem[2][1] = 3;
  elem[3][0] = 3; elem[3][1] = 4;

  monolis_initialize(&mat);
  monolis_com_initialize_by_self(&com);

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

  monolis_add_scalar_to_sparse_matrix_R(&mat, 0, 0, 1, 1, 4.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 0, 1, 1, 1, 2.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 1, 0, 1, 1, 2.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 1, 1, 1, 1, 4.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 1, 2, 1, 1, 6.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 2, 1, 1, 1, 2.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 2, 2, 1, 1, 4.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 2, 3, 1, 1, 8.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 3, 2, 1, 1, 2.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 3, 3, 1, 1, 4.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 3, 4, 1, 1,10.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 4, 3, 1, 1, 2.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 4, 4, 1, 1, 4.0);

  monolis_add_scalar_to_sparse_matrix_R(&mat, 0, 0, 2, 2, 8.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 0, 1, 2, 2, 4.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 1, 0, 2, 2, 4.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 1, 1, 2, 2, 8.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 1, 2, 2, 2,12.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 2, 1, 2, 2, 4.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 2, 2, 2, 2, 8.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 2, 3, 2, 2,16.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 3, 2, 2, 2, 4.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 3, 3, 2, 2, 8.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 3, 4, 2, 2,20.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 4, 3, 2, 2, 4.0);
  monolis_add_scalar_to_sparse_matrix_R(&mat, 4, 4, 2, 2, 8.0);

  for(i = 0; i < 15; ++i){
    rx[i] = 1.0;
  }

  for(i = 0; i < 15; ++i){
    ry[i] = 0.0;
  }

  monolis_matvec_product_R(&mat, &com, rx, ry);

  monolis_test_check_eq_R1("monolis_matvec_c_test_R33", ry[0], 3.0);
  monolis_test_check_eq_R1("monolis_matvec_c_test_R33", ry[1], 6.0);
  monolis_test_check_eq_R1("monolis_matvec_c_test_R33", ry[2],12.0);
  monolis_test_check_eq_R1("monolis_matvec_c_test_R33", ry[3], 6.0);
  monolis_test_check_eq_R1("monolis_matvec_c_test_R33", ry[4],12.0);
  monolis_test_check_eq_R1("monolis_matvec_c_test_R33", ry[5],24.0);
  monolis_test_check_eq_R1("monolis_matvec_c_test_R33", ry[6], 7.0);
  monolis_test_check_eq_R1("monolis_matvec_c_test_R33", ry[7],14.0);
  monolis_test_check_eq_R1("monolis_matvec_c_test_R33", ry[8],28.0);
  monolis_test_check_eq_R1("monolis_matvec_c_test_R33", ry[9], 8.0);
  monolis_test_check_eq_R1("monolis_matvec_c_test_R33", ry[10],16.0);
  monolis_test_check_eq_R1("monolis_matvec_c_test_R33", ry[11],32.0);
  monolis_test_check_eq_R1("monolis_matvec_c_test_R33", ry[12], 3.0);
  monolis_test_check_eq_R1("monolis_matvec_c_test_R33", ry[13], 6.0);
  monolis_test_check_eq_R1("monolis_matvec_c_test_R33", ry[14],12.0);

  monolis_finalize(&mat);
}

void monolis_matvec_c_test_C11(){
  MONOLIS mat;
  MONOLIS_COM com;
  int i;
  int n_dof;
  int n_node;
  int n_elem;
  int** elem;
  double _Complex rx[5];
  double _Complex ry[5];

  monolis_std_log_string("monolis_matvec_c_test_C11");

  elem = monolis_alloc_I_2d(elem, 4, 2);

  n_dof = 1;

  n_node = 5;

  n_elem = 4;

  elem[0][0] = 0; elem[0][1] = 1;
  elem[1][0] = 1; elem[1][1] = 2;
  elem[2][0] = 2; elem[2][1] = 3;
  elem[3][0] = 3; elem[3][1] = 4;

  monolis_initialize(&mat);
  monolis_com_initialize_by_self(&com);

  monolis_get_nonzero_pattern_by_simple_mesh_C(&mat, n_node, 2, n_dof, n_elem, elem);

  monolis_add_scalar_to_sparse_matrix_C(&mat, 0, 0, 0, 0, 2.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 0, 1, 0, 0, 1.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 1, 0, 0, 0, 1.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 1, 1, 0, 0, 2.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 1, 2, 0, 0, 3.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 2, 1, 0, 0, 1.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 2, 2, 0, 0, 2.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 2, 3, 0, 0, 4.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 3, 2, 0, 0, 1.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 3, 3, 0, 0, 2.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 3, 4, 0, 0, 5.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 4, 3, 0, 0, 1.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 4, 4, 0, 0, 2.0);

  for(i = 0; i < 5; ++i){
    rx[i] = 1.0;
  }

  for(i = 0; i < 5; ++i){
    ry[i] = 0.0;
  }

  monolis_matvec_product_C(&mat, &com, rx, ry);

  monolis_test_check_eq_C1("monolis_matvec_c_test_C11", ry[0], 3.0);
  monolis_test_check_eq_C1("monolis_matvec_c_test_C11", ry[1], 6.0);
  monolis_test_check_eq_C1("monolis_matvec_c_test_C11", ry[2], 7.0);
  monolis_test_check_eq_C1("monolis_matvec_c_test_C11", ry[3], 8.0);
  monolis_test_check_eq_C1("monolis_matvec_c_test_C11", ry[4], 3.0);

  monolis_finalize(&mat);
}

void monolis_matvec_c_test_C22(){
  MONOLIS mat;
  MONOLIS_COM com;
  int i;
  int n_dof;
  int n_node;
  int n_elem;
  int** elem;
  double _Complex rx[10];
  double _Complex ry[10];

  monolis_std_log_string("monolis_matvec_c_test_C22");

  elem = monolis_alloc_I_2d(elem, 4, 2);

  n_dof = 2;

  n_node = 5;

  n_elem = 4;

  elem[0][0] = 0; elem[0][1] = 1;
  elem[1][0] = 1; elem[1][1] = 2;
  elem[2][0] = 2; elem[2][1] = 3;
  elem[3][0] = 3; elem[3][1] = 4;

  monolis_initialize(&mat);
  monolis_com_initialize_by_self(&com);

  monolis_get_nonzero_pattern_by_simple_mesh_C(&mat, n_node, 2, n_dof, n_elem, elem);

  monolis_add_scalar_to_sparse_matrix_C(&mat, 0, 0, 0, 0, 2.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 0, 1, 0, 0, 1.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 1, 0, 0, 0, 1.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 1, 1, 0, 0, 2.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 1, 2, 0, 0, 3.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 2, 1, 0, 0, 1.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 2, 2, 0, 0, 2.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 2, 3, 0, 0, 4.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 3, 2, 0, 0, 1.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 3, 3, 0, 0, 2.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 3, 4, 0, 0, 5.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 4, 3, 0, 0, 1.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 4, 4, 0, 0, 2.0);

  monolis_add_scalar_to_sparse_matrix_C(&mat, 0, 0, 1, 1, 4.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 0, 1, 1, 1, 2.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 1, 0, 1, 1, 2.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 1, 1, 1, 1, 4.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 1, 2, 1, 1, 6.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 2, 1, 1, 1, 2.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 2, 2, 1, 1, 4.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 2, 3, 1, 1, 8.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 3, 2, 1, 1, 2.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 3, 3, 1, 1, 4.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 3, 4, 1, 1,10.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 4, 3, 1, 1, 2.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 4, 4, 1, 1, 4.0);

  for(i = 0; i < 10; ++i){
    rx[i] = 1.0;
  }

  for(i = 0; i < 10; ++i){
    ry[i] = 0.0;
  }

  monolis_matvec_product_C(&mat, &com, rx, ry);

  monolis_test_check_eq_C1("monolis_matvec_c_test_C22", ry[0], 3.0);
  monolis_test_check_eq_C1("monolis_matvec_c_test_C22", ry[1], 6.0);
  monolis_test_check_eq_C1("monolis_matvec_c_test_C22", ry[2], 6.0);
  monolis_test_check_eq_C1("monolis_matvec_c_test_C22", ry[3],12.0);
  monolis_test_check_eq_C1("monolis_matvec_c_test_C22", ry[4], 7.0);
  monolis_test_check_eq_C1("monolis_matvec_c_test_C22", ry[5],14.0);
  monolis_test_check_eq_C1("monolis_matvec_c_test_C22", ry[6], 8.0);
  monolis_test_check_eq_C1("monolis_matvec_c_test_C22", ry[7],16.0);
  monolis_test_check_eq_C1("monolis_matvec_c_test_C22", ry[8], 3.0);
  monolis_test_check_eq_C1("monolis_matvec_c_test_C22", ry[9], 6.0);

  monolis_finalize(&mat);
}

void monolis_matvec_c_test_C33(){
  MONOLIS mat;
  MONOLIS_COM com;
  int i;
  int n_dof;
  int n_node;
  int n_elem;
  int** elem;
  double _Complex rx[15];
  double _Complex ry[15];

  monolis_std_log_string("monolis_matvec_c_test_C33");

  elem = monolis_alloc_I_2d(elem, 4, 2);

  n_dof = 3;

  n_node = 5;

  n_elem = 4;

  elem[0][0] = 0; elem[0][1] = 1;
  elem[1][0] = 1; elem[1][1] = 2;
  elem[2][0] = 2; elem[2][1] = 3;
  elem[3][0] = 3; elem[3][1] = 4;

  monolis_initialize(&mat);
  monolis_com_initialize_by_self(&com);

  monolis_get_nonzero_pattern_by_simple_mesh_C(&mat, n_node, 2, n_dof, n_elem, elem);

  monolis_add_scalar_to_sparse_matrix_C(&mat, 0, 0, 0, 0, 2.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 0, 1, 0, 0, 1.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 1, 0, 0, 0, 1.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 1, 1, 0, 0, 2.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 1, 2, 0, 0, 3.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 2, 1, 0, 0, 1.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 2, 2, 0, 0, 2.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 2, 3, 0, 0, 4.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 3, 2, 0, 0, 1.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 3, 3, 0, 0, 2.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 3, 4, 0, 0, 5.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 4, 3, 0, 0, 1.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 4, 4, 0, 0, 2.0);

  monolis_add_scalar_to_sparse_matrix_C(&mat, 0, 0, 1, 1, 4.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 0, 1, 1, 1, 2.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 1, 0, 1, 1, 2.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 1, 1, 1, 1, 4.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 1, 2, 1, 1, 6.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 2, 1, 1, 1, 2.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 2, 2, 1, 1, 4.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 2, 3, 1, 1, 8.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 3, 2, 1, 1, 2.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 3, 3, 1, 1, 4.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 3, 4, 1, 1,10.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 4, 3, 1, 1, 2.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 4, 4, 1, 1, 4.0);

  monolis_add_scalar_to_sparse_matrix_C(&mat, 0, 0, 2, 2, 8.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 0, 1, 2, 2, 4.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 1, 0, 2, 2, 4.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 1, 1, 2, 2, 8.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 1, 2, 2, 2,12.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 2, 1, 2, 2, 4.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 2, 2, 2, 2, 8.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 2, 3, 2, 2,16.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 3, 2, 2, 2, 4.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 3, 3, 2, 2, 8.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 3, 4, 2, 2,20.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 4, 3, 2, 2, 4.0);
  monolis_add_scalar_to_sparse_matrix_C(&mat, 4, 4, 2, 2, 8.0);

  for(i = 0; i < 15; ++i){
    rx[i] = 1.0;
  }

  for(i = 0; i < 15; ++i){
    ry[i] = 0.0;
  }

  monolis_matvec_product_C(&mat, &com, rx, ry);

  monolis_test_check_eq_C1("monolis_matvec_c_test_C33", ry[0], 3.0);
  monolis_test_check_eq_C1("monolis_matvec_c_test_C33", ry[1], 6.0);
  monolis_test_check_eq_C1("monolis_matvec_c_test_C33", ry[2],12.0);
  monolis_test_check_eq_C1("monolis_matvec_c_test_C33", ry[3], 6.0);
  monolis_test_check_eq_C1("monolis_matvec_c_test_C33", ry[4],12.0);
  monolis_test_check_eq_C1("monolis_matvec_c_test_C33", ry[5],24.0);
  monolis_test_check_eq_C1("monolis_matvec_c_test_C33", ry[6], 7.0);
  monolis_test_check_eq_C1("monolis_matvec_c_test_C33", ry[7],14.0);
  monolis_test_check_eq_C1("monolis_matvec_c_test_C33", ry[8],28.0);
  monolis_test_check_eq_C1("monolis_matvec_c_test_C33", ry[9], 8.0);
  monolis_test_check_eq_C1("monolis_matvec_c_test_C33", ry[10],16.0);
  monolis_test_check_eq_C1("monolis_matvec_c_test_C33", ry[11],32.0);
  monolis_test_check_eq_C1("monolis_matvec_c_test_C33", ry[12], 3.0);
  monolis_test_check_eq_C1("monolis_matvec_c_test_C33", ry[13], 6.0);
  monolis_test_check_eq_C1("monolis_matvec_c_test_C33", ry[14],12.0);

  monolis_finalize(&mat);
}

void monolis_matvec_c_test(){
  if(monolis_mpi_get_global_comm_size() == 1) return

  monolis_matvec_c_test_R11();
  monolis_matvec_c_test_R22();
  monolis_matvec_c_test_R33();

  monolis_matvec_c_test_C11();
  monolis_matvec_c_test_C22();
  monolis_matvec_c_test_C33();
}
