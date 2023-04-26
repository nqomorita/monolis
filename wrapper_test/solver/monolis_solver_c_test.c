#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <complex.h>
#include "monolis.h"

void monolis_solve_c_test_R(){
  MONOLIS mat;
  MONOLIS_COM com;
  int n_node;
  int n_base;
  int n_dof;
  int n_elem;
  int i, j, k;
  int** elem;
  double val;
  double a[20];
  double b[20];

  monolis_std_log_string("monolis_solve_c_test_R");

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

  monolis_initialize(&mat);
  monolis_com_initialize_by_self(&com);

  monolis_get_nonzero_pattern_by_simple_mesh_R(&mat, n_node, n_base, n_dof, n_elem, elem);

  for(i = 0; i < n_node; ++i){
    for(j = 0; j < n_dof; ++j){
      val = rand()%1001 + 5000.0;
      monolis_set_scalar_to_sparse_matrix_R(&mat, i, i, j, j, val);
    }
  }

  for(i = 0; i < n_elem; ++i){
    for(j = 0; j < n_dof; ++j){
      for(k = 0; k < n_dof; ++k){
        val = rand()%1001;
        monolis_set_scalar_to_sparse_matrix_R(&mat, elem[i][0], elem[i][1], j, k, val);
        monolis_set_scalar_to_sparse_matrix_R(&mat, elem[i][1], elem[i][0], j, k, val);
      }
    }
  }

  for(i = 0; i < 20; ++i){
    a[i] = 1.0;
  }

  monolis_matvec_product_R(&mat, &com, a, b);

  for(i = 0; i < 20; ++i){
    a[i] = 0.0;
  }

  monolis_set_tolerance(&mat, 1.0e-10);

  monolis_solve_R(&mat, &com, b, a);

  for(i = 0; i < 20; ++i){
    monolis_test_check_eq_R1("monolis_solve_c_test R", a[i], 1.0);
  }

  monolis_finalize(&mat);
}

void monolis_solve_c_test_C(){
  MONOLIS mat;
  MONOLIS_COM com;
  int n_node;
  int n_base;
  int n_dof;
  int n_elem;
  int i, j, k;
  int** elem;
  double _Complex val;
  double _Complex a[20];
  double _Complex b[20];

  monolis_std_log_string("monolis_solve_c_test_C");

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

  monolis_initialize(&mat);
  monolis_com_initialize_by_self(&com);

  monolis_get_nonzero_pattern_by_simple_mesh_C(&mat, n_node, n_base, n_dof, n_elem, elem);

  for(i = 0; i < n_node; ++i){
    for(j = 0; j < n_dof; ++j){
      val = rand()%1001 + 5000.0;
      monolis_set_scalar_to_sparse_matrix_C(&mat, i, i, j, j, val + 0.1*val*I);
    }
  }

  for(i = 0; i < n_elem; ++i){
    for(j = 0; j < n_dof; ++j){
      for(k = 0; k < n_dof; ++k){
        val = rand()%1001;
        monolis_set_scalar_to_sparse_matrix_C(&mat, elem[i][0], elem[i][1], j, k, val + 0.1*val*I);
        monolis_set_scalar_to_sparse_matrix_C(&mat, elem[i][1], elem[i][0], j, k, val + 0.1*val*I);
      }
    }
  }

  for(i = 0; i < 20; ++i){
    a[i] = 1.0 + 1.0*I;
  }

  monolis_matvec_product_C(&mat, &com, a, b);

  for(i = 0; i < 20; ++i){
    a[i] = 0.0 + 0.0*I;
  }

  monolis_set_method(&mat, MONOLIS_ITER_COCG);

  monolis_set_tolerance(&mat, 1.0e-10);

  monolis_solve_C(&mat, &com, b, a);

  for(i = 0; i < 20; ++i){
    monolis_test_check_eq_C1("monolis_solve_c_test C", a[i], 1.0 + 1.0*I);
  }

  monolis_finalize(&mat);
}

void monolis_solve_c_test(){
  monolis_solve_c_test_R();
  monolis_solve_c_test_C();
}