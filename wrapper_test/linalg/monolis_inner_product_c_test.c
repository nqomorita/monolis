#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <complex.h>
#include "monolis.h"
#include "monolis_utils.h"

void monolis_inner_product_c_test(){
  MONOLIS mat;
  int ix[4];
  int iy[4];
  int isum;
  int n_dof;
  double rx[4];
  double ry[4];
  double rsum;
  double complex cx[4];
  double complex cy[4];
  double complex csum;

  monolis_std_log_string("monolis_inner_product_c_test");

  mat.mat.N = 2;

  mat.com.n_internal_vertex = 2;

  n_dof = 2;

  ix[0] = 1; iy[0] = 1;
  ix[1] = 1; iy[1] = 2;
  ix[2] = 1; iy[2] = 3;
  ix[3] = 1; iy[3] = 4;

  monolis_inner_product_I(&mat, n_dof, ix, iy, &isum);

  monolis_test_check_eq_I1("monolis_linalg_test 1", isum, 10);

  rx[0] = 1.0; ry[0] = 1.0;
  rx[1] = 1.0; ry[1] = 2.0;
  rx[2] = 1.0; ry[2] = 3.0;
  rx[3] = 1.0; ry[3] = 4.0;

  monolis_inner_product_R(&mat, n_dof, rx, ry, &rsum);

  monolis_test_check_eq_R1("monolis_linalg_test 2", rsum, 10.0);

  cx[0] = 1.0 + 0.0*I; cy[0] = 1.0 + 1.0*I;
  cx[1] = 1.0 + 0.0*I; cy[1] = 2.0 + 2.0*I;
  cx[2] = 1.0 + 0.0*I; cy[2] = 3.0 + 3.0*I;
  cx[3] = 1.0 + 0.0*I; cy[3] = 4.0 + 4.0*I;

  monolis_inner_product_C(&mat, n_dof, cx, cy, &csum);

  monolis_test_check_eq_R1("monolis_linalg_test 3", csum, 10.0 + 10.0*I);
}
