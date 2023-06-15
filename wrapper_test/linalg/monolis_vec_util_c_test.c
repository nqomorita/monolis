#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <complex.h>
#include "monolis.h"
#include "monolis_utils.h"

void monolis_vec_util_c_test(){
  int ix[4];
  int iy[4];
  double rx[4];
  double ry[4];
  double _Complex cx[4];
  double _Complex cy[4];

  monolis_std_log_string("monolis_vec_util_c_test");

  ix[0] = 1;
  ix[1] = 2;
  ix[2] = 3;
  ix[3] = 4;

  monolis_vec_copy_I(1, 4, ix, iy);

  monolis_test_check_eq_I1("monolis_vec_copy_I", ix[0], iy[0]);
  monolis_test_check_eq_I1("monolis_vec_copy_I", ix[1], iy[1]);
  monolis_test_check_eq_I1("monolis_vec_copy_I", ix[2], iy[2]);
  monolis_test_check_eq_I1("monolis_vec_copy_I", ix[3], iy[3]);

  rx[0] = 1.0;
  rx[1] = 2.0;
  rx[2] = 3.0;
  rx[3] = 4.0;

  monolis_vec_copy_R(1, 4, rx, ry);

  monolis_test_check_eq_R1("monolis_vec_copy_R", rx[0], ry[0]);
  monolis_test_check_eq_R1("monolis_vec_copy_R", rx[1], ry[1]);
  monolis_test_check_eq_R1("monolis_vec_copy_R", rx[2], ry[2]);
  monolis_test_check_eq_R1("monolis_vec_copy_R", rx[3], ry[3]);

  cx[0] = 1.0;
  cx[1] = 2.0;
  cx[2] = 3.0;
  cx[3] = 4.0;

  monolis_vec_copy_C(1, 4, cx, cy);

  monolis_test_check_eq_R1("monolis_vec_copy_C", cx[0], cy[0]);
  monolis_test_check_eq_R1("monolis_vec_copy_C", cx[1], cy[1]);
  monolis_test_check_eq_R1("monolis_vec_copy_C", cx[2], cy[2]);
  monolis_test_check_eq_R1("monolis_vec_copy_C", cx[3], cy[3]);
}
