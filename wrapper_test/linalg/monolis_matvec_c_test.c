#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "monolis.h"
#include "monolis_utils.h"

void monolis_matvec_c_test(){

  monolis_std_log_string("monolis_matvec_c_test");

  monolis_matvec_product_R(mat, rx, ry);

  monolis_matvec_product_C(mat, cx, cy);

}