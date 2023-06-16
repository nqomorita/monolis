#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "monolis_utils.h"

void monolis_vec_copy_I(
  int  n,
  int  n_dof,
  int* x,
  int* y)
{
  int i;
  for (i = 0; i < n*n_dof; ++i) {
    y[i] = x[i];
  }
}

void monolis_vec_copy_R(
  int     n,
  int     n_dof,
  double* x,
  double* y)
{
  int i;
  for (i = 0; i < n*n_dof; ++i) {
    y[i] = x[i];
  }
}

void monolis_vec_copy_C(
  int              n,
  int              n_dof,
  double _Complex* x,
  double _Complex* y)
{
  int i;
  for (i = 0; i < n*n_dof; ++i) {
    y[i] = x[i];
  }
}
