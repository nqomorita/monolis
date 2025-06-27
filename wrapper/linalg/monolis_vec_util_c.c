#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "monolis_utils.h"

void monolis_vec_copy_I(
  int  m,
  int* x,
  int* y)
{
  int i;
  for (i = 0; i < m; ++i) {
    y[i] = x[i];
  }
}

void monolis_vec_copy_R(
  int     m,
  double* x,
  double* y)
{
  int i;
  for (i = 0; i < m; ++i) {
    y[i] = x[i];
  }
}

void monolis_vec_copy_C(
  int              m,
  double _Complex* x,
  double _Complex* y)
{
  int i;
  for (i = 0; i < m; ++i) {
    y[i] = x[i];
  }
}

