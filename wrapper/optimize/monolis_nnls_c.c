#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "monolis_utils.h"
#include "monolis_nnls_c.h"

void monolis_optimize_nnls_R(
  double** A,
  double*  b,
  double*  x,
  int      m,
  int      n,
  int      max_iter,
  double   tol,
  double*  residual)
{
  int     i, j;
  double* A_tmp;

  A_tmp = monolis_alloc_R_1d(A_tmp, m*n);

  for (i = 0; i < n; ++i) {
    for (j = 0; j < m; ++j) {
      A_tmp[i*m + j] = A[j][i];
    }
  }

  monolis_optimize_nnls_R_c_main(
    A_tmp,
    b,
    x,
    m,
    n,
    max_iter,
    tol,
    residual);
}

void monolis_optimize_nnls_R_with_sparse_solution(
  double** A,
  double*  b,
  double*  x,
  int      m,
  int      n,
  int      max_iter,
  double   tol,
  double*  residual)
{
  int     i, j;
  double* A_tmp;

  A_tmp = monolis_alloc_R_1d(A_tmp, m*n);

  for (i = 0; i < n; ++i) {
    for (j = 0; j < m; ++j) {
      A_tmp[i*m + j] = A[j][i];
    }
  }

  monolis_optimize_nnls_R_with_sparse_solution_c_main(
    A_tmp,
    b,
    x,
    m,
    n,
    max_iter,
    tol,
    residual);
}
