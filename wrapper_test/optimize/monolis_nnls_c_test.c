#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "monolis_solver.h"
#include "monolis_nnls_c_test.h"

void monolis_optimize_nnls_c_test_1()
{
  MONOLIS_COM com;
  int    m = 4;
  int    n = 2;
  double** A;
  double*  b;
  double*  x;
  int    max_iter;
  int    i, j, k;
  double tol, residual;

  monolis_std_log_string("monolis_optimize_nnls_test");

  monolis_com_initialize_by_self(&com);
  max_iter = 10;
  tol = 1.0e-6;

  A = monolis_alloc_R_2d(A, m, n);
  x = monolis_alloc_R_1d(x, n);
  b = monolis_alloc_R_1d(b, m);

  A[0][0] = 1.0; A[0][1] = 1.0; 
  A[1][0] = 1.0; A[1][1] = 1.0; 
  A[2][0] = 1.0; A[2][1] = 1.0; 
  A[3][0] = 1.0; A[3][1] = 1.0; 

  b[0] = 1.0;
  b[1] = 1.0;
  b[2] = 1.0;
  b[3] = 1.0;

  monolis_optimize_nnls_R(A, b, x, m, n, max_iter, tol, &residual, &com);

  monolis_test_check_eq_R1("monolis_optimize_nnls_test 1a", x[0], 1.0);
  monolis_test_check_eq_R1("monolis_optimize_nnls_test 1b", x[1], 0.0);

  A[0][0] = 1.0; A[0][1] = 1.0; 
  A[1][0] = 1.0; A[1][1] = 1.0; 
  A[2][0] = 1.0; A[2][1] = 1.0; 
  A[3][0] = 1.0; A[3][1] = 1.0; 

  b[0] = 1.0;
  b[1] = 1.0;
  b[2] = 1.0;
  b[3] = 1.0;

  monolis_optimize_nnls_R_with_sparse_solution(A, b, x, m, n, max_iter, tol, &residual, &com);

  monolis_test_check_eq_R1("monolis_optimize_nnls_test 4a", x[0], 1.0);
  monolis_test_check_eq_R1("monolis_optimize_nnls_test 4b", x[1], 0.0);
}

void monolis_optimize_nnls_c_test_2()
{
  MONOLIS_COM com;
  int    m = 3;
  int    n = 2;
  double** A;
  double*  b;
  double*  x;
  int    max_iter;
  int    i, j, k;
  double tol, residual;

  monolis_com_initialize_by_self(&com);
  max_iter = 10;
  tol = 1.0e-6;

  A = monolis_alloc_R_2d(A, m, n);
  x = monolis_alloc_R_1d(x, n);
  b = monolis_alloc_R_1d(b, m);

  A[0][0] = 1.0; A[0][1] = 0.0; 
  A[1][0] = 1.0; A[1][1] = 0.0; 
  A[2][0] = 0.0; A[2][1] = 1.0; 

  b[0] = 2.0;
  b[1] = 1.0;
  b[2] = 1.0;

  monolis_optimize_nnls_R(A, b, x, m, n, max_iter, tol, &residual, &com);

  monolis_test_check_eq_R1("monolis_optimize_nnls_test 2a", x[0], 1.5);
  monolis_test_check_eq_R1("monolis_optimize_nnls_test 2b", x[1], 1.0);

  A[0][0] = 1.0; A[0][1] = 0.0; 
  A[1][0] = 1.0; A[1][1] = 0.0; 
  A[2][0] = 0.0; A[2][1] = 1.0; 

  b[0] =-1.0;
  b[1] =-1.0;
  b[2] =-1.0;

  monolis_optimize_nnls_R(A, b, x, m, n, max_iter, tol, &residual, &com);

  monolis_test_check_eq_R1("monolis_optimize_nnls_test 3a", x[0], 0.0);
  monolis_test_check_eq_R1("monolis_optimize_nnls_test 3b", x[1], 0.0);
}

void monolis_optimize_nnls_c_test()
{
  monolis_optimize_nnls_c_test_1();
  monolis_optimize_nnls_c_test_2();
}