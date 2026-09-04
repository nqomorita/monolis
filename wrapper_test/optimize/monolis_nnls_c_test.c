#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "monolis_solver.h"
#include "monolis_nnls_c_test.h"

void monolis_optimize_nnls_c_test_1()
{
  int    m = 4;
  int    n = 2;
  double** A;
  double*  b;
  double*  x;
  int    max_iter;
  int    i, j, k;
  int    self_comm;
  double tol, tol_inner, residual;

  monolis_std_log_string("monolis_optimize_nnls_test");

  max_iter = 10;
  tol = 1.0e-6;
  tol_inner = 1.0e-12;
  self_comm = monolis_mpi_get_self_comm();

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

  monolis_optimize_nnls_R(A, b, x, m, n, max_iter, tol, &residual, self_comm);

  monolis_test_check_eq_R1("monolis_optimize_nnls_test 1a", x[0], 1.0);
  monolis_test_check_eq_R1("monolis_optimize_nnls_test 1b", x[1], 0.0);

  monolis_optimize_nnls_R_with_sparse_solution(A, b, x, m, n, max_iter, tol, tol_inner, &residual, self_comm);

  monolis_test_check_eq_R1("monolis_optimize_nnls_test 4a", x[0], 1.0);
  monolis_test_check_eq_R1("monolis_optimize_nnls_test 4b", x[1], 0.0);
}

void monolis_optimize_nnls_c_test_2()
{
  int    m = 3;
  int    n = 2;
  double** A;
  double*  b;
  double*  x;
  int    max_iter;
  int    i, j, k;
  int    self_comm;
  double tol, tol_inner, residual;

  max_iter = 10;
  tol = 1.0e-6;
  tol_inner = 1.0e-12;
  self_comm = monolis_mpi_get_self_comm();

  A = monolis_alloc_R_2d(A, m, n);
  x = monolis_alloc_R_1d(x, n);
  b = monolis_alloc_R_1d(b, m);

  A[0][0] = 1.0; A[0][1] = 0.0; 
  A[1][0] = 1.0; A[1][1] = 0.0; 
  A[2][0] = 0.0; A[2][1] = 1.0; 

  b[0] = 2.0;
  b[1] = 1.0;
  b[2] = 1.0;

  monolis_optimize_nnls_R(A, b, x, m, n, max_iter, tol, &residual, self_comm);

  monolis_test_check_eq_R1("monolis_optimize_nnls_test 2a", x[0], 1.5);
  monolis_test_check_eq_R1("monolis_optimize_nnls_test 2b", x[1], 1.0);

  monolis_optimize_nnls_R_with_sparse_solution(A, b, x, m, n, max_iter, tol, tol_inner, &residual, self_comm);

  monolis_test_check_eq_R1("monolis_optimize_nnls_test 5a", x[0], 1.5);
  monolis_test_check_eq_R1("monolis_optimize_nnls_test 5b", x[1], 1.0);

  A[0][0] = 1.0; A[0][1] = 0.0; 
  A[1][0] = 1.0; A[1][1] = 0.0; 
  A[2][0] = 0.0; A[2][1] = 1.0; 

  b[0] =-1.0;
  b[1] =-1.0;
  b[2] =-1.0;

  monolis_optimize_nnls_R(A, b, x, m, n, max_iter, tol, &residual, self_comm);

  monolis_test_check_eq_R1("monolis_optimize_nnls_test 3a", x[0], 0.0);
  monolis_test_check_eq_R1("monolis_optimize_nnls_test 3b", x[1], 0.0);

  monolis_optimize_nnls_R_with_sparse_solution(A, b, x, m, n, max_iter, tol, tol_inner, &residual, self_comm);

  monolis_test_check_eq_R1("monolis_optimize_nnls_test 6a", x[0], 0.0);
  monolis_test_check_eq_R1("monolis_optimize_nnls_test 6b", x[1], 0.0);
}

void monolis_optimize_parallel_nnls_c_test()
{
  int      m = 4;
  int      n = 4;
  int      n_loc, offset;
  double** A_g;
  double** A_loc;
  double*  b;
  double*  x;
  double*  x_ref;
  int      max_iter;
  int      i, j;
  int      comm, self_comm;
  double   tol, tol_inner, residual, residual_ref;

  if(monolis_mpi_get_global_comm_size() > 2) return;

  monolis_std_log_string("monolis_optimize_parallel_nnls_test");

  comm = monolis_mpi_get_global_comm();
  self_comm = monolis_mpi_get_self_comm();
  max_iter = 100;
  tol = 1.0e-6;
  tol_inner = 1.0e-12;

  A_g = monolis_alloc_R_2d(A_g, m, n);
  b = monolis_alloc_R_1d(b, m);
  x_ref = monolis_alloc_R_1d(x_ref, n);

  /* 対角優位な 4x4 行列（全成分が active になるケース） */
  for (i = 0; i < m; ++i) {
    for (j = 0; j < n; ++j) {
      A_g[i][j] = 0.1;
    }
    A_g[i][i] = 2.0;
  }

  b[0] = 1.0;
  b[1] = 2.0;
  b[2] = 3.0;
  b[3] = 4.0;

  /* セルフコミュニケータによる参照解 */
  monolis_optimize_nnls_R_with_sparse_solution(A_g, b, x_ref, m, n, max_iter, tol, tol_inner, &residual_ref, self_comm);

  /* 列ブロック分散（1 ランク: 4 列、2 ランク: 各 2 列） */
  if (monolis_mpi_get_global_comm_size() == 1) {
    n_loc = 4;
    offset = 0;
  } else {
    n_loc = 2;
    offset = 2*monolis_mpi_get_global_my_rank();
  }

  A_loc = monolis_alloc_R_2d(A_loc, m, n_loc);
  x = monolis_alloc_R_1d(x, n_loc);

  for (i = 0; i < m; ++i) {
    for (j = 0; j < n_loc; ++j) {
      A_loc[i][j] = A_g[i][offset + j];
    }
  }

  monolis_optimize_nnls_R_with_sparse_solution(A_loc, b, x, m, n_loc, max_iter, tol, tol_inner, &residual, comm);

  for (i = 0; i < n_loc; ++i) {
    monolis_test_check_eq_R1("monolis_optimize_parallel_nnls_test 1a", x[i], x_ref[offset + i]);
  }
  monolis_test_check_eq_R1("monolis_optimize_parallel_nnls_test 1b", residual, residual_ref);

  /* 疎な解になるケース（b が第 2 列のスカラー倍） */
  for (i = 0; i < m; ++i) {
    b[i] = 2.0*A_g[i][1];
  }

  monolis_optimize_nnls_R_with_sparse_solution(A_g, b, x_ref, m, n, max_iter, tol, tol_inner, &residual_ref, self_comm);

  monolis_optimize_nnls_R_with_sparse_solution(A_loc, b, x, m, n_loc, max_iter, tol, tol_inner, &residual, comm);

  for (i = 0; i < n_loc; ++i) {
    monolis_test_check_eq_R1("monolis_optimize_parallel_nnls_test 2a", x[i], x_ref[offset + i]);
  }
  monolis_test_check_eq_R1("monolis_optimize_parallel_nnls_test 2b", residual, residual_ref);
}

void monolis_optimize_nnls_c_test()
{
  monolis_optimize_nnls_c_test_1();
  monolis_optimize_nnls_c_test_2();
  monolis_optimize_parallel_nnls_c_test();
}