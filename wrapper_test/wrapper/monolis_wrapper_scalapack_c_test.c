#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "monolis_utils.h"
#include "monolis_wrapper_scalapack_c.h"

void monolis_scalapack_test_1()
{
  int    N_loc = 4;
  int    M = 3;
  double** A;
  double** S;
  double*  V;
  double** D;
  double VD[3][3];
  double SVD[4][3];
  int    comm;
  int    i, j, k;

  monolis_std_log_string("monolis_scalapack_test");

  comm = monolis_mpi_get_global_comm();

  A = monolis_alloc_R_2d(A, 4, 3);
  S = monolis_alloc_R_2d(S, 4, 3);
  V = monolis_alloc_R_1d(V, 3);
  D = monolis_alloc_R_2d(D, 3, 3);

  if(monolis_mpi_get_global_my_rank() == 0){
    A[0][0] = 1.0;
    A[1][0] = 2.0;
    A[2][0] = 3.0;
    A[3][0] = 4.0;
    A[0][1] = 11.0;
    A[1][1] = 12.0;
    A[2][1] = 13.0;
    A[3][1] = 14.0;
    A[0][2] = 21.0;
    A[1][2] = 22.0;
    A[2][2] = 23.0;
    A[3][2] = 24.0;
  } else {
    A[0][0] = 5.0;
    A[1][0] = 6.0;
    A[2][0] = 7.0;
    A[3][0] = 8.0;
    A[0][1] = 15.0;
    A[1][1] = 16.0;
    A[2][1] = 17.0;
    A[3][1] = 18.0;
    A[0][2] = 25.0;
    A[1][2] = 26.0;
    A[2][2] = 27.0;
    A[3][2] = 28.0;
  }

  monolis_scalapack_gesvd_R(N_loc, M, A, S, V, D, comm);

  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      VD[i][j] = V[i]*D[i][j];
    }
  }

  for (i = 0; i < 4; ++i) {
    for (j = 0; j < 3; ++j) {
      SVD[i][j] = 0.0;
      for (k = 0; k < 3; ++k) {
        SVD[i][j] = SVD[i][j] + S[i][k]*VD[k][j];
      }
    }
  }

  for (i = 0; i < 4; ++i) {
    for (j = 0; j < 3; ++j) {
      monolis_test_check_eq_R1("monolis_scalapack_test R", A[i][j], SVD[i][j]);
    }
  }
}

void monolis_scalapack_test_2()
{
  int    N_loc = 2;
  int    M = 6;
  double** A;
  double** S;
  double*  V;
  double** D;
  double VD[4][6];
  double SVD[2][6];
  int    comm;
  int    i, j, k;

  monolis_std_log_string("monolis_scalapack_test");

  comm = monolis_mpi_get_global_comm();

  A = monolis_alloc_R_2d(A, 2, 6);
  S = monolis_alloc_R_2d(S, 2, 4);
  V = monolis_alloc_R_1d(V, 4);
  D = monolis_alloc_R_2d(D, 4, 6);

  if(monolis_mpi_get_global_my_rank() == 0){
    A[0][0] = 1.0;
    A[1][0] = 2.0;
    A[0][1] = 13.0;
    A[1][1] = 14.0;
    A[0][2] = 25.0;
    A[1][2] = 26.0;
    A[0][3] = 37.0;
    A[1][3] = 38.0;
    A[0][4] = 49.0;
    A[1][4] = 50.0;
    A[0][5] = 51.0;
    A[1][5] = 52.0;
  } else {
    A[0][0] = 3.0;
    A[1][0] = 4.0;
    A[0][1] = 15.0;
    A[1][1] = 16.0;
    A[0][2] = 27.0;
    A[1][2] = 28.0;
    A[0][3] = 39.0;
    A[1][3] = 40.0;
    A[0][4] = 41.0;
    A[1][4] = 42.0;
    A[0][5] = 53.0;
    A[1][5] = 54.0;
  }

  monolis_scalapack_gesvd_R(N_loc, M, A, S, V, D, comm);

  if(monolis_mpi_get_global_comm_size() == 1){
    for (i = 0; i < 2; ++i) {
      for (j = 0; j < 6; ++j) {
        VD[i][j] = V[i]*D[i][j];
      }
    }

    for (i = 0; i < 2; ++i) {
      for (j = 0; j < 6; ++j) {
        SVD[i][j] = 0.0;
        for (k = 0; k < 2; ++k) {
          SVD[i][j] = SVD[i][j] + S[i][k]*VD[k][j];
        }
      }
    }
  } else {

    for (i = 0; i < 4; ++i) {
      for (j = 0; j < 6; ++j) {
        VD[i][j] = V[i]*D[i][j];
      }
    }

    for (i = 0; i <2; ++i) {
      for (j = 0; j < 6; ++j) {
        SVD[i][j] = 0.0;
        for (k = 0; k < 4; ++k) {
          SVD[i][j] = SVD[i][j] + S[i][k]*VD[k][j];
        }
      }
    }
  }

  for (i = 0; i < 2; ++i) {
    for (j = 0; j < 6; ++j) {
      monolis_test_check_eq_R1("monolis_scalapack_test R", A[i][j], SVD[i][j]);
    }
  }
}

void monolis_scalapack_test_3()
{
  int    N_loc = 2;
  int    M = 6;
  double** A;
  double** S;
  double*  V;
  double** D;
  double VD[2][6];
  double SVD[2][6];
  int    comm;
  int    i, j, k;

  monolis_std_log_string("monolis_scalapack_test");

  comm = monolis_mpi_get_self_comm();

  A = monolis_alloc_R_2d(A, 2, 6);
  S = monolis_alloc_R_2d(S, 2, 2);
  V = monolis_alloc_R_1d(V, 2);
  D = monolis_alloc_R_2d(D, 2, 6);

  if(monolis_mpi_get_global_my_rank() == 0){
    A[0][0] = 1.0;
    A[1][0] = 2.0;
    A[0][1] = 13.0;
    A[1][1] = 14.0;
    A[0][2] = 25.0;
    A[1][2] = 26.0;
    A[0][3] = 37.0;
    A[1][3] = 38.0;
    A[0][4] = 49.0;
    A[1][4] = 50.0;
    A[0][5] = 51.0;
    A[1][5] = 52.0;
  } else {
    A[0][0] = 3.0;
    A[1][0] = 4.0;
    A[0][1] = 15.0;
    A[1][1] = 16.0;
    A[0][2] = 27.0;
    A[1][2] = 28.0;
    A[0][3] = 39.0;
    A[1][3] = 40.0;
    A[0][4] = 41.0;
    A[1][4] = 42.0;
    A[0][5] = 53.0;
    A[1][5] = 54.0;
  }

  monolis_scalapack_gesvd_R(N_loc, M, A, S, V, D, comm);

  for (i = 0; i < 2; ++i) {
    for (j = 0; j < 6; ++j) {
      VD[i][j] = V[i]*D[i][j];
    }
  }

  for (i = 0; i < 2; ++i) {
    for (j = 0; j < 6; ++j) {
      SVD[i][j] = 0.0;
      for (k = 0; k < 2; ++k) {
        SVD[i][j] = SVD[i][j] + S[i][k]*VD[k][j];
      }
    }
  }

  for (i = 0; i < 2; ++i) {
    for (j = 0; j < 6; ++j) {
      monolis_test_check_eq_R1("monolis_scalapack_test R", A[i][j], SVD[i][j]);
    }
  }
}

void monolis_scalapack_test()
{
  monolis_scalapack_test_1();
  monolis_scalapack_test_2();
  monolis_scalapack_test_3();
}