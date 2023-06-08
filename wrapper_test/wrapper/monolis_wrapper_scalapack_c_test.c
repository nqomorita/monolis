#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "monolis_utils.h"
#include "monolis_wrapper_scalapack_c.h"

void monolis_scalapack_test()
{
  int    N_loc = 4;
  int    M = 3;
  double** A;
  double** S;
  double*  V;
  double** D;
  double VD[3][3];
  double SVD[3][3];
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
      VD[i][j] = V[j]*D[i][j];
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
      printf("%f %f \n", A[i][j], SVD[i][j]);
      monolis_test_check_eq_R1("monolis_scalapack_test R", A[i][j], SVD[i][j]);
    }
  }
}
