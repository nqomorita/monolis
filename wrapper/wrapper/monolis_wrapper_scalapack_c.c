#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "monolis_utils.h"
#include "monolis_wrapper_scalapack_c.h"

void monolis_scalapack_comm_initialize(
  int      comm,
  int*     scalapack_comm)
{
  monolis_scalapack_comm_initialize_c_main(comm, scalapack_comm);
}

void monolis_scalapack_comm_finalize(
  int      scalapack_comm){
  monolis_scalapack_comm_finalize_c_main(scalapack_comm);
}

void monolis_scalapack_gesvd_R(
  int      N_loc,
  int      M,
  double** A,
  double** S,
  double*  V,
  double** D,
  int      comm,
  int      scalapack_comm)
{
  int     P;
  int     i, j;
  int     N;
  double* A_tmp;
  double* S_tmp;
  double* D_tmp;

  N = N_loc;
  monolis_allreduce_I(1, &N, MONOLIS_MPI_SUM, comm);

  if(N < M){
    P = N;
  } else {
    P = M;
  }

  A_tmp = monolis_alloc_R_1d(A_tmp, N_loc*M);
  S_tmp = monolis_alloc_R_1d(S_tmp, N_loc*P);
  D_tmp = monolis_alloc_R_1d(D_tmp, P*M);

  for (i = 0; i < M; ++i) {
    for (j = 0; j < N_loc; ++j) {
      A_tmp[i*N_loc + j] = A[j][i];
    }
  }

  monolis_scalapack_gesvd_R_c_main(
    N_loc,
    M,
    P,
    A_tmp,
    S_tmp,
    V,
    D_tmp,
    comm,
    scalapack_comm);

  for (i = 0; i < P; ++i) {
    for (j = 0; j < N_loc; ++j) {
      S[j][i] = S_tmp[i*N_loc + j];
    }
  }

  for (i = 0; i < M; ++i) {
    for (j = 0; j < P; ++j) {
      D[j][i] = D_tmp[i*P + j];
    }
  }
}

/*
void monolis_scalapack_getrf_R(
  int      N_loc,
  int      N,
  double** A,
  int*     ipiv,
  int      comm,
  int      scalapack_comm)
{
  int     i, j;
  double* A_tmp;

  A_tmp = monolis_alloc_R_1d(A_tmp, N_loc*N);

  for (i = 0; i < N; ++i) {
    for (j = 0; j < N_loc; ++j) {
      A_tmp[i*N_loc + j] = A[j][i];
    }
  }

  monolis_scalapack_getrf_R_c_main(
    N_loc,
    N,
    A_tmp,
    ipiv,
    comm,
    scalapack_comm);

  for (i = 0; i < N; ++i) {
    for (j = 0; j < N_loc; ++j) {
      A[j][i] = A_tmp[i*N_loc + j];
    }
  }
}
*/

void monolis_scalapack_getrs_R(
  int      N_loc,
  int      N,
  int      NRHS,
  double** A,
  int*     ipiv,
  double** B,
  int      comm,
  int      scalapack_comm)
{
  int     i, j;
  double* A_tmp;
  double* B_tmp;

  A_tmp = monolis_alloc_R_1d(A_tmp, N_loc*N);
  B_tmp = monolis_alloc_R_1d(B_tmp, N_loc*NRHS);

  for (i = 0; i < N; ++i) {
    for (j = 0; j < N_loc; ++j) {
      A_tmp[i*N_loc + j] = A[j][i];
    }
  }

  for (i = 0; i < NRHS; ++i) {
    for (j = 0; j < N_loc; ++j) {
      B_tmp[i*N_loc + j] = B[j][i];
    }
  }

  monolis_scalapack_getrs_R_c_main(
    N_loc,
    N,
    NRHS,
    A_tmp,
    ipiv,
    B_tmp,
    comm,
    scalapack_comm);

  for (i = 0; i < NRHS; ++i) {
    for (j = 0; j < N_loc; ++j) {
      B[j][i] = B_tmp[i*N_loc + j];
    }
  }
}
