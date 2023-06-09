#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "monolis_utils.h"
#include "monolis_wrapper_scalapack_c.h"

void monolis_scalapack_gesvd_R(
  int      N_loc,
  int      M,
  double** A,
  double** S,
  double*  V,
  double** D,
  int      comm)
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
    comm);

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
