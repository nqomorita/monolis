#include <stdio.h>
#include <stdlib.h>
#include "monolis.h"

int main(int argc, char *args[]) {
  FILE *fp;
  char fname[] = "test.mtx";
  int i, j, N, NDOF, NPU, NPL, NZ;
  int method, precond, maxiter, is_scaling;
  int *index, *item, *indexU, *indexL, *itemU, *itemL;
  double *A, *D, *AU, *AL, *X, *B;
  double tol;
  char str[16];

  printf("* monolis csr_matrix test in C\n");

  fp = fopen(fname, "r");
    fscanf(fp, "%s", str);
    fscanf(fp, "%d %d", &N, &NZ);
    index = (int *)calloc(N+1, sizeof(int));
    item = (int *)calloc(NZ, sizeof(int));
    A = (double *)calloc(NZ, sizeof(double));
    for(i=0; i<N+1; i=i+1){
      fscanf(fp, "%d", &index[i]);
    }
    for(i=0; i<NZ; i=i+1){
      fscanf(fp, "%d", &item[i]);
    }
    for(i=0; i<NZ; i=i+1){
      fscanf(fp, "%lf", &A[i]);
    }
  fclose(fp);

  NDOF = 1;

  monolis_convert_csr_get_size(N, NZ, index, item, &NPU, &NPL);

  X = (double *)calloc(N*NDOF, sizeof(double));
  B = (double *)calloc(N*NDOF, sizeof(double));
  D  = (double *)calloc(N  *NDOF*NDOF, sizeof(double));
  AU = (double *)calloc(NPU*NDOF*NDOF, sizeof(double));
  AL = (double *)calloc(NPL*NDOF*NDOF, sizeof(double));
  indexU = (int *)calloc(N+1, sizeof(int));
  indexL = (int *)calloc(N+1, sizeof(int));
  itemU  = (int *)calloc(NPU, sizeof(int));
  itemL  = (int *)calloc(NPL, sizeof(int));

  for (i=0; i<N*NDOF; i++){
    B[i] = 1.0;
  }

  monolis_convert_csr_get_index(N, NZ, index, item, NPU, NPL, indexU, itemU, indexL, itemL);

  monolis_convert_csr_update_matrix_entry(N, NZ, NDOF, A, index, item, NPU, NPL, D, AU, AL, indexU, itemU, indexL, itemL);

  method = 1;
  precond = 1;
  maxiter = 10;
  tol = 1.0e-8;
  is_scaling = 1;

  monolis_serial(N, NDOF, NPU, NPL, D, AU, AL, X, B, indexU, itemU, indexL, itemL, method, precond, maxiter, tol, is_scaling);

  /** monolis_finalize(monoPRM, monoCOM, monoMAT) **/
  free(X);
  free(B);
  free(D);
  free(AU);
  free(AL);
  free(indexU);
  free(indexL);
  free(itemU);
  free(itemL);

  return 0;
}
