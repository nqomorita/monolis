#include <stdio.h>
#include <stdlib.h>
#include "monolis.h"

void monolis_convert_alloc_matrix_c(int N, int NDOF, int NPU, int NPL,
  double *D, double *AU, double *AL, int *indexU, int *indexL, int *itemU, int *itemL, double *X, double *B)
{
  X = (double *)calloc(N*NDOF, sizeof(double));
  B = (double *)calloc(N*NDOF, sizeof(double));
  D  = (double *)calloc(N  *NDOF*NDOF, sizeof(double));
  AU = (double *)calloc(NPU*NDOF*NDOF, sizeof(double));
  AL = (double *)calloc(NPL*NDOF*NDOF, sizeof(double));
  indexU = (int *)calloc(N+1, sizeof(int));
  indexL = (int *)calloc(N+1, sizeof(int));
  itemU  = (int *)calloc(NPU, sizeof(int));
  itemL  = (int *)calloc(NPL, sizeof(int));
}

int main(int argc, char *args[]) {
  FILE *fp;
  char fname[] = "test.mtx";
  int i, j, N, NDOF, NPU, NPL, NZ;
  int method, precond, maxiter, is_scaling;
  int *indexI, *indexJ, *indexU, *indexL, *itemU, *itemL;
  double *A, *D, *AU, *AL, *X, *B;
  double tol;
  char str[16];

  printf("* monolis coo_matrix test in C\n");

  fp = fopen(fname, "r");
    fscanf(fp, "%s", str);
    fscanf(fp, "%d %d %d", &N, &N, &NZ);
    indexI = (int *)calloc(NZ, sizeof(int));
    indexJ = (int *)calloc(NZ, sizeof(int));
    A = (double *)calloc(NZ, sizeof(double));
    for(i=0; i<NZ; i=i+1){
      fscanf(fp, "%d %d %lf", &indexI[i], &indexJ[i], &A[i]);
    }
  fclose(fp);

  NDOF = 1;

  monolis_convert_coo_get_size(&N, &NZ, indexI, indexJ, &NPU, &NPL);

  //monolis_convert_alloc_matrix_c(N, NDOF, NPU, NPL, D, AU, AL, indexU, indexL, itemU, itemL, X, B);
  X = (double *)calloc(N*NDOF, sizeof(double));
  B = (double *)calloc(N*NDOF, sizeof(double));
  D  = (double *)calloc(N  *NDOF*NDOF, sizeof(double));
  AU = (double *)calloc(NPU*NDOF*NDOF, sizeof(double));
  AL = (double *)calloc(NPL*NDOF*NDOF, sizeof(double));
  indexU = (int *)calloc(N+1, sizeof(int));
  indexL = (int *)calloc(N+1, sizeof(int));
  itemU  = (int *)calloc(NPU, sizeof(int));
  itemL  = (int *)calloc(NPL, sizeof(int));

  monolis_convert_coo_get_index(&N, &NZ, indexI, indexJ, &NPU, &NPL, indexU, itemU, indexL, itemL);

  monolis_convert_coo_update_matrix_entry(&N, &NZ, &NDOF, A, indexI, indexJ, &NPU, &NPL, D, AU, AL, indexU, itemU, indexL, itemL);

  for (i=0; i<N*NDOF; i++){
    B[i] = 1.0;
  }

  method = 1;
  precond = 1;
  maxiter = 10;
  tol = 1.0e-8;
  is_scaling = 1;

  monolis_serial(&N, &NDOF, &NPU, &NPL, D, AU, AL, X, B, indexU, itemU, indexL, itemL, &method, &precond, &maxiter, &tol, &is_scaling);

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
