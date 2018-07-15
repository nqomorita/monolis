#include <stdio.h>
#include <stdlib.h>
#include "monolis_c.h"

int main(int argc, char *args[]) {
  FILE *fp;
  char fname[] = "test.mtx";
  int i, j, N, NDOF, NPU, NPL, NZ;
  int method, precond, maxiter, is_scaling;
  int *indexI, *indexJ, *indexU, *indexL, *itemU, *itemL;
  double *Af, *D, *AU, *AL, *X, *B;
  double tol;
  char str[16];

  printf("* monolis coo_matrix test in C\n");

  fp = fopen(fname, "r");
    fscanf(fp, "%s", str);
    fscanf(fp, "%d %d %d", &N, &N, &NZ);
    indexI = (int *)calloc(NZ, sizeof(int));
    indexJ = (int *)calloc(NZ, sizeof(int));
    Af =  (double *)calloc(NZ, sizeof(double));
    for(i=0; i<NZ; i=i+1){
      fscanf(fp, "%d %d %lf", &indexI[i], &indexJ[i], &Af[i]);
    }
  fclose(fp);

  NDOF = 1;

  monolis_convert_coo_get_size(&N, &NZ, indexI, indexJ, &NPU, &NPL);

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

  monolis_convert_coo_get_matrix(&N, &NZ, &NDOF, Af, indexI, indexJ, &NPU, &NPL, D, AU, AL, indexU, itemU, indexL, itemL);

  method = 1;
  precond = 1;
  maxiter = 10;
  tol = 1.0e-8;
  is_scaling = 1;

  monolis_serial(&N, &NDOF, &NPU, &NPL, D, AU, AL, X, B, indexU, itemU, indexL, itemL, &method, &precond, &maxiter, &tol, &is_scaling);

  /** monolis_finalize(monoPRM, monoCOM, monoMAT) **/

  return 0;
}
