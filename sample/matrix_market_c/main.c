#include <stdio.h>
#include <stdlib.h>
#include "monolis_c.h"

int main(int argc, char *args[]) {
  FILE *fp;
  char fname[] = "test.mtx";
  int i, j, N, NDOF, NPU, NPL, Nf, NZf, NDOFf;
  int method, precond, maxiter, is_scaling;
  int *indexI, *indexJ, *indexU, *indexL, *itemU, *itemL;
  double *Af, *D, *AU, *AL, *X, *B;
  double tol;
  char str[16];

  printf("* monolis coo_matrix test in C\n");

  fp = fopen(fname, "r");
    fscanf(fp, "%s", str);
    fscanf(fp, "%d %d %d", &Nf, &Nf, &NZf);
    indexI = (int *)calloc(NZf, sizeof(int));
    indexJ = (int *)calloc(NZf, sizeof(int));
    Af =  (double *)calloc(NZf, sizeof(double));
    for(i=0; i<NZf; i=i+1){
      fscanf(fp, "%d %d %lf", &indexI[i], &indexJ[i], &Af[i]);
    }
  fclose(fp);

  NDOFf = 1;

  monolis_convert_test(Nf);

  //monolis_convert_coo_matrix(Nf, NZf, NDOFf, Af, indexI, indexJ, N, NDOF, NPU, NPL, D, AU, AL, indexU, itemU, indexL, itemL);

  X = (double *)calloc(Nf, sizeof(double));
  B = (double *)calloc(Nf, sizeof(double));

  method = 1;
  precond = 1;
  maxiter = 1000;
  tol = 1.0e-8;
  is_scaling = 1;

  //monolis_serial(N, NDOF, NPU, NPL, D, AU, AL, X, B, indexU, itemU, indexL, itemL, method, precond, maxiter, tol, is_scaling);

  /** monolis_finalize(monoPRM, monoCOM, monoMAT) **/

  return 0;
}
