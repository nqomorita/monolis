#include <stdio.h>
#include <stdlib.h>
#include "monolis.h"

int main(int argc, char *args[]) {
  FILE *fp;
  monolis_com monoCOM;
  char fname[] = "test.mtx";
  int i, j, N, NDOF, NPU, NPL, NZ;
  int method, precond, maxiter, is_scaling, is_reordering, is_init_x;
  int show_iterlog, show_time, show_summary;
  int *index, *item;
  double *A, *X, *B;
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

  monolis_com_initialize(&monoCOM);

  X = (double *)calloc(N*NDOF, sizeof(double));
  B = (double *)calloc(N*NDOF, sizeof(double));

  for (i=0; i<N*NDOF; i++){
    B[i] = (double)i;
  }

  //X[0] = -0.241534067727462;
  //X[1] = 0.300285597715213;
  //X[2] = 0.365565075479389;
  //X[3] = 0.550795593635243;
  //X[4] = 1.077111383108927;
  //X[5] = 0.611179110567108;
  //X[6] = 1.002039983680125;
  //X[7] = 2.846185230518152;

  method = 1;
  precond = 3;
  maxiter = 10;
  tol = 1.0e-8;
  is_scaling    = 1;
  is_reordering = 1;
  is_init_x     = 1;
  show_iterlog  = 1;
  show_time     = 1;
  show_summary  = 1;

  printf("* call monolis\n");
  monolis(&monoCOM, N, N, NZ, NDOF, A, X, B, index, item,
     method, precond, maxiter, tol,
     is_scaling, is_reordering, is_init_x, show_iterlog, show_time, show_summary);

  printf("* monolis result\n");
  for (i=0; i<N; i++){
    printf("%f ", X[i]);
  }
  printf("\n");

  for (i=0; i<N*NDOF; i++){
    B[i] = 0.0;
  }
  monolis_matvec_wrapper(&monoCOM, N, N, NZ, NDOF, A, index, item, X, B);

  printf("* monolis b = Ax\n");
  printf("* monolis result\n");
  for (i=0; i<N; i++){
    printf("%f ", B[i]);
  }
  printf("\n");

  free(X);
  free(B);
  free(A);
  free(index);
  free(item);

  monolis_com_finalize(&monoCOM);

  return 0;
}
