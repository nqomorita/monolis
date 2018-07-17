#include <stdio.h>
#include <stdlib.h>
#include "monolis.h"

int main(int argc, char *args[]) {
  FILE *fp;
  char fname[] = "test.mtx";
  int i, j, N, NDOF, NPU, NPL, NZ;
  int method, precond, maxiter, is_scaling, is_reordering, is_init_x, show_iteration;
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
    B[i] = (double)i;
  }

  X[0] = -0.241534067727462;
  X[1] = 0.300285597715213;
  X[2] = 0.365565075479389;
  X[3] = 0.550795593635243;
  X[4] = 1.077111383108927;
  X[5] = 0.611179110567108;
  X[6] = 1.002039983680125;
  X[7] = 2.846185230518152;

  monolis_convert_csr_get_index(N, NZ, index, item, NPU, NPL, indexU, itemU, indexL, itemL);

  monolis_convert_csr_update_matrix_entry(N, NZ, NDOF, A, index, item, NPU, NPL, D, AU, AL, indexU, itemU, indexL, itemL);

  method = 1;
  precond = 4;
  maxiter = 10;
  tol = 1.0e-8;
  is_scaling = 1;
  is_reordering = 1;
  is_init_x = 1;
  show_iteration = 1;

/*
  printf("* indexU\n");
  for (i=0; i<N+1; i++){
    printf("%d ", indexU[i]);
  }
  printf("\n");

  printf("* indexL\n");
  for (i=0; i<N+1; i++){
    printf("%d ", indexL[i]);
  }
  printf("\n");

  printf("* itemU\n");
  for (i=0; i<NPU; i++){
    printf("%d ", itemU[i]);
  }
  printf("\n");

  printf("* itemL\n");
  for (i=0; i<NPL; i++){
    printf("%d ", itemL[i]);
  }
  printf("\n");
*/

  monolis_serial(N, NDOF, NPU, NPL, D, AU, AL, X, B, indexU, itemU, indexL, itemL, method, precond, maxiter, tol,
  is_scaling, is_reordering, is_init_x, show_iteration);

  printf("* monolis result\n");
  for (i=0; i<N; i++){
    printf("%f ", X[i]);
  }
  printf("\n");

  for (i=0; i<N*NDOF; i++){
    B[i] = 0.0;
  }

  monolis_matvec_serial(N, NDOF, NPU, NPL, D, AU, AL, indexU, itemU, indexL, itemL, X, B);

  printf("* monolis result B\n");
  for (i=0; i<N; i++){
    printf("%f ", B[i]);
  }
  printf("\n");

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