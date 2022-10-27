#include <stdio.h>
#include <stdlib.h>
#include "monolis.h"

int main(int argc, char *args[]) {
  MONOLIS mat;
  const char* dir_name;
  int nnode = 5;
  int nelem = 4;
  int **elem;

  elem = (int **)malloc(sizeof(int *) * 5);
  for (int i = 0; i < 5; i++) {
    elem[i] = (int *)malloc(sizeof(int) * 2);
  }

  // initialize
  printf("* test: monolis eigen in C\n");

  printf("* monolis_global_initialize\n");
  monolis_global_initialize();

  dir_name = "./";
  printf("* monolis_initialize\n");
  monolis_initialize(&mat, dir_name);

  // make mat
  elem[0][0] = 0; elem[0][1] = 1;
  elem[1][0] = 1; elem[1][1] = 2;
  elem[2][0] = 2; elem[2][1] = 3;
  elem[3][0] = 3; elem[3][1] = 4;

  monolis_get_nonzero_pattern(&mat, nnode, 2, 1, nelem, elem);

  // coef matri
  //  2  1  0  0  0
  //  1  2  1  0  0
  //  0  1  2  2  0
  //  0  0  1  2  1
  //  0  0  0  1  2
  monolis_add_scalar_to_sparse_matrix(&mat, 2.0, 0, 0, 0, 0);
  monolis_add_scalar_to_sparse_matrix(&mat, 1.0, 0, 1, 0, 0);
  monolis_add_scalar_to_sparse_matrix(&mat, 1.0, 1, 0, 0, 0);
  monolis_add_scalar_to_sparse_matrix(&mat, 2.0, 1, 1, 0, 0);
  monolis_add_scalar_to_sparse_matrix(&mat, 1.0, 1, 2, 0, 0);
  monolis_add_scalar_to_sparse_matrix(&mat, 1.0, 2, 1, 0, 0);
  monolis_add_scalar_to_sparse_matrix(&mat, 2.0, 2, 2, 0, 0);
  monolis_add_scalar_to_sparse_matrix(&mat, 1.0, 2, 3, 0, 0);
  monolis_add_scalar_to_sparse_matrix(&mat, 1.0, 3, 2, 0, 0);
  monolis_add_scalar_to_sparse_matrix(&mat, 2.0, 3, 3, 0, 0);
  monolis_add_scalar_to_sparse_matrix(&mat, 1.0, 3, 4, 0, 0);
  monolis_add_scalar_to_sparse_matrix(&mat, 1.0, 4, 3, 0, 0);
  monolis_add_scalar_to_sparse_matrix(&mat, 2.0, 4, 4, 0, 0);

  // eigen_value
  // 3.73205E+00
  // 3.00000E+00
  // 2.00000E+00
  // 1.00000E+00
  // 2.67949E-01

  // e_mode
  // 2.88675E-01 5.00000E-01 5.77350E-01 5.00000E-01 2.88675E-01
  //-5.00000E-01-5.00000E-01 0.0         5.00000E-01 5.00000E-01
  // 5.77350E-01 0.0        -5.77350E-01 0.0         5.77350E-01
  // 5.00000E-01-5.00000E-01 0.0         5.00000E-01-5.00000E-01
  // 2.88675E-01-5.00000E-01 5.77350E-01-5.00000E-01 2.88675E-01

  //
  int n_get_eigen = 5;
  int maxiter = 10;
  double ths = 1.0e-8;
  double *eig_val;
  double **eig_mode;
  bool *is_Dirichlet_bc;

  eig_val = (double *)malloc(sizeof(double) * 5);

  eig_mode = (double **)malloc(sizeof(double *) * 5);
  for (int i = 0; i < 5; i++) {
    eig_mode[i] = (double *)malloc(sizeof(double) * 5);
  }

  is_Dirichlet_bc = (bool *)malloc(sizeof(bool) * 5);
  for (int i = 0; i < 5; i++) {
    is_Dirichlet_bc[i] = false;
  }

  printf("* monolis_eigen_standard_lanczos\n");
  monolis_eigen_inverted_standard_lanczos(&mat, n_get_eigen, ths, maxiter, eig_val, eig_mode, is_Dirichlet_bc);

  printf("* eigen value\n");
  for (int i = 0; i < 5; i++) {
    printf("%f \n", eig_val[i]);
  }

  printf("* eigen mode\n");
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
    printf("%f ", eig_mode[i][j]);
    }
    printf("\n");
  }

  printf("* monolis_eigen_standard_lanczos\n");
  monolis_eigen_standard_lanczos(&mat, n_get_eigen, ths, maxiter, eig_val, eig_mode, is_Dirichlet_bc);

  printf("* eigen value\n");
  for (int i = 0; i < 5; i++) {
    printf("%f \n", eig_val[i]);
  }

  printf("* eigen mode\n");
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
    printf("%f ", eig_mode[i][j]);
    }
    printf("\n");
  }

  // finalize
  printf("* monolis_finalize\n");
  monolis_finalize(&mat);

  printf("* monolis_global_finalize\n");
  monolis_global_finalize();
}
