#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <stdbool.h>
#include "monolis.h"

void allocate_2d_matrix_D(
  double*** A,
  int dim1,
  int dim2)
{
  *A = (double**)calloc(dim1, sizeof(double*));
  for(int i = 0; i < dim1; i++){
    (*A)[i] = (double*)calloc(dim2, sizeof(double));
  }
  return ;
}

void allocate_2d_matrix_I(
  int*** A,
  int dim1,
  int dim2)
{
  *A = (int**)calloc(dim1, sizeof(int*));
  for(int i = 0; i < dim1; i++){
    (*A)[i] = (int*)calloc(dim2, sizeof(int));
  }
  return ;
}

void monolis_input_mesh_node_c(
  char *fname,
  int* nnode,
  double*** node)
{
  FILE *fp = NULL;
  double a,b,c;
  char buf[100];
  fp = fopen(fname, "r");
  if(fp == NULL){
    printf("%s cant open file\n", fname);
    return;
  }
  fscanf(fp, "%d", nnode);
  printf("nnode = %d\n",*nnode);
  allocate_2d_matrix_D(node, *nnode, 3);
  for(int i = 0; i < *nnode; i++){
    fscanf(fp, "%lf %lf %lf", &(*node)[i][0], &(*node)[i][1], &(*node)[i][2]);
  }
  fclose(fp);
  printf("succeed input %s\n", fname);
  return;
}

void monolis_input_mesh_elem_c(
  char *fname,
  int* nelem,
  int* nbase,
  int*** elem)
{
  FILE *fp = NULL;
  // printf("%s\n", fname);
  fp = fopen(fname, "r");
  if(fp == NULL){
    printf("%scant open file\n", fname);
    return;
  }
  fscanf(fp, "%d %d", nelem, nbase);
  printf("nbase, nelem = %d %s %d\n", *nbase, ",", *nelem);
  allocate_2d_matrix_I(elem, *nelem, *nbase);
  for(int i = 0; i < *nelem; i++){
    fscanf(fp, "%d %d", &(*elem)[i][0], &(*elem)[i][1]);  //要素毎節点数を可変にする必要あり
  }
  fclose(fp);
  printf("succeed input %s\n", fname);
  return;
}

void monolis_input_id_c(
  char *fname,
  int** global_eid)
{
  int nid, i1, i2;
  FILE *fp;
  printf("%s\n", fname);
  fp = fopen(fname, "r");
  if(fp == NULL){
    printf("%scant open file\n", fname);
    return;
  }

  fscanf(fp, "%d", &nid);
  printf("nid %d\n", nid);

  *global_eid = (int*)calloc(nid, sizeof(int));

  for(int i = 0; i < nid; i++){
    fscanf(fp, "%d %d %d", &i1, &i2, &(*global_eid)[i]);  //要素毎節点数を可変にする必要あり
  }
  fclose(fp);

  printf("succeed input %s\n", fname);
  return;
}

void monolis_input_coef_c(
  char *fname,
  int* ncoef,
  double*** coef)
{
  FILE *fp = NULL;
  double a,b,c;
  char buf[100];
  // printf("%s\n", fname);
  fp = fopen(fname, "r");
  if(fp == NULL){
    printf("%s cant open file\n", fname);
    return;
  }
  fscanf(fp, "%d", ncoef);
  allocate_2d_matrix_D(coef, *ncoef, 2);
  for(int i = 0; i < *ncoef; i++){
    fscanf(fp, "%lf", &(*coef)[i][0]);
  }
  fclose(fp);
  printf("succeed input %s\n", fname);
  return;
}

void get_coef_val_c(
  int ncoef,
  double** coef,
  int eid,
  double* val)
{
  *val = coef[eid][0];
  return;
}

void monolis_solver_COCG_test(){
  int nnode, nelem, nbase_func, i, j, id, ncoef, eid[2], my_rank, ierr;
  double real, imag;
  double val;
  int* global_eid;
  int** elem;
  double** coef;
  double** node;
  double* a;
  double* b;
  double* b_th;
  char *fname;
  MONOLIS mat;
  FILE *fp;

  monolis_initialize(&mat, "./"); // 疎行列変数の初期化

  fname = monolis_get_input_filename("node.dat");
  monolis_input_mesh_node_c(fname, &nnode, &node);

  fname = monolis_get_input_filename("elem.dat");
  monolis_input_mesh_elem_c(fname, &nelem, &nbase_func, &elem);

  if(monolis_get_global_commsize() > 1){
    fname = monolis_get_input_filename("connectivity.id");
    monolis_input_id_c(fname, &global_eid);
  }
  else{
    global_eid = (int*)calloc(nelem, sizeof(int));
    for(i = 0; i < nelem; i++){
      global_eid[i] = i;
    }
  }

  printf("nnode nelem %d %d \n", nnode, nelem);
  monolis_get_nonzero_pattern(&mat, nnode, 2, 1, nelem, elem);

  fname = "coef.dat";
  monolis_input_coef_c(fname, &ncoef, &coef);

  printf("monolis_add_scalar_to_sparse_matrix\n");
  for(i = 0; i < nelem; i++){
    for(j = 0; j < 2; j++){
      eid[j] = elem[i][j];
    }
    get_coef_val_c(ncoef, coef, global_eid[i], &val);

    if(eid[0] == eid[1]){
      monolis_add_scalar_to_sparse_matrix(&mat, val, eid[0], eid[1], 0, 0);
    }
    else{
      monolis_add_scalar_to_sparse_matrix(&mat, val, eid[0], eid[1], 0, 0);
      monolis_add_scalar_to_sparse_matrix(&mat, val, eid[1], eid[0], 0, 0);
    }
  }

  printf("monolis_matvec_product\n");
  a = (double *)calloc(nnode, sizeof(double));
  b = (double *)calloc(nnode, sizeof(double));
  b_th = (double *)calloc(nnode, sizeof(double));
  if(a == NULL) exit(0);
  if(b == NULL) exit(0);
  if(b_th == NULL) exit(0);
  real = 1.0;
  imag = 1.0;
  for(i = 0; i < nnode; i++){
    b_th[i] = 1.0;
  }

  monolis_matvec_product(&mat, b_th, b);

  monolis_set_method(&mat, monolis_iter_CG);
  monolis_set_precond(&mat, monolis_prec_DIAG);
  monolis_set_maxiter(&mat, 200000);
  monolis_set_tolerance(&mat, 1.0e-8);
  monolis_show_iterlog(&mat, false);
  monolis_show_summary(&mat, true);

  printf("get comm free\n");

  free(mat.com.recv_neib_pe);
  free(mat.com.recv_index);
  free(mat.com.recv_item);
  free(mat.com.send_neib_pe);
  free(mat.com.send_index);
  free(mat.com.send_item);

  printf("get comm\n");

  monolis_com_get_comm_table(&mat, mat.com.internal_nnode, mat.mat.NP, mat.com.global_node_id);

  printf("monolis_solve\n");

  monolis_solve(&mat, b, a);

  monolis_finalize(&mat);
  return;
}


int main()
{
  monolis_global_initialize();

  monolis_solver_COCG_test();

  monolis_global_finalize();

  return 0;
}
