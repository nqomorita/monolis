#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <stdbool.h>
#include "monolis.h"

void monolis_input_mesh_node_c(
  const char *fname,
  int* nnode,
  double*** node)
{
  int i;
  FILE *fp = NULL;
  double a,b,c;
  char buf[100];
  fp = fopen(fname, "r");
  if(fp == NULL){
    printf("%s cant open file\n", fname);
    return;
  }
  fscanf(fp, "%d", nnode);

  *node = monolis_alloc_R_2d(*node, *nnode, 3);

  for(i = 0; i < *nnode; i++){
    fscanf(fp, "%lf %lf %lf", &(*node)[i][0], &(*node)[i][1], &(*node)[i][2]);
  }
  fclose(fp);
  return;
}

void monolis_input_mesh_elem_c(
  const char *fname,
  int* nelem,
  int* nbase,
  int*** elem)
{
  int i;
  FILE *fp = NULL;
  fp = fopen(fname, "r");
  if(fp == NULL){
    printf("%scant open file\n", fname);
    return;
  }
  fscanf(fp, "%d %d", nelem, nbase);

  *elem = monolis_alloc_I_2d(*elem, *nelem, *nbase);

  for(i = 0; i < *nelem; i++){
    fscanf(fp, "%d %d", &(*elem)[i][0], &(*elem)[i][1]);
    //(*elem)[i][0] -= 1;
    //(*elem)[i][1] -= 1;
  }
  fclose(fp);
  return;
}

void monolis_input_id_c(
  const char *fname,
  int** global_eid)
{
  int i;
  int nid, kn;
  FILE *fp;
  char ctmp[1024];
  fp = fopen(fname, "r");
  if(fp == NULL){
    printf("%s cant open file\n", fname);
    return;
  }

  fscanf(fp, "%s", ctmp);
  fscanf(fp, "%d %d", &nid, &kn);

  *global_eid = monolis_alloc_I_1d(*global_eid, nid);

  for(i = 0; i < nid; i++){
    fscanf(fp, "%d", &(*global_eid)[i]);
    //(*global_eid)[i] -= 1;
  }
  fclose(fp);
  return;
}

void input_coef(
  const char *fname,
  int* n_coef,
  double** coef)
{
  int i;
  FILE *fp = NULL;
  double a,b,c;
  char buf[MONOLIS_CHARLEN];

  fp = fopen(fname, "r");
  if(fp == NULL){
    printf("%s cant open file\n", fname);
    return;
  }

  fscanf(fp, "%d", n_coef);

  *coef = monolis_alloc_R_1d(*coef, *n_coef);

  for(i = 0; i < *n_coef; i++){
    fscanf(fp, "%lf", &(*coef)[i]);
  }
  fclose(fp);
}

void monolis_solver_parallel_R_test(){
  MONOLIS mat;
  const char* fname;
  int n_node, n_elem, n_base, n_id, n_coef;
  int eid[2], iter, prec, i, j;
  int* global_eid;
  int** elem;
  double val;
  double* coef;
  double* a;
  double* b;
  double* c;
  double** node;

  fname = monolis_get_global_input_file_name("parted.0", "node.dat");
  monolis_input_mesh_node_c(fname, &n_node, &node);

  fname = monolis_get_global_input_file_name("parted.0", "elem.dat");
  monolis_input_mesh_elem_c(fname, &n_elem, &n_base, &elem);

  if(monolis_mpi_get_global_comm_size() > 1){
    fname = monolis_get_global_input_file_name("parted.0", "elem.dat.id");
    monolis_input_id_c(fname, &global_eid);
  } else {
    global_eid = monolis_alloc_I_1d(global_eid, n_elem);
    for(i = 0; i < n_elem; i++){
      global_eid[i] = i;
    }
  }
  monolis_initialize(&mat);

  monolis_get_nonzero_pattern_by_simple_mesh_R(&mat, n_node, 2, 1, n_elem, elem);

  fname = "coef.dat";
  input_coef(fname, &n_coef, &coef);

  for(i = 0; i < n_elem; i++){
    for(j = 0; j < 2; j++){
      eid[j] = elem[i][j];
    }

    val = coef[global_eid[i]];

    if(eid[0] == eid[1]){
      monolis_add_scalar_to_sparse_matrix_R(&mat, eid[0], eid[1], 0, 0, val);
    } else{
      monolis_add_scalar_to_sparse_matrix_R(&mat, eid[0], eid[1], 0, 0, val);
      monolis_add_scalar_to_sparse_matrix_R(&mat, eid[1], eid[0], 0, 0, val);
    }
  }

  a = monolis_alloc_R_1d(a, n_node);
  b = monolis_alloc_R_1d(b, n_node);
  c = monolis_alloc_R_1d(c, n_node);

  for(i = 0; i < n_node; i++){
    a[i] = 1.0;
  }

  monolis_matvec_product_R(&mat, a, c);

  monolis_set_maxiter(&mat, 1000);
  monolis_set_tolerance(&mat, 1.0e-8);
  monolis_show_timelog(&mat, true);
  monolis_show_iterlog(&mat, true);
  monolis_show_summary(&mat, true);

  for (iter = 1; iter < 9; ++iter) {
    for (prec = 0; prec < 3; ++prec) {

      for(i = 0; i < n_node; i++){
        a[i] = 0.0;
        b[i] = c[i];
      }

      monolis_set_method(&mat, iter);
      monolis_set_precond(&mat, prec);

      monolis_solve_R(&mat, b, a);

      for(i = 0; i < n_node; i++){
        monolis_test_check_eq_R1("monolis_solver_parallel_R_test Clang", 1.0, a[i]);
      }
    }
  }

  monolis_finalize(&mat);
}

int main()
{
  monolis_global_initialize();

  monolis_solver_parallel_R_test();

  //monolis_solver_parallel_C_test()

  monolis_global_finalize();

  return 0;
}
