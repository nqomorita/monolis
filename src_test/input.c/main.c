#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <stdbool.h>
#include "monolis.h"

void input_coef(
  char *fname,
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
  char* fname;
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
  //monolis_input_node(fname, &n_node, node);

  fname = monolis_get_global_input_file_name("parted.0", "elem.dat");
  //monolis_input_elem(fname, &n_elem, &n_base, elem);

  if(monolis_mpi_get_global_comm_size() > 1){
    fname = monolis_get_global_input_file_name("parted.0", "elem.dat.id");
    monolis_input_global_id(fname, &n_id, &global_eid);
  } else {
    global_eid = monolis_alloc_I_1d(global_eid, n_elem);
    for(i = 0; i < n_elem; i++){
      global_eid[i] = i;
    }
  }

  monolis_initialize(&mat);

  monolis_get_nonzero_pattern_R(&mat, n_node, 2, 1, n_elem, elem);

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

  monolis_set_method(&mat, iter);
  monolis_set_precond(&mat, prec);

  monolis_solve_R(&mat, b, a);

  monolis_finalize(&mat);
}

int main()
{
  monolis_global_initialize();

  monolis_solver_parallel_R_test();

  //monolis_solver_parallel_C_test();

  monolis_global_finalize();

  return 0;
}
