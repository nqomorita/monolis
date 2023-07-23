#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <complex.h>
#include <math.h>
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
  MONOLIS_COM com;
  const char* fname;
  int n_node, n_elem, n_base, n_id, n_coef;
  int eid[2], iter, prec, i, j, k, iter_conv;
  int* global_eid;
  int** elem;
  double val;
  double res_conv;
  double* coef;
  double* a;
  double* b;
  double* c;
  double** node;

  fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.dat");
  monolis_input_mesh_node_c(fname, &n_node, &node);

  fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "elem.dat");
  monolis_input_mesh_elem_c(fname, &n_elem, &n_base, &elem);

  if(monolis_mpi_get_global_comm_size() > 1){
    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "elem.dat.id");
    monolis_input_id_c(fname, &global_eid);
  } else {
    global_eid = monolis_alloc_I_1d(global_eid, n_elem);
    for(i = 0; i < n_elem; i++){
      global_eid[i] = i;
    }
  }

  monolis_initialize(&mat);

  monolis_com_initialize_by_parted_files(&com, monolis_mpi_get_global_comm(),
      MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.dat");

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

  monolis_matvec_product_R(&mat, &com, a, c);

  monolis_set_maxiter(&mat, 1000);
  monolis_set_tolerance(&mat, 1.0e-10);
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

      monolis_solve_R(&mat, &com, b, a);

      monolis_mpi_global_barrier();

      for(i = 0; i < n_node; i++){
        monolis_test_check_eq_R1("monolis_solver_parallel_R_test Clang", 1.0, a[i]);
      }

      monolis_get_converge_iter(&mat, &iter_conv);
      if(iter_conv <= 1){
        monolis_test_assert_fail("monolis_solver_parallel_R_test Clang", "conv iter is less than 1");
      }

      monolis_get_converge_residual(&mat, &res_conv);
      if(res_conv > 1.0e-10){
        monolis_test_assert_fail("monolis_solver_parallel_R_test", "residual is greater than ths");
      }

      monolis_mpi_global_barrier();
    }
  }


  /* eigen solver */
  monolis_std_log_string("monolis_solver_parallel_test eigen");

  int n_get_eigen = 10;
  double* eig_val1;
  double* eig_val2;
  double** eig_mode1;
  double** eig_mode2;
  bool* is_bc;

  eig_val1 = monolis_alloc_R_1d(eig_val1, n_get_eigen);
  eig_val2 = monolis_alloc_R_1d(eig_val2, n_get_eigen);
  eig_mode1 = monolis_alloc_R_2d(eig_mode1, n_get_eigen, n_node);
  eig_mode2 = monolis_alloc_R_2d(eig_mode2, n_get_eigen, n_node);
  is_bc = (bool*)calloc(n_node, sizeof(bool));;

  monolis_set_method(&mat, MONOLIS_ITER_CG);
  monolis_set_precond(&mat, MONOLIS_PREC_SOR);
  monolis_show_timelog(&mat, false);
  monolis_show_iterlog(&mat, false);
  monolis_show_summary(&mat, false);

  monolis_eigen_standard_lanczos_R(&mat, &com, &n_get_eigen, 1.0e-6, 100, eig_val1, eig_mode1, is_bc);

  monolis_eigen_inverted_standard_lanczos_R(&mat, &com, &n_get_eigen, 1.0e-6, 100, eig_val2, eig_mode2, is_bc);

  for (int i = 0; i < n_get_eigen; ++i) {
    j = n_get_eigen - i - 1;
    monolis_test_check_eq_R1("monolis_solver_parallel_R_test eig value", eig_val1[i], eig_val2[j]);
    for (int k = 0; k < n_node; ++k) {
      monolis_test_check_eq_R1("monolis_solver_parallel_R_test eig mode", fabs(eig_mode1[i][k]), fabs(eig_mode2[j][k]));
    }
  }

  monolis_finalize(&mat);
}

void monolis_solver_parallel_C_test(){
  MONOLIS mat;
  MONOLIS_COM com;
  const char* fname;
  int n_node, n_elem, n_base, n_id, n_coef;
  int eid[2], iter, prec, i, j, iter_conv;
  int* global_eid;
  int** elem;
  double res_conv;
  double _Complex val;
  double* coef;
  double _Complex* a;
  double _Complex* b;
  double _Complex* c;
  double** node;

  fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.dat");
  monolis_input_mesh_node_c(fname, &n_node, &node);

  fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "elem.dat");
  monolis_input_mesh_elem_c(fname, &n_elem, &n_base, &elem);

  if(monolis_mpi_get_global_comm_size() > 1){
    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "elem.dat.id");
    monolis_input_id_c(fname, &global_eid);
  } else {
    global_eid = monolis_alloc_I_1d(global_eid, n_elem);
    for(i = 0; i < n_elem; i++){
      global_eid[i] = i;
    }
  }

  monolis_initialize(&mat);

  monolis_com_initialize_by_parted_files(&com, monolis_mpi_get_global_comm(),
      MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.dat");

  monolis_get_nonzero_pattern_by_simple_mesh_C(&mat, n_node, 2, 1, n_elem, elem);

  fname = "coef.dat";
  input_coef(fname, &n_coef, &coef);

  for(i = 0; i < n_elem; i++){
    for(j = 0; j < 2; j++){
      eid[j] = elem[i][j];
    }

    val = coef[global_eid[i]] + coef[global_eid[i]]*I;

    if(eid[0] == eid[1]){
      monolis_add_scalar_to_sparse_matrix_C(&mat, eid[0], eid[1], 0, 0, val);
    } else{
      monolis_add_scalar_to_sparse_matrix_C(&mat, eid[0], eid[1], 0, 0, val);
      monolis_add_scalar_to_sparse_matrix_C(&mat, eid[1], eid[0], 0, 0, val);
    }
  }

  a = monolis_alloc_C_1d(a, n_node);
  b = monolis_alloc_C_1d(b, n_node);
  c = monolis_alloc_C_1d(c, n_node);

  for(i = 0; i < n_node; i++){
    a[i] = 1.0 + 1.0*I;
  }

  monolis_matvec_product_C(&mat, &com, a, c);

  monolis_set_maxiter(&mat, 1000);
  monolis_set_tolerance(&mat, 1.0e-10);
  monolis_show_timelog(&mat, true);
  monolis_show_iterlog(&mat, true);
  monolis_show_summary(&mat, true);

  for (iter = 9; iter < 10; ++iter) {
    for (prec = 0; prec < 3; ++prec) {

      for(i = 0; i < n_node; i++){
        a[i] = 0.0 + 0.0*I;
        b[i] = c[i];
      }

      monolis_set_method(&mat, iter);
      monolis_set_precond(&mat, prec);

      monolis_solve_C(&mat, &com, b, a);

      monolis_mpi_global_barrier();

      for(i = 0; i < n_node; i++){
        monolis_test_check_eq_C1("monolis_solver_parallel_C_test Clang", 1.0 + 1.0*I, a[i]);
      }

      monolis_get_converge_iter(&mat, &iter_conv);
      if(iter_conv <= 1){
        monolis_test_assert_fail("monolis_solver_parallel_C_test Clang", "conv iter is less than 1");
      }

      monolis_get_converge_residual(&mat, &res_conv);
      if(res_conv > 1.0e-10){
        monolis_test_assert_fail("monolis_solver_parallel_C_test", "residual is greater than ths");
      }

      monolis_mpi_global_barrier();
    }
  }

  monolis_finalize(&mat);
}

int main()
{
  monolis_global_initialize();

  monolis_solver_parallel_R_test();

  monolis_solver_parallel_C_test();

  monolis_global_finalize();

  return 0;
}
