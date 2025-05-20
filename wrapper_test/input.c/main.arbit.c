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
    printf("%s cant open file\n", fname);
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
  int* global_nid;
  int* n_dof_list;
  int** elem;
  double val;
  double res_conv;
  double* coef;
  double* a;
  double* b;
  double* c;
  double** node;
  double condition_number_lanczos = 0.0;

  fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.dat");
  monolis_input_mesh_node_c(fname, &n_node, &node);

  fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "elem.dat");
  monolis_input_mesh_elem_c(fname, &n_elem, &n_base, &elem);

  if(monolis_mpi_get_global_comm_size() > 1){
    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.dat.id");
    monolis_input_id_c(fname, &global_nid);
    
    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "elem.dat.id");
    monolis_input_id_c(fname, &global_eid);
  } else {
    global_nid = monolis_alloc_I_1d(global_nid, n_node);
    for(i = 0; i < n_node; i++){
      global_nid[i] = i; 
    }
    
    global_eid = monolis_alloc_I_1d(global_eid, n_elem);
    for(i = 0; i < n_elem; i++){
      global_eid[i] = i;
    }
  }

  monolis_initialize(&mat);

  monolis_com_initialize_by_parted_files(&com, monolis_mpi_get_global_comm(),
      MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.dat");

  // 可変自由度リストを作成
  n_dof_list = monolis_alloc_I_1d(n_dof_list, n_node);
  for(i = 0; i < n_node; i++){
    j = global_nid[i];
    if(j % 2 != 0){
      n_dof_list[i] = 1;
    } else {
      n_dof_list[i] = 2;
    }
  }

  // 可変自由度に対応した非ゼロパターン取得関数を使用
  monolis_get_nonzero_pattern_by_simple_mesh_V_R(&mat, n_node, 2, n_dof_list, n_elem, elem);

  fname = "coef.dat";
  input_coef(fname, &n_coef, &coef);

  for(i = 0; i < n_elem; i++){
    for(j = 0; j < 2; j++){
      eid[j] = elem[i][j];
    }

    val = coef[global_eid[i]]; 

    if(eid[0] == eid[1]){
      j = global_nid[eid[0]]; 
      if(j % 2 != 0){
        monolis_add_scalar_to_sparse_matrix_R(&mat, eid[0], eid[1], 0, 0, val);
      } else {
        monolis_add_scalar_to_sparse_matrix_R(&mat, eid[0], eid[1], 0, 0, val);
        monolis_add_scalar_to_sparse_matrix_R(&mat, eid[0], eid[1], 1, 1, val);
        monolis_add_scalar_to_sparse_matrix_R(&mat, eid[0], eid[1], 0, 1, 0.25*val);
        monolis_add_scalar_to_sparse_matrix_R(&mat, eid[0], eid[1], 1, 0, 0.25*val);
      }
    } else {
      monolis_add_scalar_to_sparse_matrix_R(&mat, eid[0], eid[1], 0, 0, val);
      monolis_add_scalar_to_sparse_matrix_R(&mat, eid[1], eid[0], 0, 0, val);
    }
  }

  int n = mat.mat.n_dof_index[n_node];
  
  a = monolis_alloc_R_1d(a, n);
  b = monolis_alloc_R_1d(b, n);
  c = monolis_alloc_R_1d(c, n);

  for(i = 0; i < n; i++){
    a[i] = 1.0;
  }

  monolis_matvec_product_R(&mat, &com, a, c);

  monolis_set_maxiter(&mat, 1000);
  monolis_set_tolerance(&mat, 1.0e-10);
  monolis_show_timelog(&mat, true);
  monolis_show_iterlog(&mat, true);
  monolis_show_summary(&mat, true);

  for (iter = MONOLIS_ITER_CG; iter < MONOLIS_ITER_IDRS + 1; ++iter) {
    for (prec = MONOLIS_PREC_NONE; prec < MONOLIS_PREC_SOR + 1; ++prec) {
      // いくつかの組み合わせはスキップ
      if(iter == MONOLIS_ITER_PIPECG && prec == MONOLIS_PREC_SOR) continue;
      if(iter == MONOLIS_ITER_PIPECR && prec == MONOLIS_PREC_SOR) continue;
      if(iter == MONOLIS_ITER_PIPEBICGSTAB && prec == MONOLIS_PREC_SOR) continue;

      for(i = 0; i < n; i++){
        a[i] = 0.0;
        b[i] = c[i];
      }

      if(monolis_mpi_get_global_my_rank() == 0){
        printf("iter %d, prec %d\n", iter, prec);
      }

      monolis_set_method(&mat, iter);
      monolis_set_precond(&mat, prec);

      monolis_solve_R(&mat, &com, b, a);

      monolis_mpi_global_barrier();

      for(i = 0; i < n; i++){
        monolis_test_check_eq_R1("monolis_solver_parallel_R_test", 1.0, a[i]);
      }

      monolis_get_converge_iter(&mat, &iter_conv);
      if(iter_conv <= 1){
        monolis_test_assert_fail("monolis_solver_parallel_R_test", "conv iter is less than 1");
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

  // 固有値計算用に行列の対角成分を変更
  for(i = 0; i < n_node; i++){
    j = global_nid[i];
    if(j % 2 != 0){
      monolis_add_scalar_to_sparse_matrix_R(&mat, i, i, 0, 0, (double)j);
    } else {
      monolis_add_scalar_to_sparse_matrix_R(&mat, i, i, 0, 0, (double)j);
      monolis_add_scalar_to_sparse_matrix_R(&mat, i, i, 1, 1, (double)j + 0.5);
    }
  }

  int n_get_eigen = 15;
  double* eig_val1;
  double* eig_val2;
  double** eig_mode1;
  double** eig_mode2;
  bool* is_bc;

  eig_val1 = monolis_alloc_R_1d(eig_val1, n_get_eigen);
  eig_val2 = monolis_alloc_R_1d(eig_val2, n_get_eigen);
  eig_mode1 = monolis_alloc_R_2d(eig_mode1, n, n_get_eigen);
  eig_mode2 = monolis_alloc_R_2d(eig_mode2, n, n_get_eigen);
  is_bc = (bool*)calloc(n, sizeof(bool));

  monolis_set_method(&mat, MONOLIS_ITER_CG);
  monolis_set_precond(&mat, MONOLIS_PREC_SOR);
  monolis_show_timelog(&mat, false);
  monolis_show_iterlog(&mat, false);
  monolis_show_summary(&mat, false);

  monolis_eigen_standard_lanczos_R(&mat, &com, &n_get_eigen, 1.0e-6, 100, eig_val1, eig_mode1, is_bc);

  monolis_eigen_inverted_standard_lanczos_R(&mat, &com, &n_get_eigen, 1.0e-6, 100, eig_val2, eig_mode2, is_bc);

  condition_number_lanczos = eig_val1[0]/eig_val1[14];  // C では 0 から始まるインデックス

  for (int i = 0; i < n_get_eigen; ++i) {
    j = n_get_eigen - i - 1;
    monolis_test_check_eq_R1("monolis_solver_parallel_R_test eig value", eig_val1[i], eig_val2[j]);
    for (int k = 0; k < n; ++k) {
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
  int* global_nid;
  int* global_eid;
  int* n_dof_list;
  int** elem;
  double res_conv;
  double r;
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
    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.dat.id");
    monolis_input_id_c(fname, &global_nid);
    
    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "elem.dat.id");
    monolis_input_id_c(fname, &global_eid);
  } else {
    global_nid = monolis_alloc_I_1d(global_nid, n_node);
    for(i = 0; i < n_node; i++){
      global_nid[i] = i; 
    }
    
    global_eid = monolis_alloc_I_1d(global_eid, n_elem);
    for(i = 0; i < n_elem; i++){
      global_eid[i] = i;
    }
  }

  monolis_initialize(&mat);

  monolis_com_initialize_by_parted_files(&com, monolis_mpi_get_global_comm(),
      MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.dat");

  // 可変自由度リストを作成
  n_dof_list = monolis_alloc_I_1d(n_dof_list, n_node);
  for(i = 0; i < n_node; i++){
    j = global_nid[i];
    if(j % 2 != 0){
      n_dof_list[i] = 1;
    } else {
      n_dof_list[i] = 2;
    }
  }

  // 可変自由度に対応した非ゼロパターン取得関数を使用
  monolis_get_nonzero_pattern_by_simple_mesh_V_C(&mat, n_node, 2, n_dof_list, n_elem, elem);

  fname = "coef.dat";
  input_coef(fname, &n_coef, &coef);

  for(i = 0; i < n_elem; i++){
    for(j = 0; j < 2; j++){
      eid[j] = elem[i][j];
    }

    r = coef[global_eid[i]];
    val = r + r*I;

    if(eid[0] == eid[1]){
      j = global_nid[eid[0]]; 
      if(j % 2 != 0){
        monolis_add_scalar_to_sparse_matrix_C(&mat, eid[0], eid[1], 0, 0, val);
      } else {
        monolis_add_scalar_to_sparse_matrix_C(&mat, eid[0], eid[1], 0, 0, val);
        monolis_add_scalar_to_sparse_matrix_C(&mat, eid[0], eid[1], 1, 1, val);
        monolis_add_scalar_to_sparse_matrix_C(&mat, eid[0], eid[1], 0, 1, 0.25*val);
        monolis_add_scalar_to_sparse_matrix_C(&mat, eid[0], eid[1], 1, 0, 0.25*val);
      }
    } else {
      monolis_add_scalar_to_sparse_matrix_C(&mat, eid[0], eid[1], 0, 0, val);
      monolis_add_scalar_to_sparse_matrix_C(&mat, eid[1], eid[0], 0, 0, val);
    }
  }

  // 自由度の合計数を取得
  int n = mat.mat.n_dof_index[n_node];
  
  
  a = monolis_alloc_C_1d(a, n);
  b = monolis_alloc_C_1d(b, n);
  c = monolis_alloc_C_1d(c, n);

  for(i = 0; i < n; i++){
    a[i] = 1.0 + 1.0*I;
  }

  monolis_matvec_product_C(&mat, &com, a, c);

  monolis_set_maxiter(&mat, 1000);
  monolis_set_tolerance(&mat, 1.0e-10);
  monolis_show_timelog(&mat, true);
  monolis_show_iterlog(&mat, true);
  monolis_show_summary(&mat, true);

  for (iter = MONOLIS_ITER_COCG; iter < MONOLIS_ITER_COCG + 1; ++iter) {
    for (prec = MONOLIS_PREC_NONE; prec < MONOLIS_PREC_SOR + 1; ++prec) {

      for(i = 0; i < n; i++){
        a[i] = 0.0 + 0.0*I;
        b[i] = c[i];
      }

      monolis_set_method(&mat, iter);
      monolis_set_precond(&mat, prec);

      monolis_solve_C(&mat, &com, b, a);

      monolis_mpi_global_barrier();

      for(i = 0; i < n; i++){
        monolis_test_check_eq_C1("monolis_solver_parallel_C_test", 1.0 + 1.0*I, a[i]);
      }

      monolis_get_converge_iter(&mat, &iter_conv);
      if(iter_conv <= 1){
        monolis_test_assert_fail("monolis_solver_parallel_C_test", "conv iter is less than 1");
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

void monolis_condition_number_R_test(){
  MONOLIS mat;
  MONOLIS_COM com;
  const char* fname;
  int n_node, n_elem, n_base, n_id, n_coef;
  int eid[2], i, j, k;
  int* global_nid;
  int* global_eid;
  int* n_dof_list;
  int** elem;
  double val;
  double** dense;
  double* coef;
  double** node;
  double rmax, rmin, condition_number;

  fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.dat");
  monolis_input_mesh_node_c(fname, &n_node, &node);

  fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "elem.dat");
  monolis_input_mesh_elem_c(fname, &n_elem, &n_base, &elem);

  if(monolis_mpi_get_global_comm_size() > 1){
    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.dat.id");
    monolis_input_id_c(fname, &global_nid);
    
    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "elem.dat.id");
    monolis_input_id_c(fname, &global_eid);
  } else {
    global_nid = monolis_alloc_I_1d(global_nid, n_node);
    for(i = 0; i < n_node; i++){
      global_nid[i] = i + 1;
    }
    
    global_eid = monolis_alloc_I_1d(global_eid, n_elem);
    for(i = 0; i < n_elem; i++){
      global_eid[i] = i + 1;
    }
  }

  monolis_initialize(&mat);

  monolis_com_initialize_by_parted_files(&com, monolis_mpi_get_global_comm(),
      MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.dat");

  // 可変自由度リストを作成
  n_dof_list = monolis_alloc_I_1d(n_dof_list, n_node);
  for(i = 0; i < n_node; i++){
    j = global_nid[i];
    if(j % 2 != 0){
      n_dof_list[i] = 1;
    } else {
      n_dof_list[i] = 2;
    }
  }

  // 可変自由度に対応した非ゼロパターン取得関数を使用
  monolis_get_nonzero_pattern_by_simple_mesh_V_R(&mat, n_node, 2, n_dof_list, n_elem, elem);

  fname = "coef.dat";
  input_coef(fname, &n_coef, &coef);

  for(i = 0; i < n_elem; i++){
    for(j = 0; j < 2; j++){
      eid[j] = elem[i][j];
    }

    val = coef[global_eid[i] - 1];

    if(eid[0] == eid[1]){
      j = global_nid[eid[0] - 1];
      if(j % 2 != 0){
        monolis_add_scalar_to_sparse_matrix_R(&mat, eid[0], eid[1], 0, 0, val);
      } else {
        monolis_add_scalar_to_sparse_matrix_R(&mat, eid[0], eid[1], 0, 0, val);
        monolis_add_scalar_to_sparse_matrix_R(&mat, eid[0], eid[1], 1, 1, val);
        monolis_add_scalar_to_sparse_matrix_R(&mat, eid[0], eid[1], 0, 1, 0.25*val);
        monolis_add_scalar_to_sparse_matrix_R(&mat, eid[0], eid[1], 1, 0, 0.25*val);
      }
    } else {
      monolis_add_scalar_to_sparse_matrix_R(&mat, eid[0], eid[1], 0, 0, val);
      monolis_add_scalar_to_sparse_matrix_R(&mat, eid[1], eid[0], 0, 0, val);
    }
  }

  // 疎行列を密行列に変換
  int n = mat.mat.n_dof_index[n_node];
  //monolis_convert_sparse_matrix_to_dense_matrix_R(&mat, &com, &dense);

  // if(monolis_mpi_get_global_comm_size() == 1){
  //   for(i = 0; i < 10; i++){
  //     //monolis_test_check_eq_R1("monolis_convert_sparse_matrix_to_dense_matrix_R test 1", dense[i][i], 4.0);
  //   }
  //   for(i = 0; i < 9; i++){
  //     //monolis_test_check_eq_R1("monolis_convert_sparse_matrix_to_dense_matrix_R test 1a", dense[i][i+1], 1.0);
  //     //monolis_test_check_eq_R1("monolis_convert_sparse_matrix_to_dense_matrix_R test 1b", dense[i+1][i], 1.0);
  //   }
  // }

  // ノード番号に基づいて行列の要素を追加 (固有値計算のため)
  for(i = 0; i < n_node; i++){
    j = global_nid[i];
    if(j % 2 != 0){
      monolis_add_scalar_to_sparse_matrix_R(&mat, i, i, 0, 0, (double)j);
    } else {
      monolis_add_scalar_to_sparse_matrix_R(&mat, i, i, 0, 0, (double)j);
      monolis_add_scalar_to_sparse_matrix_R(&mat, i, i, 1, 1, (double)j + 0.5);
    }
  }

  // 条件数を計算
  monolis_get_condition_number_R(&mat, &com, &rmax, &rmin);
  condition_number = rmax/rmin;

  printf("Condition Number: %.15f\n", condition_number);

  monolis_finalize(&mat);
}

int main()
{
  monolis_global_initialize();

  monolis_solver_parallel_R_test();

  //monolis_solver_parallel_C_test();

  //monolis_condition_number_R_test();

  monolis_global_finalize();

  return 0;
}
