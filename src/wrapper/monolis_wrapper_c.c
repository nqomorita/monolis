#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include "monolis.h"
#include "metis.h"

void monolis_param_set_method   (MONOLIS* mat, int    flag) {mat->prm.method    = flag;}
void monolis_param_set_precond  (MONOLIS* mat, int    flag) {mat->prm.precond   = flag;}
void monolis_param_set_maxiter  (MONOLIS* mat, int    flag) {mat->prm.maxiter   = flag;}
void monolis_param_set_tolerance(MONOLIS* mat, double flag) {mat->prm.tol       = flag;}
void monolis_param_set_performance_measurement (MONOLIS* mat, bool   flag) {mat->prm.is_measurement = flag;}
void monolis_param_set_check_diag (MONOLIS* mat, bool flag) {mat->prm.is_check_diag = flag;}
void monolis_param_set_init_x (MONOLIS* mat, bool flag) {mat->prm.is_init_x = flag;}

void monolis_param_set_show_iterlog (MONOLIS* mat, bool   flag) {mat->prm.show_iterlog = flag;}
void monolis_param_set_show_timelog (MONOLIS* mat, bool   flag) {mat->prm.show_timelog = flag;}
void monolis_param_set_show_timelog_statistics (MONOLIS* mat, bool   flag) {mat->prm.show_timelog_statistics = flag;}
void monolis_param_set_show_summary (MONOLIS* mat, bool   flag) {mat->prm.show_summary = flag;}

void monolis_get_time_solver       (MONOLIS* mat, double* value) {*value = mat->prm.tsol;};
void monolis_get_time_preparing    (MONOLIS* mat, double* value) {*value = mat->prm.tprep;};
void monolis_get_time_spmv         (MONOLIS* mat, double* value) {*value = mat->prm.tspmv;};
void monolis_get_time_inner_product(MONOLIS* mat, double* value) {*value = mat->prm.tdotp;};
void monolis_get_time_precondition (MONOLIS* mat, double* value) {*value = mat->prm.tprec;};
void monolis_get_time_comm_inner_product (MONOLIS* mat, double* value) {*value = mat->prm.tcomm_dotp;};
void monolis_get_time_comm_spmv    (MONOLIS* mat, double* value) {*value = mat->prm.tcomm_spmv;};

void monolis_get_converge_iter     (MONOLIS* mat, int* value) {*value = mat->prm.curiter;};
void monolis_get_converge_residual (MONOLIS* mat, double* value) {*value = mat->prm.curresid;};

/* to be deleted >>>*/
void monolis_set_method   (MONOLIS* mat, int    flag) {mat->prm.method    = flag;}
void monolis_set_precond  (MONOLIS* mat, int    flag) {mat->prm.precond   = flag;}
void monolis_set_maxiter  (MONOLIS* mat, int    flag) {mat->prm.maxiter   = flag;}
void monolis_set_tolerance(MONOLIS* mat, double flag) {mat->prm.tol       = flag;}
void monolis_set_performance_measurement (MONOLIS* mat, bool   flag) {mat->prm.is_measurement = flag;}
void monolis_set_check_diag (MONOLIS* mat, bool flag) {mat->prm.is_check_diag = flag;}
void monolis_set_init_x (MONOLIS* mat, bool flag) {mat->prm.is_init_x = flag;}

void monolis_show_iterlog (MONOLIS* mat, bool   flag) {mat->prm.show_iterlog = flag;}
void monolis_show_timelog (MONOLIS* mat, bool   flag) {mat->prm.show_timelog = flag;}
void monolis_show_timelog_statistics (MONOLIS* mat, bool   flag) {mat->prm.show_timelog_statistics = flag;}
void monolis_show_summary (MONOLIS* mat, bool   flag) {mat->prm.show_summary = flag;}
/* to be deleted <<<*/

/* initializer */
FILE* monolis_open_comm_table(
    FILE*       fp,
    const char* input_file_dir,
    const char* filename,
    int         myrank)
{
  const char* dirname = "parted.0";
  const int BUFFER_SIZE = 1024;
  char fname[BUFFER_SIZE];
  snprintf(fname, BUFFER_SIZE, "%s/%s/%s.%d", input_file_dir, dirname, filename, myrank);

  fp = fopen(fname, "r");
  if( fp == NULL ){
    printf("** error: file name %s\n", fname);
    printf("** error: monolis_open_comm_table\n");
    exit(monolis_fail);
  }
  return fp;
}

void monolis_com_input_comm_table(
  MONOLIS* mat,
  const char* input_file_dir)
{
  FILE* fp;

  if(mat->com.commsize <= 1){
    mat->com.commsize = 1;
    mat->com.send_n_neib = 0;
    mat->com.recv_n_neib = 0;
    mat->com.send_neib_pe = (int*)calloc(1, sizeof(int));
    mat->com.send_index = (int*)calloc(2, sizeof(int));
    mat->com.send_item = (int*)calloc(1, sizeof(int));
    mat->com.recv_neib_pe = (int*)calloc(1, sizeof(int));
    mat->com.recv_index = (int*)calloc(2, sizeof(int));
    mat->com.recv_item = (int*)calloc(1, sizeof(int));
    mat->com.global_node_id = (int*)calloc(1, sizeof(int));
    mat->com.global_elem_id = (int*)calloc(1, sizeof(int));
    return;
  }

  int nneib, nitem, nin, ein;
  fp = monolis_open_comm_table(fp, input_file_dir, "monolis.send", mat->com.myrank);
    fscanf(fp, "%d %d", &nneib, &nitem);
    mat->com.send_n_neib = nneib;
    mat->com.send_neib_pe = (int*)calloc(nneib, sizeof(int));
    for(int i=0; i<nneib; i++){
      fscanf(fp, "%d", &(mat->com.send_neib_pe[i]));
    }

    mat->com.send_index = (int*)calloc(nneib+1, sizeof(int));
    mat->com.send_item = (int*)calloc(nitem, sizeof(int));
    for(int i=0; i<nneib+1; i++){
      fscanf(fp, "%d", &(mat->com.send_index[i]));
    }
    for(int i=0; i<nitem; i++){
      fscanf(fp, "%d", &(mat->com.send_item[i]));
    }
  fclose(fp);

  fp = monolis_open_comm_table(fp, input_file_dir, "monolis.recv", mat->com.myrank);
    fscanf(fp, "%d %d", &nneib, &nitem);
    mat->com.recv_n_neib = nneib;
    mat->com.recv_neib_pe = (int*)calloc(nneib, sizeof(int));
    for(int i=0; i<nneib; i++){
      fscanf(fp, "%d", &(mat->com.recv_neib_pe[i]));
    }

    mat->com.recv_index = (int*)calloc(nneib+1, sizeof(int));
    mat->com.recv_item = (int*)calloc(nitem, sizeof(int));
    for(int i=0; i<nneib+1; i++){
      fscanf(fp, "%d", &(mat->com.recv_index[i]));
    }
    for(int i=0; i<nitem; i++){
      fscanf(fp, "%d", &(mat->com.recv_item[i]));
    }

  fclose(fp);

  fp = monolis_open_comm_table(fp, input_file_dir, "node.id", mat->com.myrank);
    fscanf(fp, "%d", &nitem);
    mat->com.nnode = nitem;
    mat->com.global_node_id = (int*)calloc(nitem, sizeof(int));
    for(int i=0; i<nitem; i++){
      fscanf(fp, "%d %d %d", &nin, &nin, &(mat->com.global_node_id[i]));
    }
  fclose(fp);

  fp = monolis_open_comm_table(fp, input_file_dir, "connectivity.id", mat->com.myrank);
    fscanf(fp, "%d", &nitem);
    mat->com.nelem = nitem;
    mat->com.global_elem_id = (int*)calloc(nitem, sizeof(int));
    for(int i=0; i<nitem; i++){
      fscanf(fp, "%d %d %d", &nin, &nin, &(mat->com.global_elem_id[i]));
    }
  fclose(fp);

  fp = monolis_open_comm_table(fp, input_file_dir, "node.n_internal", mat->com.myrank);
    fscanf(fp, "%d", &nin);
    mat->com.internal_nnode = nin;
  fclose(fp);

  fp = monolis_open_comm_table(fp, input_file_dir, "connectivity.n_internal", mat->com.myrank);
    fscanf(fp, "%d", &nin);
    mat->com.internal_nelem = nin;
  fclose(fp);
}

void monolis_prm_initialize(
  MONOLIS*    mat)
{
    mat->prm.method = 1;
    mat->prm.precond = 1;
    mat->prm.maxiter = 1000;
    mat->prm.curiter = 0;
    mat->prm.ierr = -1;
    mat->prm.tol = 1.0e-8;
    mat->prm.curresid = 0.0;
    mat->prm.is_scaling = false;
    mat->prm.is_reordering = false;
    mat->prm.is_init_x = true;
    mat->prm.is_debug = false;
    mat->prm.is_check_diag = true;
    mat->prm.show_iterlog = false;
    mat->prm.show_timelog_statistics = false;
    mat->prm.show_timelog = false;
    mat->prm.show_summary = true;

    mat->prm.tsol  = 0.0;
    mat->prm.tprep = 0.0;
    mat->prm.tspmv = 0.0;
    mat->prm.tdotp = 0.0;
    mat->prm.tprec = 0.0;
    mat->prm.tcomm_dotp = 0.0;
    mat->prm.tcomm_spmv = 0.0;

}

void monolis_com_initialize(
  MONOLIS*    mat)
{
    mat->com.recv_n_neib = 0;
    mat->com.send_n_neib = 0;
    mat->com.internal_nnode = 0;
    mat->com.internal_nelem = 0;
    mat->com.myrank = monolis_get_global_myrank();
    mat->com.comm = monolis_get_global_comm();
    mat->com.commsize = monolis_get_global_commsize();
}

void monolis_mat_initialize(
  MONOLIS*    mat)
{
    mat->mat.A = NULL;
    mat->mat.X = NULL;
    mat->mat.B = NULL;
    mat->mat.index = NULL;
    mat->mat.item = NULL;
    mat->mat.indexR = NULL;
    mat->mat.itemR = NULL;
    mat->mat.permR = NULL;
    mat->mat.N = 0;
    mat->mat.NP = 0;
    mat->mat.NZ = 0;
    mat->mat.NDOF = 0;
}

void monolis_initialize(
  MONOLIS*    mat,
  const char* input_file_dir)
{
    monolis_prm_initialize(mat);
    monolis_com_initialize(mat);
    monolis_mat_initialize(mat);
    monolis_com_input_comm_table(mat, input_file_dir);
}

void monolis_finalize(
  MONOLIS* mat)
{
  monolis_free_int_1d(mat->com.send_neib_pe);
  monolis_free_int_1d(mat->com.send_index);
  monolis_free_int_1d(mat->com.send_item);
  monolis_free_int_1d(mat->com.recv_neib_pe);
  monolis_free_int_1d(mat->com.recv_index);
  monolis_free_int_1d(mat->com.recv_item);
  monolis_free_int_1d(mat->com.global_node_id);
  monolis_free_int_1d(mat->com.global_elem_id);

  monolis_free_double_1d(mat->mat.A);
  monolis_free_double_1d(mat->mat.X);
  monolis_free_double_1d(mat->mat.B);
  monolis_free_int_1d(mat->mat.index);
  monolis_free_int_1d(mat->mat.item);
  monolis_free_int_1d(mat->mat.indexR);
  monolis_free_int_1d(mat->mat.itemR);
  monolis_free_int_1d(mat->mat.permR);
}

void monolis_free_int_1d(
  int* array)
{
  if(array == NULL) return;
  free(array);
}

void monolis_free_double_1d(
  double* array)
{
  if(array == NULL) return;
  free(array);
}

/* mat clear */

void monolis_clear(
  MONOLIS* mat)
{
  monolis_clear_mat(mat);
  monolis_clear_rhs(mat);
  monolis_clear_solution(mat);
}

void monolis_clear_mat(
  MONOLIS* mat)
{
  int ndof2 = mat->mat.NDOF*mat->mat.NDOF;

  for(int i=0; i<(mat->mat.NZ*ndof2); i++) {
    mat->mat.A[i] = 0.0;
  }
}

void monolis_clear_rhs(
  MONOLIS* mat)
{
  int ndof = mat->mat.NDOF;

  for(int i=0; i<(mat->mat.NP*ndof); i++) {
    mat->mat.B[i] = 0.0;
  }
}

void monolis_clear_solution(
  MONOLIS* mat)
{
  int ndof = mat->mat.NDOF;

  for(int i=0; i<(mat->mat.NP*ndof); i++) {
    mat->mat.X[i] = 0.0;
  }
}

/* mat copy */

void monolis_copy_all(
  MONOLIS* in,
  MONOLIS* out)
{
  int ndof = in->mat.NDOF;
  int np = in->mat.NP;
  int nz = in->mat.NZ;

  monolis_copy_nonzero_pattern(in, out);

  for(int i=0; i<ndof*ndof*nz; i++) {
    out->mat.A[i] = in->mat.A[i];
  }

  for(int i=0; i<ndof*np; i++) {
    out->mat.X[i] = in->mat.X[i];
  }

  for(int i=0; i<ndof*np; i++) {
    out->mat.B[i] = in->mat.B[i];
  }
}

void monolis_copy_nonzero_pattern(
  MONOLIS* in,
  MONOLIS* out)
{
  int n = in->mat.N;
  int np = in->mat.NP;
  int ndof = in->mat.NDOF;
  int nz = in->mat.NZ;

  monolis_copy_param(in, out);

  out->mat.N = n;
  out->mat.NP = np;
  out->mat.NDOF = ndof;
  out->mat.NZ = nz;

  out->mat.A = (double*)calloc(ndof*ndof*nz, sizeof(double));
  for(int i=0; i<ndof*ndof*nz; i++) {
    out->mat.A[i] = 0.0;
  }

  out->mat.X = (double*)calloc(ndof*np, sizeof(double));
  for(int i=0; i<ndof*np; i++) {
    out->mat.X[i] = 0.0;
  }

  out->mat.B = (double*)calloc(ndof*np, sizeof(double));
  for(int i=0; i<ndof*np; i++) {
    out->mat.B[i] = 0.0;
  }

  out->mat.index = (int* )calloc(np+1, sizeof(int));
  for(int i=0; i<np+1; i++) {
    out->mat.index[i] = in->mat.index[i];
  }

  out->mat.item = (int*)calloc(nz, sizeof(int));
  for(int i=0; i<nz; i++) {
    out->mat.item[i] = in->mat.item[i];
  }

  out->mat.indexR = (int*)calloc(np+1, sizeof(int));
  for(int i=0; i<np+1; i++) {
    out->mat.indexR[i] = in->mat.indexR[i];
  }

  out->mat.itemR = (int*)calloc(nz, sizeof(int));
  for(int i=0; i<nz; i++) {
    out->mat.itemR[i] = in->mat.itemR[i];
  }

  out->mat.permR = (int*)calloc(nz, sizeof(int));
  for(int i=0; i<nz; i++) {
    out->mat.permR[i] = in->mat.permR[i];
  }
}

void monolis_copy_param(
  MONOLIS* in,
  MONOLIS* out)
{
  out->prm.method = in->prm.method;
  out->prm.precond = in->prm.precond;
  out->prm.maxiter = in->prm.maxiter;
  out->prm.curiter = in->prm.curiter;
  out->prm.ierr = in->prm.ierr;
  out->prm.tol = in->prm.tol;
  out->prm.curresid = in->prm.curresid;
  out->prm.is_scaling = in->prm.is_scaling;
  out->prm.is_check_diag = in->prm.is_check_diag;
  out->prm.is_reordering = in->prm.is_reordering;
  out->prm.is_init_x = in->prm.is_init_x;
  out->prm.is_debug = in->prm.is_debug;
  out->prm.show_iterlog = in->prm.show_iterlog;
  out->prm.show_timelog = in->prm.show_timelog;
  out->prm.show_timelog_statistics = in->prm.show_timelog_statistics;
  out->prm.show_summary = in->prm.show_summary;

  out->prm.tsol = in->prm.tsol;
  out->prm.tprep = in->prm.tprep;
  out->prm.tspmv = in->prm.tspmv;
  out->prm.tdotp = in->prm.tdotp;
  out->prm.tprec = in->prm.tprec;
  out->prm.tcomm_dotp = in->prm.tcomm_dotp;
  out->prm.tcomm_spmv = in->prm.tcomm_spmv;

  out->com.myrank = in->com.myrank;
  out->com.comm = in->com.comm;
  out->com.commsize = in->com.commsize;
  out->com.recv_n_neib = in->com.recv_n_neib;
  out->com.send_n_neib = in->com.send_n_neib;
}

const char* monolis_get_input_filename(
  const char* filename_body)
{
  int commsize;
  int myrank;
  int buf_size = 1000;
  char* filename;

  commsize = monolis_get_global_commsize();
  myrank = monolis_get_global_myrank();
  filename = (char*)malloc(sizeof(char)*buf_size);

  if(commsize > 1){
    snprintf(filename, buf_size, "%s/%s.%d", "parted.0", filename_body, myrank);
  } else {
    snprintf(filename, buf_size, "%s", filename_body);
  }
  return filename;
}

const char* monolis_get_output_filename(
  const char* filename_body)
{
  int commsize;
  int myrank;
  int buf_size = 1000;
  char* head;
  char* post;
  char ctmp[buf_size], body[buf_size];
  char* filename;

  commsize = monolis_get_global_commsize();
  myrank = monolis_get_global_myrank();
  filename = (char*)malloc(sizeof(char)*buf_size);

  if(commsize > 1){
    snprintf(body, buf_size, "%s", filename_body);
    head = strtok(body, ".");
    post = strtok(NULL, ".");
    snprintf(filename, buf_size, "%s.%d.%s", head, myrank, post);
  } else {
    snprintf(filename, buf_size, "%s", filename_body);
  }
  return filename;
}

void monolis_convert_mesh_to_connectivity(
  int      nelem,
  int      nbase_func,
  int**    elem,
  idx_t*  conn_index,
  idx_t*  con)
{
  for(int i=0; i<nelem+1; i++) {
    conn_index[i] = i*nbase_func;
  }

  for(int i=0; i<nelem; i++) {
    for(int j=0; j<nbase_func; j++) {
      con[nbase_func*i+j] = elem[i][j];
    }
  }
}

void monolis_convert_connectivity_to_nodal_graph(
  int      nnode,
  int      nelem,
  idx_t*   conn_index,
  idx_t*   con,
  idx_t** index,
  idx_t** item)
{
#ifdef WITH_METIS
  idx_t numflag = 0;
  int ierr = METIS_MeshToNodal(
    &nelem,
    &nnode,
    conn_index,
    con,
    &numflag,
    index,
    item);
#else

#endif
}

void monolis_get_nonzero_pattern_by_nodal_graph(
  MONOLIS* mat,
  int      nnode,
  int      ndof,
  int*     index,
  int*     item)
{
  mat->mat.N = nnode;
  mat->mat.NP = nnode;
  mat->mat.NDOF = ndof;
  mat->mat.X = (double*)calloc(ndof*nnode, sizeof(double));
  mat->mat.B = (double*)calloc(ndof*nnode, sizeof(double));

  mat->mat.index = (int* )calloc(nnode+1, sizeof(int));
  for(int i=1; i<nnode+1; i++) {
    mat->mat.index[i] = index[i] + i;
  }

  int nz = mat->mat.index[nnode];
  mat->mat.NZ = nz;
  mat->mat.A = (double*)calloc(ndof*ndof*nz, sizeof(double));
  mat->mat.item = (int*)calloc(nz, sizeof(int));

  for(int i=0; i<nnode; i++) {
    int jS = mat->mat.index[i];
    int jE = mat->mat.index[i+1];
    mat->mat.item[jS] = i+1;
    for(int j=jS+1; j<jE; j++){
      mat->mat.item[j] = item[j-i-1] + 1;
    }
    monolis_qsort_int(
      &(mat->mat.item[jS]),
      1,
      jE-jS);
  }

  mat->mat.indexR = (int*)calloc(nnode+1, sizeof(int));
  mat->mat.itemR = (int*)calloc(nz, sizeof(int));
  mat->mat.permR = (int*)calloc(nz, sizeof(int));
  monolis_get_CRR_format(
    nnode,
    nz,
    mat->mat.index,
    mat->mat.item,
    mat->mat.indexR,
    mat->mat.itemR,
    mat->mat.permR);
}

void monolis_get_nonzero_pattern(
  MONOLIS* mat,
  int      nnode,
  int      nbase_func,
  int      ndof,
  int      nelem,
  int**    elem)
{
  idx_t* conn_index;
  idx_t* con;
  idx_t* index;
  idx_t* item;

  conn_index = (idx_t*)calloc(nelem+1, sizeof(idx_t));
  con = (idx_t*)calloc(nelem*nbase_func, sizeof(idx_t));

  monolis_convert_mesh_to_connectivity(
    nelem,
    nbase_func,
    elem,
    conn_index,
    con);

  monolis_convert_connectivity_to_nodal_graph(
    nnode,
    nelem,
    conn_index,
    con,
    &index,
    &item);

  monolis_get_nonzero_pattern_by_nodal_graph(
    mat,
    nnode,
    ndof,
    index,
    item);

  free(conn_index);
  free(con);
  free(index);
  free(item);
}

void monolis_add_scalar_to_sparse_matrix(
  MONOLIS* mat,
  double   val,
  int      i,
  int      j,
  int      submat_i,
  int      submat_j)
{
  int nnode = mat->mat.NP;
  int ndof = mat->mat.NDOF;
  int nz = mat->mat.NZ;

  monolis_add_scalar_to_sparse_matrix_c_main(
    nnode,
    nz,
    ndof,
    mat->mat.index,
    mat->mat.item,
    mat->mat.A,
    i,
    j,
    submat_i,
    submat_j,
    val);
}

void monolis_get_scalar_from_sparse_matrix(
  MONOLIS* mat,
  int      i,
  int      j,
  int      submat_i,
  int      submat_j,
  double*  val,
  bool*    is_find)
{
  int n = mat->mat.NP;
  int nz = mat->mat.NZ;
  int ndof = mat->mat.NDOF;
  int is_find_t = 0;

  monolis_get_scalar_from_sparse_matrix_c_main(
    n,
    nz,
    ndof,
    mat->mat.index,
    mat->mat.item,
    mat->mat.A,
    i,
    j,
    submat_i,
    submat_j,
    val,
    &is_find_t);

  *is_find = false;
  if(is_find_t == 1){
    *is_find = true;
  }
}

/* set BCSR information */
void monolis_set_matrix_BCSR(
  MONOLIS* mat,
  int      n,
  int      np,
  int      ndof,
  int      nz,
  double*  A,
  int*     index,
  int*     item){

  mat->mat.N = n;
  mat->mat.NP = np;
  mat->mat.NDOF = ndof;
  mat->mat.NZ = nz;
  mat->mat.A = (double*)calloc(ndof*ndof*nz, sizeof(double));
  mat->mat.X = (double*)calloc(ndof*np, sizeof(double));
  mat->mat.B = (double*)calloc(ndof*np, sizeof(double));
  mat->mat.index = (int* )calloc(np+1, sizeof(int));
  mat->mat.item = (int*)calloc(nz, sizeof(int));

  int i;
  for(i = 0; i < np + 1; i++) {
    mat->mat.index[i] = index[i];
  }

  for(i = 0; i < nz; i++) {
    mat->mat.item[i] = item[i] + 1;
  }

  for(i = 0; i < ndof*ndof*nz; i++) {
    mat->mat.A[i] = A[i];
  }

  mat->mat.indexR = (int*)calloc(np+1, sizeof(int));
  mat->mat.itemR = (int*)calloc(nz, sizeof(int));
  mat->mat.permR = (int*)calloc(nz, sizeof(int));

  monolis_get_CRR_format(
    np,
    nz,
    mat->mat.index,
    mat->mat.item,
    mat->mat.indexR,
    mat->mat.itemR,
    mat->mat.permR);
}

void monolis_set_matrix_BCSR_mat_val(
  MONOLIS* mat,
  int      ndof,
  int      nz,
  double*  A){

  int i;
  for(i = 0; i < ndof*ndof*nz; i++) {
    mat->mat.A[i] = A[i];
  }
}


void monolis_set_Dirichlet_bc(
  MONOLIS* mat,
  double*  b,
  int      node_id,
  int      ndof_bc,
  double   val)
{
  int nnode = mat->mat.NP;
  int ndof = mat->mat.NDOF;
  int nz = mat->mat.NZ;

  monolis_set_Dirichlet_bc_c_main(
    nnode,
    nz,
    ndof,
    mat->mat.index,
    mat->mat.item,
    mat->mat.indexR,
    mat->mat.itemR,
    mat->mat.permR,
    mat->mat.A,
    b,
    node_id,
    ndof_bc,
    val);
}

/* get communication table */
void monolis_com_get_comm_table(
  MONOLIS* mat,
  int      n,
  int      np,
  int*     glonal_node_id){
  int n_neib_recv;
  int n_recv_item;
  int n_neib_send;
  int n_send_item;

  monolis_com_get_comm_table_analysis_c_main(
    n,
    np,
    glonal_node_id,
    &n_neib_recv,
    &n_recv_item,
    &n_neib_send,
    &n_send_item,
    mat->com.comm);

  mat->com.recv_n_neib = n_neib_recv;
  mat->com.recv_neib_pe = (int*)calloc(n_neib_recv, sizeof(int));
  mat->com.recv_index = (int*)calloc(n_neib_recv + 1, sizeof(int));
  mat->com.recv_item = (int*)calloc(n_recv_item, sizeof(int));

  mat->com.send_n_neib = n_neib_send;
  mat->com.send_neib_pe = (int*)calloc(n_neib_send, sizeof(int));
  mat->com.send_index = (int*)calloc(n_neib_send + 1, sizeof(int));
  mat->com.send_item = (int*)calloc(n_send_item, sizeof(int));

  monolis_com_get_comm_table_set_c_main(
    n,
    np,
    glonal_node_id,
    mat->com.comm,
    mat->com.recv_n_neib,
    n_recv_item,
    mat->com.recv_neib_pe,
    mat->com.recv_index,
    mat->com.recv_item,
    mat->com.send_n_neib,
    n_send_item,
    mat->com.send_neib_pe,
    mat->com.send_index,
    mat->com.send_item);
}

void monolis_matvec_product(
  MONOLIS* mat,
  double*  x,
  double*  y)
{
  int n = mat->mat.N;
  if(mat->com.commsize > 1 ) n = mat->com.internal_nnode;
  int np = mat->mat.NP;
  int ndof = mat->mat.NDOF;
  int nz = mat->mat.NZ;
  int recv_nitem = mat->com.recv_index[mat->com.recv_n_neib];
  int send_nitem = mat->com.send_index[mat->com.send_n_neib];

  monolis_matvec_product_c_main(
    /* mat */
    n,
    np,
    nz,
    ndof,
    mat->mat.A,
    x,
    y,
    mat->mat.index,
    mat->mat.item,
    /* comm */
    mat->com.myrank,
    mat->com.comm,
    mat->com.commsize,
    mat->com.recv_n_neib,
    recv_nitem,
    mat->com.recv_neib_pe,
    mat->com.recv_index,
    mat->com.recv_item,
    mat->com.send_n_neib,
    send_nitem,
    mat->com.send_neib_pe,
    mat->com.send_index,
    mat->com.send_item);
}

void monolis_inner_product(
  MONOLIS* mat,
  double*  x,
  double*  y,
  double*  sum)
{
  int n = mat->mat.NP;
  int ndof = mat->mat.NDOF;
  int comm = mat->com.comm;

  monolis_inner_product_c_main(n, ndof, x, y, sum, comm);
}

double monolis_allreduce_double_scalar(
  MONOLIS* mat,
  int      tag,
  double   val)
{
  val = monolis_allreduce_double_scalar_c_main(val, tag, mat->com.comm);
  return val;
}

void monolis_solve(
  MONOLIS* mat,
  double*  b,
  double*  x)
{
  int n = mat->mat.N;
  if(mat->com.commsize > 1 ) n = mat->com.internal_nnode;
  int np = mat->mat.NP;
  int ndof = mat->mat.NDOF;
  int nz = mat->mat.NZ;
  int iterlog = 0;
  int timelog = 0;
  int timelog_statistics = 0;
  int summary = 0;
  int is_check_diag = 1;
  int is_measurement = 0;
  int is_init_x = 1;
  int recv_nitem = mat->com.recv_index[mat->com.recv_n_neib];
  int send_nitem = mat->com.send_index[mat->com.send_n_neib];
  int curiter;
  double curresid;
  double time[7];

  if(mat->prm.show_iterlog) iterlog = 1;
  if(mat->prm.show_timelog) timelog = 1;
  if(mat->prm.show_timelog_statistics) timelog_statistics = 1;
  if(mat->prm.show_summary) summary = 1;
  if(!mat->prm.is_check_diag) is_check_diag = 0;
  if(!mat->prm.is_init_x) is_init_x = 0;
  if(mat->prm.is_measurement) is_measurement = 1;

  monolis_solve_c_main(
    /* mat */
    n,
    np,
    nz,
    ndof,
    mat->mat.A,
    x,
    b,
    mat->mat.index,
    mat->mat.item,
    /* comm */
    mat->com.myrank,
    mat->com.comm,
    mat->com.commsize,
    mat->com.recv_n_neib,
    recv_nitem,
    mat->com.recv_neib_pe,
    mat->com.recv_index,
    mat->com.recv_item,
    mat->com.send_n_neib,
    send_nitem,
    mat->com.send_neib_pe,
    mat->com.send_index,
    mat->com.send_item,
    /* parameter */
    mat->prm.method,
    mat->prm.precond,
    mat->prm.maxiter,
    mat->prm.tol,
    iterlog,
    timelog,
    timelog_statistics,
    summary,
    is_check_diag,
    is_measurement,
    is_init_x,
    &curiter,
    &curresid,
    time);

    mat->prm.tsol = time[0];
    mat->prm.tprep = time[1];
    mat->prm.tspmv = time[2];
    mat->prm.tdotp = time[3];
    mat->prm.tprec = time[4];
    mat->prm.tcomm_dotp = time[5];
    mat->prm.tcomm_spmv = time[6];
    mat->prm.curiter = curiter;
    mat->prm.curresid = curresid;
}


/* eigen solver (inverted Lnaczos)*/
void monolis_eigen_inverted_standard_lanczos(
  MONOLIS* mat,
  int      n_get_eigen,
  double   ths,
  int      maxiter,
  double*  eigen_value,
  double** eigen_mode,
  bool*    is_Dirichlet_bc)
{
  int n = mat->mat.N;
  if(mat->com.commsize > 1 ) n = mat->com.internal_nnode;
  int np = mat->mat.NP;
  int ndof = mat->mat.NDOF;
  int nz = mat->mat.NZ;
  int iterlog = 0;
  int timelog = 0;
  int timelog_statistics = 0;
  int summary = 0;
  int is_check_diag = 1;
  int is_measurement = 0;
  int is_init_x = 1;
  int recv_nitem = mat->com.recv_index[mat->com.recv_n_neib];
  int send_nitem = mat->com.send_index[mat->com.send_n_neib];
  int curiter, i, j;
  int* is_Dirichlet_bc_int;
  double* eigen_mode_tmp;
  double curresid;
  double time[7];

  eigen_mode_tmp = (double*)calloc(np*ndof*n_get_eigen, sizeof(double));
  is_Dirichlet_bc_int = (int*)calloc(np*ndof, sizeof(int));

  for(int i = 0; i < np*ndof; i++){
    if (is_Dirichlet_bc[i]) {
      is_Dirichlet_bc_int[i] = 1;
    } else {
      is_Dirichlet_bc_int[i] = 0;
    }
  }

  if(mat->prm.show_iterlog) iterlog = 1;
  if(mat->prm.show_timelog) timelog = 1;
  if(mat->prm.show_timelog_statistics) timelog_statistics = 1;
  if(mat->prm.show_summary) summary = 1;
  if(!mat->prm.is_check_diag) is_check_diag = 0;
  if(!mat->prm.is_init_x) is_init_x = 0;
  if(mat->prm.is_measurement) is_measurement = 1;

  monolis_eigen_inverted_standard_lanczos_c_main(
    n,
    np,
    nz,
    ndof,
    mat->mat.A,
    mat->mat.index,
    mat->mat.item,
    /* comm */
    mat->com.myrank,
    mat->com.comm,
    mat->com.commsize,
    mat->com.recv_n_neib,
    recv_nitem,
    mat->com.recv_neib_pe,
    mat->com.recv_index,
    mat->com.recv_item,
    mat->com.send_n_neib,
    send_nitem,
    mat->com.send_neib_pe,
    mat->com.send_index,
    mat->com.send_item,
    /* parameter */
    mat->prm.method,
    mat->prm.precond,
    mat->prm.maxiter,
    mat->prm.tol,
    iterlog,
    timelog,
    timelog_statistics,
    summary,
    is_check_diag,
    is_measurement,
    is_init_x,
    &curiter,
    &curresid,
    time,
    n_get_eigen,
    ths,
    maxiter,
    eigen_value,
    eigen_mode_tmp,
    is_Dirichlet_bc_int);

    mat->prm.tsol = time[0];
    mat->prm.tprep = time[1];
    mat->prm.tspmv = time[2];
    mat->prm.tdotp = time[3];
    mat->prm.tprec = time[4];
    mat->prm.tcomm_dotp = time[5];
    mat->prm.tcomm_spmv = time[6];
    mat->prm.curiter = curiter;
    mat->prm.curresid = curresid;

    for(int i = 0; i < n_get_eigen; i++){
      for(int j = 0; j < np*ndof; j++){
        eigen_mode[i][j] = eigen_mode_tmp[np*ndof*i + j];
      }
    }

    free(eigen_mode_tmp);
    free(is_Dirichlet_bc_int);
}

/* eigen solver (forwared Lnaczos) */
void monolis_eigen_standard_lanczos(
  MONOLIS* mat,
  int      n_get_eigen,
  double   ths,
  int      maxiter,
  double*  eigen_value,
  double** eigen_mode,
  bool*    is_Dirichlet_bc)
{
  int n = mat->mat.N;
  if(mat->com.commsize > 1 ) n = mat->com.internal_nnode;
  int np = mat->mat.NP;
  int ndof = mat->mat.NDOF;
  int nz = mat->mat.NZ;
  int iterlog = 0;
  int timelog = 0;
  int timelog_statistics = 0;
  int summary = 0;
  int is_check_diag = 1;
  int is_measurement = 0;
  int is_init_x = 1;
  int recv_nitem = mat->com.recv_index[mat->com.recv_n_neib];
  int send_nitem = mat->com.send_index[mat->com.send_n_neib];
  int curiter, i, j;
  int* is_Dirichlet_bc_int;
  double* eigen_mode_tmp;
  double curresid;
  double time[7];

  eigen_mode_tmp = (double*)calloc(np*ndof*n_get_eigen, sizeof(double));
  is_Dirichlet_bc_int = (int*)calloc(np*ndof, sizeof(int));

  for(int i = 0; i < np*ndof; i++){
    if (is_Dirichlet_bc[i]) {
      is_Dirichlet_bc_int[i] = 1;
    } else {
      is_Dirichlet_bc_int[i] = 0;
    }
  }

  if(mat->prm.show_iterlog) iterlog = 1;
  if(mat->prm.show_timelog) timelog = 1;
  if(mat->prm.show_timelog_statistics) timelog_statistics = 1;
  if(mat->prm.show_summary) summary = 1;
  if(!mat->prm.is_check_diag) is_check_diag = 0;
  if(!mat->prm.is_init_x) is_init_x = 0;
  if(mat->prm.is_measurement) is_measurement = 1;

  monolis_eigen_standard_lanczos_c_main(
    n,
    np,
    nz,
    ndof,
    mat->mat.A,
    mat->mat.index,
    mat->mat.item,
    /* comm */
    mat->com.myrank,
    mat->com.comm,
    mat->com.commsize,
    mat->com.recv_n_neib,
    recv_nitem,
    mat->com.recv_neib_pe,
    mat->com.recv_index,
    mat->com.recv_item,
    mat->com.send_n_neib,
    send_nitem,
    mat->com.send_neib_pe,
    mat->com.send_index,
    mat->com.send_item,
    /* parameter */
    mat->prm.method,
    mat->prm.precond,
    mat->prm.maxiter,
    mat->prm.tol,
    iterlog,
    timelog,
    timelog_statistics,
    summary,
    is_check_diag,
    is_measurement,
    is_init_x,
    &curiter,
    &curresid,
    time,
    n_get_eigen,
    ths,
    maxiter,
    eigen_value,
    eigen_mode_tmp,
    is_Dirichlet_bc_int);

    mat->prm.tsol = time[0];
    mat->prm.tprep = time[1];
    mat->prm.tspmv = time[2];
    mat->prm.tdotp = time[3];
    mat->prm.tprec = time[4];
    mat->prm.tcomm_dotp = time[5];
    mat->prm.tcomm_spmv = time[6];
    mat->prm.curiter = curiter;
    mat->prm.curresid = curresid;

    for(int i = 0; i < n_get_eigen; i++){
      for(int j = 0; j < np*ndof; j++){
        eigen_mode[i][j] = eigen_mode_tmp[np*ndof*i + j];
      }
    }

    free(eigen_mode_tmp);
    free(is_Dirichlet_bc_int);
}

void monolis_get_internal_node_number(
  MONOLIS* mat,
  int*     nnode_internal)
{
  if (mat->com.commsize > 1) {
    *nnode_internal = mat->com.internal_nnode;
  } else {
    *nnode_internal = mat->mat.N;
  }
}

void monolis_get_internal_elem_1d_bool(
  MONOLIS* mat,
  int      nelem,
  bool*    is_internal_elem)
{
  if (mat->com.commsize > 1) {
    for(int i=0; i<nelem; i++){
      is_internal_elem[i] = false;
    }
    for(int i=0; i<mat->com.internal_nelem; i++){
      is_internal_elem[i] = true;
    }
  } else {
    for(int i=0; i<nelem; i++){
      is_internal_elem[i] = true;
    }
  }
}

void monolis_barrier(
  MONOLIS* mat)
{
  monolis_barrier_c_main(mat->com.comm);
}
