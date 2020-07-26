#include <stdio.h>
#include <stdlib.h>
#include "monolis.h"
#include "metis.h"

/* header */
void monolis_add_sparse_matrix_c_main();
void monolis_set_Dirichlet_bc_c_main();
void monolis_solve_c_main();

/* body */

void monolis_initialize(
  MONOLIS* mat)
{
    // prm
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
    mat->prm.show_iterlog = true;
    mat->prm.show_time = true;
    mat->prm.show_summary = true;

    mat->prm.tsol  = 0.0;
    mat->prm.tprep = 0.0;
    mat->prm.tspmv = 0.0;
    mat->prm.tdotp = 0.0;
    mat->prm.tprec = 0.0;
    mat->prm.tcomm = 0.0;

    // mat
    mat->mat.N = 0;
    mat->mat.NP = 0;
    mat->mat.NZ = 0;
    mat->mat.NDOF = 0;

    // comm
    mat->com.myrank = 0;
    mat->com.comm = 0;
    mat->com.commsize = 1;
    mat->com.recv_n_neib = 0;
    mat->com.send_n_neib = 0;
}

void monolis_finalize(
  MONOLIS* mat)
{

}

void monolis_get_CRR_format(
  MONOLIS* mat,
  int      nnode,
  int      nbase_func,
  int      ndof,
  int      nelem,
  int**    elem)
{


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

void monolis_convert_connectivity_to_nodal(
  int      nnode,
  int      nelem,
  idx_t*   conn_index,
  idx_t*   con,
  idx_t** index,
  idx_t** item)
{
#ifdef WITH_METIS
  idx_t numflag = 0;
  int ierr = METIS_MeshToNodal(&nelem, &nnode, conn_index, con, &numflag, index, item);
#else

#endif
}

void monolis_get_nonzero_pattern_by_nodal(
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
  mat->mat.A = (double*)calloc(ndof*ndof*nz, sizeof(double));
  mat->mat.item = (int*)calloc(nz, sizeof(int));

  for(int i=0; i<nnode; i++) {
    int jS = mat->mat.index[i];
    int jE = mat->mat.index[i+1];
    mat->mat.item[jS] = i+1;
    for(int j=jS+1; j<jE+1; j++){
      mat->mat.item[j] = item[j-i-1] + 1;
    }
    monolis_qsort_int(&(mat->mat.item[jS]), 1, jE-jS);
  }

  //monolis_get_CRR_format();
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

  monolis_convert_connectivity_to_nodal(
    nnode,
    nelem,
    conn_index,
    con,
    &index,
    &item);

  monolis_get_nonzero_pattern_by_nodal(
    mat,
    nnode,
    ndof,
    index,
    item);

  free(conn_index);
  free(con);
  //free(index);
  //free(item);
}

void monolis_add_sparse_matrix(
  MONOLIS* mat,
  int      nbase_func,
  int      *connectivity,
  double** local_mat)
{
  monolis_add_sparse_matrix_c_main();
}

void monolis_set_Dirichlet_bc(
  MONOLIS* mat,
  double*  b,
  int      node_id,
  int      ndof_bc,
  double   val)
{
  monolis_set_Dirichlet_bc_c_main();
}

void monolis_solve(
  MONOLIS* mat,
  double*  b,
  double*  x)
{
  monolis_solve_c_main();
}
