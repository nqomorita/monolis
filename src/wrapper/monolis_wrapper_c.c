#include <stdio.h>
#include <stdlib.h>
#include "monolis.h"

/* header */

void monolis_initialize_c_main(
    // prm
    int method,
    int precond,
    int curiter,
    int maxiter,
    int ierr,
    double tol,
    double curresid,
    // mat
    int N,
    int NP,
    int NZ,
    int NDOF,
    // com
    int myrank,
    int comm,
    int commsize,
    int recv_n_neib,
    int send_n_neib);

void monolis_finalize_c_main();
void monolis_get_nonzero_pattern_c_main(
    int nnode,
    int nbase_func,
    int ndof,
    int nelem,
    int* elem_t,
    int* index,
    int* item,
    int* A,
    double* B,
    double* X);

void monolis_add_sparse_matrix_c_main();
void monolis_set_Dirichlet_bc_c_main();
void monolis_solve_c_main();

/* body */

void monolis_initialize(
  MONOLIS* mat)
{
  monolis_initialize_c_main(
    // prm
    mat->prm.method,
    mat->prm.precond,
    mat->prm.curiter,
    mat->prm.maxiter,
    mat->prm.ierr,
    mat->prm.tol,
    mat->prm.curresid,
    //mat->prm.is_scaling,
    //mat->prm.is_reordering,
    //mat->prm.is_init_x,
    //mat->prm.is_debug,
    //mat->prm.show_iterlog,
    //mat->prm.show_time,
    //mat->prm.show_summary
    // mat
    mat->mat.N,
    mat->mat.NP,
    mat->mat.NZ,
    mat->mat.NDOF,
    // comm
    mat->com.myrank,
    mat->com.comm,
    mat->com.commsize,
    mat->com.recv_n_neib,
    mat->com.send_n_neib
    );
}

void monolis_finalize(
  MONOLIS* mat)
{
  monolis_finalize_c_main();
}

void monolis_get_nonzero_pattern(
  MONOLIS* mat,
  int      nnode,
  int      nbase_func,
  int      ndof,
  int      nelem,
  int**    elem)
{
  int* elem_t = (int*)calloc(nbase_func*nelem, sizeof(int));

  for(int i=0; i<nelem; i++) {
    for(int j=0; j<nbase_func; j++) {
      elem_t[nelem*i+j] = elem[i][j];
    }
  }

  monolis_get_nonzero_pattern_c_main(
    nnode,
    nbase_func,
    ndof,
    nelem,
    elem_t,
    mat->mat.index,
    mat->mat.item,
    mat->mat.A,
    mat->mat.B,
    mat->mat.X
    );

  free(elem_t);
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
