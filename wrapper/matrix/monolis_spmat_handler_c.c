#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "monolis_utils.h"
#include "monolis_spmat_handler_util_c.h"
#include "monolis_spmat_nzpattern_util_c.h"
#include "monolis_def_struc_c.h"

void monolis_set_scalar_to_sparse_matrix_R(
  MONOLIS* mat,
  int      i,
  int      j,
  int      submat_i,
  int      submat_j,
  double   val)
{
  int n_node = mat->mat.NP;
  int n_dof = mat->mat.NDOF;
  int nz = mat->mat.CSR.index[n_node];

  monolis_set_scalar_to_sparse_matrix_R_c_main(
    n_node,
    nz,
    n_dof,
    mat->mat.CSR.index,
    mat->mat.CSR.item,
    mat->mat.R.A,
    i,
    j,
    submat_i,
    submat_j,
    val);
}

void monolis_add_scalar_to_sparse_matrix_R(
  MONOLIS* mat,
  int      i,
  int      j,
  int      submat_i,
  int      submat_j,
  double   val)
{
  int n_node = mat->mat.NP;
  int n_dof = mat->mat.NDOF;
  int nz = mat->mat.CSR.index[n_node];

  monolis_add_scalar_to_sparse_matrix_R_c_main(
    n_node,
    nz,
    n_dof,
    mat->mat.CSR.index,
    mat->mat.CSR.item,
    mat->mat.R.A,
    i,
    j,
    submat_i,
    submat_j,
    val);
}

void monolis_get_scalar_from_sparse_matrix_R(
  MONOLIS* mat,
  int      i,
  int      j,
  int      submat_i,
  int      submat_j,
  double*  val,
  bool*    is_find)
{
  int n = mat->mat.NP;
  int n_node = mat->mat.NP;
  int nz = mat->mat.CSR.index[n_node];
  int n_dof = mat->mat.NDOF;
  int is_find_t = 0;

  monolis_get_scalar_from_sparse_matrix_R_c_main(
    n,
    nz,
    n_dof,
    mat->mat.CSR.index,
    mat->mat.CSR.item,
    mat->mat.R.A,
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

void monolis_add_matrix_to_sparse_matrix_R(
  MONOLIS* mat,
  int      n_base,
  int*     connectivity,
  double** val)
{
  int n_node = mat->mat.NP;
  int n_dof = mat->mat.NDOF;
  int nz = mat->mat.CSR.index[n_node];

  monolis_add_matrix_to_sparse_matrix_main_R_c_main(
    n_node,
    nz,
    n_dof,
    mat->mat.CSR.index,
    mat->mat.CSR.item,
    mat->mat.R.A,
    n_base,
    n_base,
    connectivity,
    connectivity,
    val);
}

void monolis_add_matrix_to_sparse_matrix_offdiag_R(
  MONOLIS* mat,
  int      n_base1,
  int      n_base2,
  int*     connectivity1,
  int*     connectivity2,
  double** val)
{
  int n_node = mat->mat.NP;
  int n_dof = mat->mat.NDOF;
  int nz = mat->mat.CSR.index[n_node];

  monolis_add_matrix_to_sparse_matrix_main_R_c_main(
    n_node,
    nz,
    n_dof,
    mat->mat.CSR.index,
    mat->mat.CSR.item,
    mat->mat.R.A,
    n_base1,
    n_base2,
    connectivity1,
    connectivity2,
    val);
}

void monolis_set_matrix_BCSR_R(
  MONOLIS* mat,
  int      n,
  int      np,
  int      n_dof,
  int      nz,
  double*  A,
  int*     index,
  int*     item){

  mat->mat.N = n;
  mat->mat.NP = np;
  mat->mat.NDOF = n_dof;
  //mat->mat.A = (double*)calloc(n_dof*n_dof*nz, sizeof(double));
  //mat->mat.X = (double*)calloc(n_dof*np, sizeof(double));
  //mat->mat.B = (double*)calloc(n_dof*np, sizeof(double));
  //mat->mat.index = (int* )calloc(np+1, sizeof(int));
  //mat->mat.item = (int*)calloc(nz, sizeof(int));

  int i;
  for(i = 0; i < np + 1; i++) {
    mat->mat.CSR.index[i] = index[i];
  }

  for(i = 0; i < nz; i++) {
    mat->mat.CSR.item[i] = item[i] + 1;
  }

  for(i = 0; i < n_dof*n_dof*nz; i++) {
    mat->mat.R.A[i] = A[i];
  }

  //mat->mat.indexR = (int*)calloc(np+1, sizeof(int));
  //mat->mat.itemR = (int*)calloc(nz, sizeof(int));
  //mat->mat.permR = (int*)calloc(nz, sizeof(int));

  monolis_get_CSC_format(
    n,
    np,
    nz,
    mat->mat.CSR.index,
    mat->mat.CSR.item,
    mat->mat.CSC.index,
    mat->mat.CSC.item,
    mat->mat.CSC.perm);
}

void monolis_set_matrix_BCSR_mat_val_R(
  MONOLIS* mat,
  int      n_dof,
  int      nz,
  double*  A){

  int i;
  for(i = 0; i < n_dof*n_dof*nz; i++) {
    mat->mat.R.A[i] = A[i];
  }
}

void monolis_set_Dirichlet_bc_R(
  MONOLIS* mat,
  double*  b,
  int      node_id,
  int      n_dof_bc,
  double   val)
{
  int n_node = mat->mat.NP;
  int n_dof = mat->mat.NDOF;
  int nz = mat->mat.CSR.index[n_node];

  monolis_set_Dirichlet_bc_main_R_c(
    n_node,
    nz,
    n_dof,
    mat->mat.CSR.index,
    mat->mat.CSR.item,
    mat->mat.CSC.index,
    mat->mat.CSC.item,
    mat->mat.CSC.perm,
    mat->mat.R.A,
    b,
    node_id,
    n_dof_bc,
    val);
}
