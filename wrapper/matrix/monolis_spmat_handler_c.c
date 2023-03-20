#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <complex.h>
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
  double* val_t;

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
    val_t);
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
  double* val_t;

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
    val_t);
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

  monolis_alloc_I_1d(mat->mat.CSR.index, np + 1);
  monolis_alloc_I_1d(mat->mat.CSR.item, nz);

  monolis_alloc_R_1d(mat->mat.R.A, n_dof*n_dof*nz);
  monolis_alloc_R_1d(mat->mat.R.X, n_dof*np);
  monolis_alloc_R_1d(mat->mat.R.B, n_dof*np);

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

  monolis_alloc_I_1d(mat->mat.CSC.index, np + 1);
  monolis_alloc_I_1d(mat->mat.CSC.item, nz);
  monolis_alloc_I_1d(mat->mat.CSC.perm, nz);

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

  monolis_set_Dirichlet_bc_R_c_main(
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

void monolis_set_scalar_to_sparse_matrix_C(
  MONOLIS*       mat,
  int            i,
  int            j,
  int            submat_i,
  int            submat_j,
  double complex val)
{
  int n_node = mat->mat.NP;
  int n_dof = mat->mat.NDOF;
  int nz = mat->mat.CSR.index[n_node];

  monolis_set_scalar_to_sparse_matrix_C_c_main(
    n_node,
    nz,
    n_dof,
    mat->mat.CSR.index,
    mat->mat.CSR.item,
    mat->mat.C.A,
    i,
    j,
    submat_i,
    submat_j,
    val);
}

void monolis_add_scalar_to_sparse_matrix_C(
  MONOLIS*       mat,
  int            i,
  int            j,
  int            submat_i,
  int            submat_j,
  double complex val)
{
  int n_node = mat->mat.NP;
  int n_dof = mat->mat.NDOF;
  int nz = mat->mat.CSR.index[n_node];

  monolis_add_scalar_to_sparse_matrix_C_c_main(
    n_node,
    nz,
    n_dof,
    mat->mat.CSR.index,
    mat->mat.CSR.item,
    mat->mat.C.A,
    i,
    j,
    submat_i,
    submat_j,
    val);
}

void monolis_get_scalar_from_sparse_matrix_C(
  MONOLIS*        mat,
  int             i,
  int             j,
  int             submat_i,
  int             submat_j,
  double complex* val,
  bool*           is_find)
{
  int n = mat->mat.NP;
  int n_node = mat->mat.NP;
  int nz = mat->mat.CSR.index[n_node];
  int n_dof = mat->mat.NDOF;
  int is_find_t = 0;

  monolis_get_scalar_from_sparse_matrix_C_c_main(
    n,
    nz,
    n_dof,
    mat->mat.CSR.index,
    mat->mat.CSR.item,
    mat->mat.C.A,
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

void monolis_add_matrix_to_sparse_matrix_C(
  MONOLIS*         mat,
  int              n_base,
  int*             connectivity,
  double complex** val)
{
  int n_node = mat->mat.NP;
  int n_dof = mat->mat.NDOF;
  int nz = mat->mat.CSR.index[n_node];

  monolis_add_matrix_to_sparse_matrix_main_C_c_main(
    n_node,
    nz,
    n_dof,
    mat->mat.CSR.index,
    mat->mat.CSR.item,
    mat->mat.C.A,
    n_base,
    n_base,
    connectivity,
    connectivity,
    val);
}

void monolis_add_matrix_to_sparse_matrix_offdiag_C(
  MONOLIS*         mat,
  int              n_base1,
  int              n_base2,
  int*             connectivity1,
  int*             connectivity2,
  double complex** val)
{
  int n_node = mat->mat.NP;
  int n_dof = mat->mat.NDOF;
  int nz = mat->mat.CSR.index[n_node];

  monolis_add_matrix_to_sparse_matrix_main_C_c_main(
    n_node,
    nz,
    n_dof,
    mat->mat.CSR.index,
    mat->mat.CSR.item,
    mat->mat.C.A,
    n_base1,
    n_base2,
    connectivity1,
    connectivity2,
    val);
}

void monolis_set_matrix_BCSR_C(
  MONOLIS*        mat,
  int             n,
  int             np,
  int             n_dof,
  int             nz,
  double complex* A,
  int*            index,
  int*            item){

  mat->mat.N = n;
  mat->mat.NP = np;
  mat->mat.NDOF = n_dof;

  monolis_alloc_I_1d(mat->mat.CSR.index, np + 1);
  monolis_alloc_I_1d(mat->mat.CSR.item, nz);

  monolis_alloc_C_1d(mat->mat.C.A, n_dof*n_dof*nz);
  monolis_alloc_C_1d(mat->mat.C.X, n_dof*np);
  monolis_alloc_C_1d(mat->mat.C.B, n_dof*np);

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

  monolis_alloc_I_1d(mat->mat.CSC.index, np + 1);
  monolis_alloc_I_1d(mat->mat.CSC.item, nz);
  monolis_alloc_I_1d(mat->mat.CSC.perm, nz);

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

void monolis_set_matrix_BCSR_mat_val_C(
  MONOLIS*        mat,
  int             n_dof,
  int             nz,
  double complex* A){

  int i;
  for(i = 0; i < n_dof*n_dof*nz; i++) {
    mat->mat.C.A[i] = A[i];
  }
}

void monolis_set_Dirichlet_bc_C(
  MONOLIS*        mat,
  double complex* b,
  int             node_id,
  int             n_dof_bc,
  double complex  val)
{
  int n_node = mat->mat.NP;
  int n_dof = mat->mat.NDOF;
  int nz = mat->mat.CSR.index[n_node];

  monolis_set_Dirichlet_bc_C_c_main(
    n_node,
    nz,
    n_dof,
    mat->mat.CSR.index,
    mat->mat.CSR.item,
    mat->mat.CSC.index,
    mat->mat.CSC.item,
    mat->mat.CSC.perm,
    mat->mat.C.A,
    b,
    node_id,
    n_dof_bc,
    val);
}