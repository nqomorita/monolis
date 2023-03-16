#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "monolis_utils.h"
#include "monolis_def_solver_c.h"
#include "monolis_def_struc_c.h"
#include "monolis_spmat_nzpattern_util_c.h"

void monolis_get_nonzero_pattern_by_nodal_graph_main(
  MONOLIS_MAT* mat,
  int          n_node,
  int          n_dof,
  int*         index,
  int*         item)
{
  int i, j;
  int nz, jS, jE;

  mat->N = n_node;
  mat->NP = n_node;
  mat->NDOF = n_dof;

  monolis_alloc_I_1d(mat->CSR.index, n_node + 1);

  for(i = 1; i < n_node + 1; i++) {
    mat->CSR.index[i] = index[i] + i;
  }

  nz = mat->CSR.index[n_node];

  monolis_alloc_I_1d(mat->CSR.item, nz);

  for(i = 0; i < n_node; i++) {
    jS = mat->CSR.index[i];
    jE = mat->CSR.index[i + 1];
    mat->CSR.item[jS] = i + 1;
    for(j = jS + 1; j < jE; j++){
      mat->CSR.item[j] = item[j-i-1] + 1;
    }
    //monolis_qsort_int(
    //  &(mat->CSR.item[jS]),
    //  1,
    //  jE - jS);
  }

  monolis_alloc_I_1d(mat->CSC.index, n_node + 1);
  monolis_alloc_I_1d(mat->CSC.item, nz);
  monolis_alloc_I_1d(mat->CSC.perm, nz);

  monolis_get_CSC_format(
    n_node,
    n_node,
    nz,
    mat->CSR.index,
    mat->CSR.item,
    mat->CSC.index,
    mat->CSC.item,
    mat->CSC.perm);
}

void monolis_alloc_nonzero_pattern_mat_val_R(
  MONOLIS_MAT* mat)
{
  int n_node, nz, n_dof;

  n_node = mat->N;
  n_dof = mat->NDOF;
  nz = mat->CSR.index[n_node];

  monolis_alloc_R_1d(mat->R.A, n_dof*n_dof*nz);
  monolis_alloc_R_1d(mat->R.X, n_dof*n_node);
  monolis_alloc_R_1d(mat->R.B, n_dof*n_node);
}

void monolis_alloc_nonzero_pattern_mat_val_C(
  MONOLIS_MAT* mat)
{
  int n_node, nz, n_dof;

  n_node = mat->N;
  n_dof = mat->NDOF;
  nz = mat->CSR.index[n_node];

  monolis_alloc_C_1d(mat->C.A, n_dof*n_dof*nz);
  monolis_alloc_C_1d(mat->C.X, n_dof*n_node);
  monolis_alloc_C_1d(mat->C.B, n_dof*n_node);
}
