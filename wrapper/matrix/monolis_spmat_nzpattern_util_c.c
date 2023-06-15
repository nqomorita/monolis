#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "monolis_utils.h"
#include "monolis_def_solver_prm_c.h"
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

  mat->CSR.index = monolis_alloc_I_1d(mat->CSR.index, n_node + 1);

  for(i = 1; i < n_node + 1; i++) {
    mat->CSR.index[i] = index[i] + i;
  }

  nz = mat->CSR.index[n_node];

  mat->CSR.item = monolis_alloc_I_1d(mat->CSR.item, nz);

  for(i = 0; i < n_node; i++) {
    jS = mat->CSR.index[i];
    jE = mat->CSR.index[i + 1];
    mat->CSR.item[jS] = i + 1;
    for(j = jS + 1; j < jE; j++){
      mat->CSR.item[j] = item[j-i-1] + 1;
    }
    monolis_qsort_I_1d(
      &(mat->CSR.item[jS]),
      1,
      jE - jS);
  }

  mat->CSC.index = monolis_alloc_I_1d(mat->CSC.index, n_node + 1);
  mat->CSC.item = monolis_alloc_I_1d(mat->CSC.item, nz);
  mat->CSC.perm = monolis_alloc_I_1d(mat->CSC.perm, nz);

  monolis_get_CSC_format(
    n_node,
    n_node,
    nz,
    mat->CSR.index,
    mat->CSR.item,
    mat->CSC.index,
    mat->CSC.item,
    mat->CSC.perm);

  mat->SCSR.indexU = NULL;
  mat->SCSR.itemU = NULL;
  mat->SCSR.indexL = NULL;
  mat->SCSR.itemL = NULL;
}

void monolis_alloc_nonzero_pattern_mat_val_R(
  MONOLIS_MAT* mat)
{
  int n_node, nz, n_dof;

  n_node = mat->N;
  n_dof = mat->NDOF;
  nz = mat->CSR.index[n_node];

  mat->R.A = monolis_alloc_R_1d(mat->R.A, n_dof*n_dof*nz);
  mat->R.X = monolis_alloc_R_1d(mat->R.X, n_dof*n_node);
  mat->R.B = monolis_alloc_R_1d(mat->R.B, n_dof*n_node);

  mat->R.L = NULL;
  mat->R.U = NULL;
  mat->R.D = NULL;
}

void monolis_alloc_nonzero_pattern_mat_val_C(
  MONOLIS_MAT* mat)
{
  int n_node, nz, n_dof;

  n_node = mat->N;
  n_dof = mat->NDOF;
  nz = mat->CSR.index[n_node];

  mat->C.A = monolis_alloc_C_1d(mat->C.A, n_dof*n_dof*nz);
  mat->C.X = monolis_alloc_C_1d(mat->C.X, n_dof*n_node);
  mat->C.B = monolis_alloc_C_1d(mat->C.B, n_dof*n_node);

  mat->C.L = NULL;
  mat->C.U = NULL;
  mat->C.D = NULL;
}
