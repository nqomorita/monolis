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
      jE - jS,
      0,
      jE - jS - 1);
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

void monolis_alloc_nonzero_pattern_mat_val_V_R(
  MONOLIS_MAT* mat)
{
  int n, nz, in;

  n = mat->n_dof_index[mat->NP];
  in = mat->CSR.index[mat->NP];
  nz = mat->n_dof_index2[in];

  mat->R.A = monolis_alloc_R_1d(mat->R.A, nz);
  mat->R.X = monolis_alloc_R_1d(mat->R.X, n);
  mat->R.B = monolis_alloc_R_1d(mat->R.B, n);

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

void monolis_alloc_nonzero_pattern_mat_val_V_C(
  MONOLIS_MAT* mat)
{
  int n, nz, in;

  n = mat->n_dof_index[mat->NP];
  in = mat->CSR.index[mat->NP];
  nz = mat->n_dof_index2[in];

  mat->C.A = monolis_alloc_C_1d(mat->C.A, nz);
  mat->C.X = monolis_alloc_C_1d(mat->C.X, n);
  mat->C.B = monolis_alloc_C_1d(mat->C.B, n);

  mat->C.L = NULL;
  mat->C.U = NULL;
  mat->C.D = NULL;
}

void monolis_set_n_dof_index(
  MONOLIS_MAT* mat,
  int*         n_dof_list)
{
  int np, nz;

  np = mat->NP;
  nz = mat->CSR.index[np];

  mat->n_dof_list = monolis_alloc_I_1d(mat->n_dof_list, np);
  mat->n_dof_index = monolis_alloc_I_1d(mat->n_dof_index, np + 1);
  mat->n_dof_index2 = monolis_alloc_I_1d(mat->n_dof_index2, nz + 1);

  for (int i = 0; i < np; ++i) {
    mat->n_dof_list[i] = n_dof_list[i];
  }

  for (int i = 0; i < np; ++i) {
    mat->n_dof_index[i + 1] = mat->n_dof_index[i] + n_dof_list[i];
  }

  for (int i = 0; i < np; ++i) {
    int jS = mat->CSR.index[i];
    int jE = mat->CSR.index[i + 1];
    for (int j = jS; j < jE; ++j) {
      int in = mat->CSR.item[j] - 1;
      mat->n_dof_index2[j + 1] = mat->n_dof_index2[j] + n_dof_list[i]*n_dof_list[in];
    }
  }
}
