#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "monolis_utils.h"
#include "monolis_def_struc_c.h"
#include "monolis_spmat_copy_c.h"
#include "monolis_vec_util_c.h"

void monolis_copy_mat_R(
  MONOLIS* mat_in,
  MONOLIS* mat_out)
{
  monolis_finalize(mat_out);

  monolis_copy_mat_nonzero_pattern_R(mat_in, mat_out);

  monolis_copy_mat_value_R(mat_in, mat_out);
}

void monolis_copy_mat_C(
  MONOLIS* mat_in,
  MONOLIS* mat_out)
{
  monolis_finalize(mat_out);

  monolis_copy_mat_nonzero_pattern_C(mat_in, mat_out);
  monolis_copy_mat_value_C(mat_in, mat_out);
}

void monolis_copy_mat_nonzero_pattern_R(
  MONOLIS* mat_in,
  MONOLIS* mat_out)
{
  int NP, NZ, NZU, NZL;

  mat_out->mat.N = mat_in->mat.N;
  mat_out->mat.NP = mat_in->mat.NP;
  mat_out->mat.NDOF = mat_in->mat.NDOF;

  monolis_copy_mat_nonzero_pattern_CSR (mat_in->mat.NP, &mat_in->mat.CSR,  &mat_out->mat.CSR);
  monolis_copy_mat_nonzero_pattern_CSC (mat_in->mat.NP, &mat_in->mat.CSC,  &mat_out->mat.CSC);
  monolis_copy_mat_nonzero_pattern_SCSR(mat_in->mat.NP, &mat_in->mat.SCSR, &mat_out->mat.SCSR);

  NP = mat_in->mat.NP;
  NZ = mat_in->mat.CSR.index[NP];

  NZU = 0;
  if(mat_in->mat.SCSR.indexU != NULL){
    NZU = mat_in->mat.SCSR.indexU[NP];
  }

  NZL = 0;
  if(mat_in->mat.SCSR.indexL != NULL){
    NZL = mat_in->mat.SCSR.indexL[NP];
  }

  monolis_copy_mat_nonzero_pattern_val_R(
    mat_in->mat.NP,
    mat_in->mat.NDOF,
    NZ, NZU, NZL,
    &mat_in->mat.R,
    &mat_out->mat.R);
}

void monolis_copy_mat_nonzero_pattern_C(
  MONOLIS* mat_in,
  MONOLIS* mat_out)
{
  int NP, NZ, NZU, NZL;

  mat_out->mat.N = mat_in->mat.N;
  mat_out->mat.NP = mat_in->mat.NP;
  mat_out->mat.NDOF = mat_in->mat.NDOF;

  monolis_copy_mat_nonzero_pattern_CSR (mat_in->mat.NP, &mat_in->mat.CSR,  &mat_out->mat.CSR);
  monolis_copy_mat_nonzero_pattern_CSC (mat_in->mat.NP, &mat_in->mat.CSC,  &mat_out->mat.CSC);
  monolis_copy_mat_nonzero_pattern_SCSR(mat_in->mat.NP, &mat_in->mat.SCSR, &mat_out->mat.SCSR);

  NP = mat_in->mat.NP;
  NZ = mat_in->mat.CSR.index[NP];

  NZU = 0;
  if(mat_in->mat.SCSR.indexU != NULL){
    NZU = mat_in->mat.SCSR.indexU[NP];
  }

  NZL = 0;
  if(mat_in->mat.SCSR.indexL != NULL){
    NZL = mat_in->mat.SCSR.indexL[NP];
  }

  monolis_copy_mat_nonzero_pattern_val_C(
    mat_in->mat.NP,
    mat_in->mat.NDOF,
    NZ, NZU, NZL,
    &mat_in->mat.C,
    &mat_out->mat.C);
}

void monolis_copy_mat_nonzero_pattern_val_R(
  int      NP,
  int      NDOF,
  int      NZ,
  int      NZU,
  int      NZL,
  MONOLIS_MAT_VAL_R* mat_in,
  MONOLIS_MAT_VAL_R* mat_out)
{
  monolis_dealloc_R_1d(&mat_out->X);
  mat_out->X = monolis_alloc_R_1d(mat_out->X, NP*NDOF);
  monolis_vec_copy_R(NP, NDOF, mat_in->X, mat_out->X);

  monolis_dealloc_R_1d(&mat_out->B);
  mat_out->B = monolis_alloc_R_1d(mat_out->B, NP*NDOF);
  monolis_vec_copy_R(NP, NDOF, mat_in->B, mat_out->B);

  monolis_dealloc_R_1d(&mat_out->A);
  mat_out->A = monolis_alloc_R_1d(mat_out->A, NZ*NDOF*NDOF);
  monolis_vec_copy_R(NZ, NDOF*NDOF, mat_in->A, mat_out->A);

  if(mat_in->U != NULL){
    monolis_dealloc_R_1d(&mat_out->U);
    mat_out->U = monolis_alloc_R_1d(mat_out->U, NZU*NDOF*NDOF);
    monolis_vec_copy_R(NZU, NDOF*NDOF, mat_in->U, mat_out->U);
  }

  if(mat_in->L != NULL){
    monolis_dealloc_R_1d(&mat_out->L);
    mat_out->L = monolis_alloc_R_1d(mat_out->L, NZ*NDOF*NDOF);
    monolis_vec_copy_R(NZL, NDOF*NDOF, mat_in->L, mat_out->L);
  }

  if(mat_in->D != NULL){
    monolis_dealloc_R_1d(&mat_out->D);
    mat_out->D = monolis_alloc_R_1d(mat_out->D, NP*NDOF*NDOF);
    monolis_vec_copy_R(NP, NDOF*NDOF, mat_in->D, mat_out->D);
  }
}

void monolis_copy_mat_nonzero_pattern_val_C(
  int      NP,
  int      NDOF,
  int      NZ,
  int      NZU,
  int      NZL,
  MONOLIS_MAT_VAL_C* mat_in,
  MONOLIS_MAT_VAL_C* mat_out)
{
  monolis_dealloc_C_1d(&mat_out->X);
  mat_out->X = monolis_alloc_C_1d(mat_out->X, NP*NDOF);
  monolis_vec_copy_C(NP, NDOF, mat_in->X, mat_out->X);

  monolis_dealloc_C_1d(&mat_out->B);
  mat_out->B = monolis_alloc_C_1d(mat_out->B, NP*NDOF);
  monolis_vec_copy_C(NP, NDOF, mat_in->B, mat_out->B);

  monolis_dealloc_C_1d(&mat_out->A);
  mat_out->A = monolis_alloc_C_1d(mat_out->A, NZ*NDOF*NDOF);
  monolis_vec_copy_C(NZ, NDOF*NDOF, mat_in->A, mat_out->A);

  if(mat_in->U != NULL){
    monolis_dealloc_C_1d(&mat_out->U);
    mat_out->U = monolis_alloc_C_1d(mat_out->U, NZU*NDOF*NDOF);
    monolis_vec_copy_C(NZU, NDOF*NDOF, mat_in->U, mat_out->U);
  }

  if(mat_in->L != NULL){
    monolis_dealloc_C_1d(&mat_out->L);
    mat_out->L = monolis_alloc_C_1d(mat_out->L, NZ*NDOF*NDOF);
    monolis_vec_copy_C(NZL, NDOF*NDOF, mat_in->L, mat_out->L);
  }

  if(mat_in->D != NULL){
    monolis_dealloc_C_1d(&mat_out->D);
    mat_out->D = monolis_alloc_C_1d(mat_out->D, NP*NDOF*NDOF);
    monolis_vec_copy_C(NP, NDOF*NDOF, mat_in->D, mat_out->D);
  }
}

void monolis_copy_mat_nonzero_pattern_CSR(
  int      NP,
  MONOLIS_MAT_CSR* mat_in,
  MONOLIS_MAT_CSR* mat_out)
{
  int NZ;
  NZ = mat_in->index[NP];

  monolis_dealloc_I_1d(&mat_out->index);
  mat_out->index = monolis_alloc_I_1d(mat_out->index, NP + 1);
  monolis_vec_copy_I(NP + 1, 1, mat_in->index, mat_out->index);

  monolis_dealloc_I_1d(&mat_out->item);
  mat_out->item = monolis_alloc_I_1d(mat_out->item, NZ);
  monolis_vec_copy_I(NZ, 1, mat_in->item, mat_out->item);
}

void monolis_copy_mat_nonzero_pattern_CSC(
  int      NP,
  MONOLIS_MAT_CSC* mat_in,
  MONOLIS_MAT_CSC* mat_out)
{
  int NZ;
  NZ = mat_in->index[NP];

  monolis_dealloc_I_1d(&mat_out->index);
  mat_out->index = monolis_alloc_I_1d(mat_out->index, NP + 1);
  monolis_vec_copy_I(NP + 1, 1, mat_in->index, mat_out->index);

  monolis_dealloc_I_1d(&mat_out->item);
  mat_out->item = monolis_alloc_I_1d(mat_out->item, NZ);
  monolis_vec_copy_I(NZ, 1, mat_in->item, mat_out->item);

  monolis_dealloc_I_1d(&mat_out->perm);
  mat_out->perm = monolis_alloc_I_1d(mat_out->perm, NZ);
  monolis_vec_copy_I(NZ, 1, mat_in->perm, mat_out->perm);
}

void monolis_copy_mat_nonzero_pattern_SCSR(
  int      NP,
  MONOLIS_MAT_SEPARATED_CSR* mat_in,
  MONOLIS_MAT_SEPARATED_CSR* mat_out)
{
  int NZ;

  if(mat_in->indexU != NULL){
    monolis_dealloc_I_1d(&mat_out->indexU);
    mat_out->indexU = monolis_alloc_I_1d(mat_out->indexU, NP + 1);
    monolis_vec_copy_I(NP + 1, 1, mat_in->indexU, mat_out->indexU);
  }

  if(mat_in->itemU != NULL){
    NZ = mat_in->indexU[NP];
    monolis_dealloc_I_1d(&mat_out->itemU);
    mat_out->itemU = monolis_alloc_I_1d(mat_out->itemU, NZ);
    monolis_vec_copy_I(NZ, 1, mat_in->itemU, mat_out->itemU);
  }

  if(mat_in->indexL != NULL){
    monolis_dealloc_I_1d(&mat_out->indexL);
    mat_out->indexL = monolis_alloc_I_1d(mat_out->indexL, NP + 1);
    monolis_vec_copy_I(NP + 1, 1, mat_in->indexL, mat_out->indexL);
  }

  if(mat_in->itemL != NULL){
    NZ = mat_in->indexL[NP];
    monolis_dealloc_I_1d(&mat_out->itemL);
    mat_out->itemL = monolis_alloc_I_1d(mat_out->itemL, NZ);
    monolis_vec_copy_I(NZ, 1, mat_in->itemL, mat_out->itemL);
  }
}

void monolis_copy_mat_value_R(
  MONOLIS* mat_in,
  MONOLIS* mat_out)
{
  monolis_copy_mat_value_matrix_R(mat_in, mat_out);
  monolis_copy_mat_value_solution_R(mat_in, mat_out);
  monolis_copy_mat_value_rhs_R(mat_in, mat_out);
}

void monolis_copy_mat_value_matrix_R(
  MONOLIS* mat_in,
  MONOLIS* mat_out)
{
  int i;
  int n_dof = mat_in->mat.NDOF;
  int nz = mat_in->mat.CSR.index[mat_in->mat.NP];

  for(i = 0; i < n_dof*n_dof*nz; i++) {
    mat_out->mat.R.A[i] = mat_in->mat.R.A[i];
  }
}

void monolis_copy_mat_value_solution_R(
  MONOLIS* mat_in,
  MONOLIS* mat_out)
{
  int i;
  int n_dof = mat_in->mat.NDOF;
  int np = mat_in->mat.NP;

  for(i = 0; i < n_dof*np; i++) {
    mat_out->mat.R.X[i] = mat_in->mat.R.X[i];
  }
}

void monolis_copy_mat_value_rhs_R(
  MONOLIS* mat_in,
  MONOLIS* mat_out)
{
  int i;
  int n_dof = mat_in->mat.NDOF;
  int np = mat_in->mat.NP;

  for(i = 0; i < n_dof*np; i++) {
    mat_out->mat.R.B[i] = mat_in->mat.R.B[i];
  }
}

void monolis_copy_mat_value_C(
  MONOLIS* mat_in,
  MONOLIS* mat_out)
{
  monolis_copy_mat_value_matrix_C(mat_in, mat_out);
  monolis_copy_mat_value_solution_C(mat_in, mat_out);
  monolis_copy_mat_value_rhs_C(mat_in, mat_out);
}

void monolis_copy_mat_value_matrix_C(
  MONOLIS* mat_in,
  MONOLIS* mat_out)
{
  int i;
  int n_dof = mat_in->mat.NDOF;
  int nz = mat_in->mat.CSR.index[mat_in->mat.N];

  for(i = 0; i < n_dof*n_dof*nz; i++) {
    mat_out->mat.C.A[i] = mat_in->mat.C.A[i];
  }
}

void monolis_copy_mat_value_solution_C(
  MONOLIS* mat_in,
  MONOLIS* mat_out)
{
  int i;
  int n_dof = mat_in->mat.NDOF;
  int np = mat_in->mat.NP;

  for(i = 0; i < n_dof*np; i++) {
    mat_out->mat.C.X[i] = mat_in->mat.C.X[i];
  }
}

void monolis_copy_mat_value_rhs_C(
  MONOLIS* mat_in,
  MONOLIS* mat_out)
{
  int i;
  int n_dof = mat_in->mat.NDOF;
  int np = mat_in->mat.NP;

  for(i = 0; i < n_dof*np; i++) {
    mat_out->mat.C.B[i] = mat_in->mat.C.B[i];
  }
}

void monolis_clear_mat_value_R(
  MONOLIS* mat)
{
  monolis_clear_mat_value_matrix_R(mat);
  monolis_clear_mat_value_solution_R(mat);
  monolis_clear_mat_value_rhs_R(mat);
}

void monolis_clear_mat_value_matrix_R(
  MONOLIS* mat)
{
  int i;
  int n_dof = mat->mat.NDOF;
  int nz = mat->mat.CSR.index[mat->mat.N];

  for(i = 0; i < n_dof*n_dof*nz; i++) {
    mat->mat.R.A[i] = 0.0;
  }
}

void monolis_clear_mat_value_solution_R(
  MONOLIS* mat)
{
  int i;
  int n_dof = mat->mat.NDOF;
  int np = mat->mat.NP;

  for(i = 0; i < n_dof*np; i++) {
    mat->mat.R.X[i] = 0.0;
  }
}

void monolis_clear_mat_value_rhs_R(
  MONOLIS* mat)
{
  int i;
  int n_dof = mat->mat.NDOF;
  int np = mat->mat.NP;

  for(i = 0; i < n_dof*np; i++) {
    mat->mat.R.B[i] = 0.0;
  }
}

void monolis_clear_mat_value_C(
  MONOLIS* mat)
{
  monolis_clear_mat_value_matrix_C(mat);
  monolis_clear_mat_value_solution_C(mat);
  monolis_clear_mat_value_rhs_C(mat);
}

void monolis_clear_mat_value_matrix_C(
  MONOLIS* mat)
{
  int i;
  int n_dof = mat->mat.NDOF;
  int nz = mat->mat.CSR.index[mat->mat.N];

  for(i = 0; i < n_dof*n_dof*nz; i++) {
    mat->mat.C.A[i] = 0.0;
  }
}

void monolis_clear_mat_value_solution_C(
  MONOLIS* mat)
{
  int i;
  int n_dof = mat->mat.NDOF;
  int np = mat->mat.NP;

  for(i = 0; i < n_dof*np; i++) {
    mat->mat.C.X[i] = 0.0;
  }
}

void monolis_clear_mat_value_rhs_C(
  MONOLIS* mat)
{
  int i;
  int n_dof = mat->mat.NDOF;
  int np = mat->mat.NP;

  for(i = 0; i < n_dof*np; i++) {
    mat->mat.C.B[i] = 0.0;
  }
}
