#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "monolis_utils.h"
#include "monolis_def_struc_c.h"
#include "monolis_spmat_copy_c.h"

void monolis_copy_mat_value_R(
  MONOLIS* mat_in,
  MONOLIS* mat_out)
{
  monolis_copy_mat_value_A_R(mat_in, mat_out);
  monolis_copy_mat_value_X_R(mat_in, mat_out);
  monolis_copy_mat_value_B_R(mat_in, mat_out);
}

void monolis_copy_mat_value_A_R(
  MONOLIS* mat_in,
  MONOLIS* mat_out)
{
  int i;
  int n_dof = mat_in->mat.NDOF;
  int nz = mat_in->mat.CSR.index[mat_in->mat.N];

  for(i = 0; i < n_dof*n_dof*nz; i++) {
    mat_out->mat.R.A[i] = mat_in->mat.R.A[i];
  }
}

void monolis_copy_mat_value_X_R(
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

void monolis_copy_mat_value_B_R(
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
  monolis_copy_mat_value_A_C(mat_in, mat_out);
  monolis_copy_mat_value_X_C(mat_in, mat_out);
  monolis_copy_mat_value_B_C(mat_in, mat_out);
}

void monolis_copy_mat_value_A_C(
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

void monolis_copy_mat_value_X_C(
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

void monolis_copy_mat_value_B_C(
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
  monolis_clear_mat_value_A_R(mat);
  monolis_clear_mat_value_X_R(mat);
  monolis_clear_mat_value_B_R(mat);
}

void monolis_clear_mat_value_A_R(
  MONOLIS* mat)
{
  int i;
  int n_dof = mat->mat.NDOF;
  int nz = mat->mat.CSR.index[mat->mat.N];

  for(i = 0; i < n_dof*n_dof*nz; i++) {
    mat->mat.R.A[i] = 0.0;
  }
}

void monolis_clear_mat_value_X_R(
  MONOLIS* mat)
{
  int i;
  int n_dof = mat->mat.NDOF;
  int np = mat->mat.NP;

  for(i = 0; i < n_dof*np; i++) {
    mat->mat.R.X[i] = 0.0;
  }
}

void monolis_clear_mat_value_B_R(
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
  monolis_clear_mat_value_A_C(mat);
  monolis_clear_mat_value_X_C(mat);
  monolis_clear_mat_value_B_C(mat);
}

void monolis_clear_mat_value_A_C(
  MONOLIS* mat)
{
  int i;
  int n_dof = mat->mat.NDOF;
  int nz = mat->mat.CSR.index[mat->mat.N];

  for(i = 0; i < n_dof*n_dof*nz; i++) {
    mat->mat.C.A[i] = 0.0;
  }
}

void monolis_clear_mat_value_X_C(
  MONOLIS* mat)
{
  int i;
  int n_dof = mat->mat.NDOF;
  int np = mat->mat.NP;

  for(i = 0; i < n_dof*np; i++) {
    mat->mat.C.X[i] = 0.0;
  }
}

void monolis_clear_mat_value_B_C(
  MONOLIS* mat)
{
  int i;
  int n_dof = mat->mat.NDOF;
  int np = mat->mat.NP;

  for(i = 0; i < n_dof*np; i++) {
    mat->mat.C.B[i] = 0.0;
  }
}
