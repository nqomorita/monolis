#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "monolis_utils.h"
#include "monolis_def_mat_c.h"

void monolis_mat_initialize_val_R(
  MONOLIS_MAT_VAL_R* mat)
{
  mat->A = NULL;
  mat->U = NULL;
  mat->D = NULL;
  mat->L = NULL;
  mat->X = NULL;
  mat->B = NULL;
}

void monolis_mat_initialize_val_C(
  MONOLIS_MAT_VAL_C* mat)
{
  mat->A = NULL;
  mat->U = NULL;
  mat->D = NULL;
  mat->L = NULL;
  mat->X = NULL;
  mat->B = NULL;
}

void monolis_mat_initialize_SCSR(
  MONOLIS_MAT_SEPARATED_CSR* mat)
{
  mat->indexU = NULL;
  mat->indexL = NULL;
  mat->itemU = NULL;
  mat->itemL = NULL;
}

void monolis_mat_initialize_CSR(
  MONOLIS_MAT_CSR* mat)
{
  mat->index = NULL;
  mat->item = NULL;
}

void monolis_mat_initialize_CSC(
  MONOLIS_MAT_CSC* mat)
{
  mat->index = NULL;
  mat->item = NULL;
  mat->perm = NULL;
}

void monolis_mat_initialize(
  MONOLIS_MAT* mat)
{
  mat->N = 0;
  mat->NP = 0;
  mat->NDOF = 0;
  mat->n_dof_list = NULL;
  mat->n_dof_index = NULL;
  mat->n_dof_index2 = NULL;
  monolis_mat_initialize_val_R(&mat->R);
  monolis_mat_initialize_val_C(&mat->C);
  monolis_mat_initialize_SCSR(&mat->SCSR);
  monolis_mat_initialize_CSR(&mat->CSR);
  monolis_mat_initialize_CSC(&mat->CSC);
}

void monolis_mat_finalize_val_R(
  MONOLIS_MAT_VAL_R* mat)
{
  monolis_dealloc_R_1d(&mat->A);
  monolis_dealloc_R_1d(&mat->U);
  monolis_dealloc_R_1d(&mat->D);
  monolis_dealloc_R_1d(&mat->L);
  monolis_dealloc_R_1d(&mat->X);
  monolis_dealloc_R_1d(&mat->B);
}

void monolis_mat_finalize_val_C(
  MONOLIS_MAT_VAL_C* mat)
{
  monolis_dealloc_C_1d(&mat->A);
  monolis_dealloc_C_1d(&mat->U);
  monolis_dealloc_C_1d(&mat->D);
  monolis_dealloc_C_1d(&mat->L);
  monolis_dealloc_C_1d(&mat->X);
  monolis_dealloc_C_1d(&mat->B);
}

void monolis_mat_finalize_SCSR(
  MONOLIS_MAT_SEPARATED_CSR* mat)
{
  monolis_dealloc_I_1d(&mat->indexU);
  monolis_dealloc_I_1d(&mat->indexL);
  monolis_dealloc_I_1d(&mat->itemU);
  monolis_dealloc_I_1d(&mat->itemL);
}

void monolis_mat_finalize_CSR(
  MONOLIS_MAT_CSR* mat)
{
  monolis_dealloc_I_1d(&mat->index);
  monolis_dealloc_I_1d(&mat->item);
}

void monolis_mat_finalize_CSC(
  MONOLIS_MAT_CSC* mat)
{
  monolis_dealloc_I_1d(&mat->index);
  monolis_dealloc_I_1d(&mat->item);
  monolis_dealloc_I_1d(&mat->perm);
}

void monolis_mat_finalize(
  MONOLIS_MAT* mat)
{
  mat->N = 0;
  mat->NP = 0;
  mat->NDOF = 0;
  monolis_dealloc_I_1d(&mat->n_dof_list);
  monolis_dealloc_I_1d(&mat->n_dof_index);
  monolis_dealloc_I_1d(&mat->n_dof_index2);
  monolis_mat_finalize_val_R(&mat->R);
  monolis_mat_finalize_val_C(&mat->C);
  monolis_mat_finalize_SCSR(&mat->SCSR);
  monolis_mat_finalize_CSR(&mat->CSR);
  monolis_mat_finalize_CSC(&mat->CSC);
}

void monolis_get_vec_size(
  int N,
  int NP,
  int NDOF,
  int NZ,
  int* n_dof_index,
  int* n_dof_index2,
  int* N_size,
  int* NP_size,
  int* NZ_size)
{
  if(NDOF == -1){
    *N_size  = n_dof_index[N];
    *NP_size = n_dof_index[NP];
    *NZ_size = n_dof_index2[NZ];
  } else {
    *N_size  = N *NDOF;
    *NP_size = NP*NDOF;
    *NZ_size = NZ*NDOF*NDOF;
  }
}

void monolis_get_mat_size(
  MONOLIS_MAT* mat,
  int* NZ,
  int* NZD,
  int* NZU,
  int* NZL)
{
  int NP   = mat->NP;
  int NDOF = mat->NDOF;

  if(NDOF == -1){
    int in = mat->CSR.index[NP];
    *NZ = mat->n_dof_index2[in];

    *NZD = 0;
    *NZU = 0;
    *NZL = 0;
    for (int i = 0; i < NP; ++i) {
      int jS = mat->CSR.index[i];
      int jE = mat->CSR.index[i + 1];
      for (int j = jS; j < jE; ++j) {
        in = mat->CSR.item[j];
        if (i == in) {
          *NZD += mat->n_dof_index[i] * mat->n_dof_index[in];
        } else if (i < in) {
          *NZU += mat->n_dof_index[i] * mat->n_dof_index[in];
        } else {
          *NZL += mat->n_dof_index[i] * mat->n_dof_index[in];
        }
      }
    }

  } else {
    *NZ = mat->CSR.index[NP];
    *NZ = *NZ*NDOF*NDOF;

    *NZU = 0;
    if(mat->SCSR.indexU != NULL){
      *NZU = mat->SCSR.indexU[NP];
    }

    *NZL = 0;
    if(mat->SCSR.indexL != NULL){
      *NZL = mat->SCSR.indexL[NP];
    }

    *NZD = NP  *NDOF*NDOF;
    *NZU = *NZU*NDOF*NDOF;
    *NZL = *NZL*NDOF*NDOF;
  }
}
