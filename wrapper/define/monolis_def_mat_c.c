#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void monolis_mat_initialize_val_R(
  MONOLIS_MAT_VAL_R* mat)
{

}

void monolis_mat_initialize_val_C(
  MONOLIS_MAT_VAL_C* mat)
{

}

void monolis_mat_initialize_SCSR(
  MONOLIS_MAT_SEPARATED_CSR* mat)
{

}

void monolis_mat_initialize_CSR(
  MONOLIS_MAT_CSR* mat)
{

}

void monolis_mat_initialize_CSC(
  MONOLIS_MAT_CSC* mat)
{

}

void monolis_mat_initialize(
  MONOLIS* mat)
{
  MONOLIS->N = 0;
  MONOLIS->NP = 0;
  MONOLIS->NDOF = 0;
  monolis_mat_initialize_val_R(MONOLIS->R);
  monolis_mat_initialize_val_C(MONOLIS->C);
  monolis_mat_initialize_SCSR(MONOLIS->SCSR);
  monolis_mat_initialize_CSR(MONOLIS->CSR);
  monolis_mat_initialize_CSC(MONOLIS->CSC);
}

void monolis_mat_finalize_val_R(
  MONOLIS_MAT_VAL_R* mat)
{

}

void monolis_mat_finalize_val_C(
  MONOLIS_MAT_VAL_C* mat)
{

}

void monolis_mat_finalize_SCSR(
  MONOLIS_MAT_SEPARATED_CSR* mat)
{

}

void monolis_mat_finalize_CSR(
  MONOLIS_MAT_CSR* mat)
{

}

void monolis_mat_finalize_CSC(
  MONOLIS_MAT_CSC* mat)
{

}

void monolis_mat_finalize(
  MONOLIS* mat)
{
  MONOLIS->N = 0;
  MONOLIS->NP = 0;
  MONOLIS->NDOF = 0;
  monolis_mat_finalize_val_R(MONOLIS->R);
  monolis_mat_finalize_val_C(MONOLIS->C);
  monolis_mat_finalize_SCSR(MONOLIS->SCSR);
  monolis_mat_finalize_CSR(MONOLIS->CSR);
  monolis_mat_finalize_CSC(MONOLIS->CSC);
}