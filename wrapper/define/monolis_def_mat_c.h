/* monolis_def_struc.h */
#ifndef MONOLIS_DEF_SOLVER_H
#define MONOLIS_DEF_SOLVER_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  double* A;
  double* U;
  double* D;
  double* L;
  double* X;
  double* B;
} MONOLIS_MAT_VAL_R;

typedef struct {
  double complex* A;
  double complex* U;
  double complex* D;
  double complex* L;
  double complex* X;
  double complex* B;
} MONOLIS_MAT_VAL_C;

typedef struct {
  int* indexU;
  int* indexL;
  int* itemU;
  int* itemL;
} MONOLIS_MAT_SEPARATED_CSR;

typedef struct {
  int* index;
  int* item;
} MONOLIS_MAT_CSR;

typedef struct {
  int* index;
  int* item;
  int* perm;
} MONOLIS_MAT_CSC;

typedef struct {
  int N;
  int NP;
  int NDOF;
  MONOLIS_MAT_VAL_R R;
  MONOLIS_MAT_VAL_C C;
  MONOLIS_MAT_SEPARATED_CSR SCSR;
  MONOLIS_MAT_CSR CSR;
  MONOLIS_MAT_CSC CSC;
} MONOLIS_MAT;

void monolis_mat_initialize_val_R(
  MONOLIS_MAT_VAL_R* mat);

void monolis_mat_initialize_val_C(
  MONOLIS_MAT_VAL_C* mat);

void monolis_mat_initialize_SCSR(
  MONOLIS_MAT_SEPARATED_CSR* mat);

void monolis_mat_initialize_CSR(
  MONOLIS_MAT_CSR* mat);

void monolis_mat_initialize_CSC(
  MONOLIS_MAT_CSC* mat);

void monolis_mat_initialize(
  MONOLIS* mat);

void monolis_mat_finalize_val_R(
  MONOLIS_MAT_VAL_R* mat);

void monolis_mat_finalize_val_C(
  MONOLIS_MAT_VAL_C* mat);

void monolis_mat_finalize_SCSR(
  MONOLIS_MAT_SEPARATED_CSR* mat);

void monolis_mat_finalize_CSR(
  MONOLIS_MAT_CSR* mat);

void monolis_mat_finalize_CSC(
  MONOLIS_MAT_CSC* mat);

void monolis_mat_finalize(
  MONOLIS* mat);

#ifdef __cplusplus
}
#endif

#endif
