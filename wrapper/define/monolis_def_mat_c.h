/* monolis_def_mat.h */
#ifndef MONOLIS_DEF_MAT_H
#define MONOLIS_DEF_MAT_H

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
  double _Complex* A;
  double _Complex* U;
  double _Complex* D;
  double _Complex* L;
  double _Complex* X;
  double _Complex* B;
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

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup def_mat_init
 */
void monolis_mat_initialize_val_R(
  MONOLIS_MAT_VAL_R* mat);

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup def_mat_init
 */
void monolis_mat_initialize_val_C(
  MONOLIS_MAT_VAL_C* mat);

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup def_mat_init
 */
void monolis_mat_initialize_SCSR(
  MONOLIS_MAT_SEPARATED_CSR* mat);

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup def_mat_init
 */
void monolis_mat_initialize_CSR(
  MONOLIS_MAT_CSR* mat);

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup def_mat_init
 */
void monolis_mat_initialize_CSC(
  MONOLIS_MAT_CSC* mat);

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup def_mat_init
 */
void monolis_mat_initialize(
  MONOLIS_MAT* mat);

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup def_mat_init
 */
void monolis_mat_finalize_val_R(
  MONOLIS_MAT_VAL_R* mat);

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup def_mat_init
 */
void monolis_mat_finalize_val_C(
  MONOLIS_MAT_VAL_C* mat);

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup def_mat_init
 */
void monolis_mat_finalize_SCSR(
  MONOLIS_MAT_SEPARATED_CSR* mat);

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup def_mat_init
 */
void monolis_mat_finalize_CSR(
  MONOLIS_MAT_CSR* mat);

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup def_mat_init
 */
void monolis_mat_finalize_CSC(
  MONOLIS_MAT_CSC* mat);

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup def_mat_init
 */
void monolis_mat_finalize(
  MONOLIS_MAT* mat);

#ifdef __cplusplus
}
#endif

#endif
