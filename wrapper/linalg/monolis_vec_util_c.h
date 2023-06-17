/* monolis_vec_util.h */
#ifndef MONOLIS_VEC_UTIL_H
#define MONOLIS_VEC_UTIL_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup linalg
 */
void monolis_vec_copy_I(
  int  n,
  int  n_dof,
  int* x,
  int* y);

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup linalg
 */
void monolis_vec_copy_R(
  int     n,
  int     n_dof,
  double* x,
  double* y);

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup linalg
 */
void monolis_vec_copy_C(
  int              n,
  int              n_dof,
  double _Complex* x,
  double _Complex* y);

#ifdef __cplusplus
}
#endif

#endif
