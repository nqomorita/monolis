/* monolis_vec_util.h */
#ifndef MONOLIS_VEC_UTIL_H
#define MONOLIS_VEC_UTIL_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief ベクトル配列コピー（整数型）
 * @param[in] n 領域の全計算点数
 * @param[in] n_dof 計算点あたりの自由度
 * @param[in] x ベクトル 1 (コピー元)
 * @param[out] y ベクトル 2 (コピー先)
 * @ingroup linalg
 */
void monolis_vec_copy_I(
  int  n,
  int  n_dof,
  int* x,
  int* y);

/**
 * @brief ベクトル配列コピー（実数型）
 * @param[in] n 領域の全計算点数
 * @param[in] n_dof 計算点あたりの自由度
 * @param[in] x ベクトル 1 (コピー元)
 * @param[out] y ベクトル 2 (コピー先)
 * @ingroup linalg
 */
void monolis_vec_copy_R(
  int     n,
  int     n_dof,
  double* x,
  double* y);

/**
 * @brief ベクトル配列コピー（複素数型）
 * @param[in] n 領域の全計算点数
 * @param[in] n_dof 計算点あたりの自由度
 * @param[in] x ベクトル 1 (コピー元)
 * @param[out] y ベクトル 2 (コピー先)
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
