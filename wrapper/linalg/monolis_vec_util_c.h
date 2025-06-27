/* monolis_vec_util.h */
#ifndef MONOLIS_VEC_UTIL_H
#define MONOLIS_VEC_UTIL_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief ベクトル配列コピー（整数型）
 * @param[in] m 配列長さ
 * @param[in] x ベクトル 1 (コピー元)
 * @param[out] y ベクトル 2 (コピー先)
 * @ingroup linalg
 */
void monolis_vec_copy_I(
  int  m,
  int* x,
  int* y);

/**
 * @brief ベクトル配列コピー（実数型）
 * @param[in] m 配列長さ
 * @param[in] x ベクトル 1 (コピー元)
 * @param[out] y ベクトル 2 (コピー先)
 * @ingroup linalg
 */
void monolis_vec_copy_R(
  int     m,
  double* x,
  double* y);

/**
 * @brief ベクトル配列コピー（複素数型）
 * @param[in] m 配列長さ
 * @param[in] x ベクトル 1 (コピー元)
 * @param[out] y ベクトル 2 (コピー先)
 * @ingroup linalg
 */
void monolis_vec_copy_C(
  int              m,
  double _Complex* x,
  double _Complex* y);

#ifdef __cplusplus
}
#endif

#endif
