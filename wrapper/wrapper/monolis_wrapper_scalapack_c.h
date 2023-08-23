/* monolis_wrap_scalapack.h */
#ifndef MONOLIS_WRAP_SCALAPACK_H
#define MONOLIS_WRAP_SCALAPACK_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief PDGESVD 関数（実数型）
 * @param[in] N_loc 行列の大きさ（行数 N）
 * @param[in] M 行列の大きさ（列数 M）
 * @param[in] A 入力行列（N_loc x M）
 * @param[out] S 左特異行列（N_loc x P）
 * @param[out] V 特異値（P）
 * @param[out] D 右特異行列（P x M）
 * @param[in] comm コミュニケータ
 * @ingroup wrapper
 */
void monolis_scalapack_gesvd_R(
  int      N_loc,
  int      M,
  double** A,
  double** S,
  double*  V,
  double** D,
  int      comm);

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[in] N_loc 行列の大きさ（行数 N）
 * @param[in] M 行列の大きさ（列数 M）
 * @param[in] P 行列の大きさの最小値
 * @param[in] A 入力行列（N_loc x M）
 * @param[out] S 左特異行列（N_loc x P）
 * @param[out] V 特異値（P）
 * @param[out] D 右特異行列（P x M）
 * @param[in] comm コミュニケータ
 * @ingroup wrapper
 */
void monolis_scalapack_gesvd_R_c_main(
  int      N_loc,
  int      M,
  int      P,
  double*  A,
  double*  S,
  double*  V,
  double*  D,
  int      comm);

#ifdef __cplusplus
}
#endif

#endif
