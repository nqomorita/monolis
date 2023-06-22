/* monolis_inner_product.h */
#ifndef MONOLIS_INNER_PRODUCT_H
#define MONOLIS_INNER_PRODUCT_H

#ifdef __cplusplus
extern "C" {
#endif

#include "monolis_def_struc_c.h"

/**
 * @brief ベクトル内積（整数型）
 * @param[in] mat monolis 構造体
 * @param[in] com COM 構造体
 * @param[in] n_dof 計算点が持つ自由度
 * @param[in] x ベクトル 1
 * @param[in] y ベクトル 2
 * @param[out] sum 内積結果
 * @ingroup linalg
 */
void monolis_inner_product_I(
  MONOLIS*     mat,
  MONOLIS_COM* com,
  int          n_dof,
  int*         x,
  int*         y,
  int*         sum);

/**
 * @brief ベクトル内積（実数型）
 * @param[in] mat monolis 構造体
 * @param[in] com COM 構造体
 * @param[in] n_dof 計算点が持つ自由度
 * @param[in] x ベクトル 1
 * @param[in] y ベクトル 2
 * @param[out] sum 内積結果
 * @ingroup linalg
 */
void monolis_inner_product_R(
  MONOLIS*     mat,
  MONOLIS_COM* com,
  int          n_dof,
  double*      x,
  double*      y,
  double*      sum);

/**
 * @brief ベクトル内積（複素数型）
 * @param[in] mat monolis 構造体
 * @param[in] com COM 構造体
 * @param[in] n_dof 計算点が持つ自由度
 * @param[in] x ベクトル 1
 * @param[in] y ベクトル 2
 * @param[out] sum 内積結果
 * @ingroup linalg
 */
void monolis_inner_product_C(
  MONOLIS*        mat,
  MONOLIS_COM*    com,
  int             n_dof,
  double _Complex* x,
  double _Complex* y,
  double _Complex* sum);

/**
 * @brief ベクトル内積（整数型）
 * @param[in] n 計算点数
 * @param[in] n_dof 計算点あたりの自由度
 * @param[in] x ベクトル 1
 * @param[in] y ベクトル 2
 * @param[out] sum 内積結果
 * @param[in] comm MPI コミュニケータ
 * @ingroup dev_linalg
 */
void monolis_inner_product_I_c_main(
  int  n,
  int  n_dof,
  int* x,
  int* y,
  int* sum,
  int  comm);

/**
 * @brief ベクトル内積（実数型）
 * @param[in] n 計算点数
 * @param[in] n_dof 計算点あたりの自由度
 * @param[in] x ベクトル 1
 * @param[in] y ベクトル 2
 * @param[out] sum 内積結果
 * @param[in] comm MPI コミュニケータ
 * @ingroup dev_linalg
 */
void monolis_inner_product_R_c_main(
  int     n,
  int     n_dof,
  double* x,
  double* y,
  double* sum,
  int     comm);

/**
 * @brief ベクトル内積（複素数型）
 * @param[in] n 計算点数
 * @param[in] n_dof 計算点あたりの自由度
 * @param[in] x ベクトル 1
 * @param[in] y ベクトル 2
 * @param[out] sum 内積結果
 * @param[in] comm MPI コミュニケータ
 * @ingroup dev_linalg
 */
void monolis_inner_product_C_c_main(
  int             n,
  int             n_dof,
  double _Complex* x,
  double _Complex* y,
  double _Complex* sum,
  int             comm);

#ifdef __cplusplus
}
#endif

#endif
