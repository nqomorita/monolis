/* monolis_spmat_nonzero_pattern_util.h */
#ifndef MONOLIS_SPMAT_NONZERO_PATTERN_UTIL_H
#define MONOLIS_SPMAT_NONZERO_PATTERN_UTIL_H

#ifdef __cplusplus
extern "C" {
#endif

#include "monolis_def_mat_c.h"

/**
 * @brief 節点グラフから疎行列パターンを決定（メイン関数）
 * @param[inout] mat monolis 構造体
 * @param[in] n_node 節点数
 * @param[in] n_dof 計算点が持つ自由度
 * @param[in] index index 配列
 * @param[in] item item 配列
 * @ingroup dev_matrix
 */
void monolis_get_nonzero_pattern_by_nodal_graph_main(
  MONOLIS_MAT* mat,
  int          n_node,
  int          n_dof,
  int*         index,
  int*         item);

/**
 * @brief 疎行列の行列成分のメモリ確保（実数型）
 * @param[inout] mat monolis 構造体
 * @ingroup dev_matrix
 */
void monolis_alloc_nonzero_pattern_mat_val_R(
  MONOLIS_MAT* mat);

/**
 * @brief 疎行列の行列成分のメモリ確保（複素数型）
 * @param[inout] mat monolis 構造体
 * @ingroup dev_matrix
 */
void monolis_alloc_nonzero_pattern_mat_val_C(
  MONOLIS_MAT* mat);

/**
 * @brief CSR 形式から CSC 形式のデータを生成
 * @param[in] ncol 行数
 * @param[in] nrow 列数
 * @param[in] nz 非零要素数
 * @param[in] index index 配列（CSR 形式）
 * @param[in] item item 配列（CSR 形式）
 * @param[out] indexR index 配列（CSC 形式）
 * @param[inout] itemR index 配列（CSC 形式）
 * @param[inout] permR index 配列（CSC 形式）
 * @ingroup dev_matrix
 */
void monolis_get_CSC_format(
  int  ncol,
  int  nrow,
  int  nz,
  int* index,
  int* item,
  int* indexR,
  int* itemR,
  int* permR);

#ifdef __cplusplus
}
#endif

#endif
