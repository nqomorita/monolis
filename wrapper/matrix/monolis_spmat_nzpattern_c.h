/* monolis_spmat_nonzero_pattern.h */
#ifndef MONOLIS_SPMAT_NONZERO_PATTERN_H
#define MONOLIS_SPMAT_NONZERO_PATTERN_H

#ifdef __cplusplus
extern "C" {
#endif

#include "monolis_def_struc_c.h"

/**
 * @brief 単一メッシュデータから疎行列パターンを決定（実数型）
 * @param[inout] mat monolis 構造体
 * @param[in] n_node 節点数
 * @param[in] n_base 要素を構成する形状関数の数
 * @param[in] n_dof 計算点が持つ自由度
 * @param[in] n_elem 要素数
 * @param[in] elem 要素コネクティビティ
 * @ingroup nzpattern
 */
void monolis_get_nonzero_pattern_by_simple_mesh_R(
  MONOLIS* mat,
  int      n_node,
  int      n_base,
  int      n_dof,
  int      n_elem,
  int**    elem);

/**
 * @brief コネクティビティグラフから疎行列パターンを決定（実数型）
 * @param[inout] mat monolis 構造体
 * @param[in] n_node 節点数
 * @param[in] n_base 要素を構成する形状関数の数
 * @param[in] n_dof 計算点が持つ自由度
 * @param[in] n_elem 要素数
 * @param[in] conn_index コネクティビティ index 配列
 * @param[in] conn_item コネクティビティ item 配列
 * @ingroup nzpattern
 */
void monolis_get_nonzero_pattern_by_connectivity_R(
  MONOLIS* mat,
  int      n_node,
  int      n_dof,
  int      n_elem,
  int*     conn_index,
  int*     conn_item);

/**
 * @brief 節点グラフから疎行列パターンを決定（実数型）
 * @param[inout] mat monolis 構造体
 * @param[in] n_node 節点数
 * @param[in] n_dof 計算点が持つ自由度
 * @param[in] index index 配列
 * @param[in] item item 配列
 * @ingroup nzpattern
 */
void monolis_get_nonzero_pattern_by_nodal_graph_R(
  MONOLIS* mat,
  int      n_node,
  int      n_dof,
  int*     index,
  int*     item);

/**
 * @brief 単一メッシュデータから疎行列パターンを決定（複素数型）
 * @param[inout] mat monolis 構造体
 * @param[in] n_node 節点数
 * @param[in] n_base 要素を構成する形状関数の数
 * @param[in] n_dof 計算点が持つ自由度
 * @param[in] n_elem 要素数
 * @param[in] elem 要素コネクティビティ
 * @ingroup nzpattern
 */
void monolis_get_nonzero_pattern_by_simple_mesh_C(
  MONOLIS* mat,
  int      n_node,
  int      n_base,
  int      n_dof,
  int      n_elem,
  int**    elem);

/**
 * @brief コネクティビティグラフから疎行列パターンを決定（複素数型）
 * @param[inout] mat monolis 構造体
 * @param[in] n_node 節点数
 * @param[in] n_base 要素を構成する形状関数の数
 * @param[in] n_dof 計算点が持つ自由度
 * @param[in] n_elem 要素数
 * @param[in] conn_index コネクティビティ index 配列
 * @param[in] conn_item コネクティビティ item 配列
 * @ingroup nzpattern
 */
void monolis_get_nonzero_pattern_by_connectivity_C(
  MONOLIS* mat,
  int      n_node,
  int      n_dof,
  int      n_elem,
  int*     conn_index,
  int*     conn_item);

/**
 * @brief 節点グラフから疎行列パターンを決定（複素数型）
 * @param[inout] mat monolis 構造体
 * @param[in] n_node 節点数
 * @param[in] n_dof 計算点が持つ自由度
 * @param[in] index index 配列
 * @param[in] item item 配列
 * @ingroup nzpattern
 */
void monolis_get_nonzero_pattern_by_nodal_graph_C(
  MONOLIS* mat,
  int      n_node,
  int      n_dof,
  int*     index,
  int*     item);

#ifdef __cplusplus
}
#endif

#endif
