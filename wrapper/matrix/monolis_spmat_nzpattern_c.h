/* monolis_spmat_nonzero_pattern.h */
#ifndef MONOLIS_SPMAT_NONZERO_PATTERN_H
#define MONOLIS_SPMAT_NONZERO_PATTERN_H

#ifdef __cplusplus
extern "C" {
#endif

#include "monolis_def_struc_c.h"

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
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
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup nzpattern
 */
void monolis_get_nonzero_pattern_by_connectivity_R(
  MONOLIS* mat,
  int      n_node,
  int      n_base,
  int      n_dof,
  int      n_elem,
  int*     conn_index,
  int*     conn_item);

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup nzpattern
 */
void monolis_get_nonzero_pattern_by_nodal_graph_R(
  MONOLIS* mat,
  int      n_node,
  int      n_dof,
  int*     index,
  int*     item);

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
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
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup nzpattern
 */
void monolis_get_nonzero_pattern_by_connectivity_C(
  MONOLIS* mat,
  int      n_node,
  int      n_base,
  int      n_dof,
  int      n_elem,
  int*     conn_index,
  int*     conn_item);

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
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
