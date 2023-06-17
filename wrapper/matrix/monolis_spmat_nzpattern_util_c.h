/* monolis_spmat_nonzero_pattern_util.h */
#ifndef MONOLIS_SPMAT_NONZERO_PATTERN_UTIL_H
#define MONOLIS_SPMAT_NONZERO_PATTERN_UTIL_H

#ifdef __cplusplus
extern "C" {
#endif

#include "monolis_def_mat_c.h"

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup dev_matrix
 */
void monolis_get_nonzero_pattern_by_nodal_graph_main(
  MONOLIS_MAT* mat,
  int          n_node,
  int          n_dof,
  int*         index,
  int*         item);

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup dev_matrix
 */
void monolis_alloc_nonzero_pattern_mat_val_R(
  MONOLIS_MAT* mat);

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup dev_matrix
 */
void monolis_alloc_nonzero_pattern_mat_val_C(
  MONOLIS_MAT* mat);

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
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
