/* monolis_spmat_handler.h */
#ifndef MONOLIS_SPMAT_HANDLER_UTIL_H
#define MONOLIS_SPMAT_HANDLER_UTIL_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup dev_matrix
 */
void monolis_set_scalar_to_sparse_matrix_R_c_main(
  int      n_node,
  int      nz,
  int      n_dof,
  int*     index,
  int*     item,
  double*  A,
  int      i,
  int      j,
  int      submat_i,
  int      submat_j,
  double   val);

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup dev_matrix
 */
void monolis_add_scalar_to_sparse_matrix_R_c_main(
  int      n_node,
  int      nz,
  int      n_dof,
  int*     index,
  int*     item,
  double*  A,
  int      i,
  int      j,
  int      submat_i,
  int      submat_j,
  double   val);

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup dev_matrix
 */
void monolis_get_scalar_from_sparse_matrix_R_c_main(
  int      n_node,
  int      nz,
  int      n_dof,
  int*     index,
  int*     item,
  double*  A,
  int      i,
  int      j,
  int      submat_i,
  int      submat_j,
  double*  val,
  int*     is_find);

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup dev_matrix
 */
void monolis_add_matrix_to_sparse_matrix_main_R_c_main(
  int     n_node,
  int     nz,
  int     n_dof,
  int*    index,
  int*    item,
  double* A,
  int     n_base1,
  int     n_base2,
  int*    connectivity1,
  int*    connectivity2,
  double* val);

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup dev_matrix
 */
void monolis_set_Dirichlet_bc_R_c_main(
  int      n_node,
  int      nz,
  int      n_dof,
  int*     index,
  int*     item,
  int*     indexR,
  int*     itemR,
  int*     permR,
  double*  A,
  double*  b,
  int      node_id,
  int      n_dof_bc,
  double   val);

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup dev_matrix
 */
void monolis_set_scalar_to_sparse_matrix_C_c_main(
  int             n_node,
  int             nz,
  int             n_dof,
  int*            index,
  int*            item,
  double _Complex* A,
  int             i,
  int             j,
  int             submat_i,
  int             submat_j,
  double _Complex  val);

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup dev_matrix
 */
void monolis_add_scalar_to_sparse_matrix_C_c_main(
  int             n_node,
  int             nz,
  int             n_dof,
  int*            index,
  int*            item,
  double _Complex* A,
  int             i,
  int             j,
  int             submat_i,
  int             submat_j,
  double _Complex  val);

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup dev_matrix
 */
void monolis_get_scalar_from_sparse_matrix_C_c_main(
  int             n_node,
  int             nz,
  int             n_dof,
  int*            index,
  int*            item,
  double _Complex* A,
  int             i,
  int             j,
  int             submat_i,
  int             submat_j,
  double _Complex* val,
  int*            is_find);

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup dev_matrix
 */
void monolis_add_matrix_to_sparse_matrix_main_C_c_main(
  int              n_node,
  int              nz,
  int              n_dof,
  int*             index,
  int*             item,
  double _Complex*  A,
  int              n_base1,
  int              n_base2,
  int*             connectivity1,
  int*             connectivity2,
  double _Complex*  val);

/**
 * @brief 1 次元整数配列のメモリ確保
 * @param[inout] var メモリ確保する配列
 * @ingroup dev_matrix
 */
void monolis_set_Dirichlet_bc_C_c_main(
  int             n_node,
  int             nz,
  int             n_dof,
  int*            index,
  int*            item,
  int*            indexR,
  int*            itemR,
  int*            permR,
  double _Complex* A,
  double _Complex* b,
  int             node_id,
  int             n_dof_bc,
  double _Complex  val);

#ifdef __cplusplus
}
#endif

#endif
