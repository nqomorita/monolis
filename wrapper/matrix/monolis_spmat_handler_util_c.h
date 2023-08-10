/* monolis_spmat_handler.h */
#ifndef MONOLIS_SPMAT_HANDLER_UTIL_H
#define MONOLIS_SPMAT_HANDLER_UTIL_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief スカラ値を疎行列に設定（メイン関数、実数型）
 * @param[in] n_node 
 * @param[in] nz
 * @param[in] n_dof ブロック自由度
 * @param[in] index index 配列
 * @param[in] item item 配列
 * @param[inout] A 係数行列
 * @param[in] i 行番号
 * @param[in] j 列番号
 * @param[in] submat_i ブロック中の行番号
 * @param[in] submat_j ブロック中の列番号
 * @param[in] val 設定値
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
 * @brief スカラ値を疎行列に足込（メイン関数、実数型）
 * @param[in] n_node 
 * @param[in] nz
 * @param[in] n_dof ブロック自由度
 * @param[in] index index 配列
 * @param[in] item item 配列
 * @param[inout] A 係数行列
 * @param[in] i 行番号
 * @param[in] j 列番号
 * @param[in] submat_i ブロック中の行番号
 * @param[in] submat_j ブロック中の列番号
 * @param[in] val 設定値
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
 * @brief スカラ値を疎行列から取得（メイン関数、実数型）
 * @param[in] n_node 
 * @param[in] nz
 * @param[in] n_dof ブロック自由度
 * @param[in] index index 配列
 * @param[in] item item 配列
 * @param[in] A 係数行列
 * @param[in] i 行番号
 * @param[in] j 列番号
 * @param[in] submat_i ブロック中の行番号
 * @param[in] submat_j ブロック中の列番号
 * @param[out] val 設定値
 * @param[out] is_find 取得判定フラグ
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
 * @brief 行列値を疎行列に足込（メイン関数、実数型）
 * @param[in] n_node 
 * @param[in] nz
 * @param[in] n_dof ブロック自由度
 * @param[in] index index 配列
 * @param[in] item item 配列
 * @param[inout] A 係数行列
 * @param[in] n_base1
 * @param[in] n_base2
 * @param[in] connectivity1
 * @param[in] connectivity2
 * @param[in] val 設定値
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
 * @brief 境界条件処理（実数型、メイン関数）
 * @param[in] n_node 
 * @param[in] nz
 * @param[in] n_dof ブロック自由度
 * @param[in] index index 配列
 * @param[in] item item 配列
 * @param[in] indexR index 配列（CSC 形式）
 * @param[in] itemR item 配列（CSC 形式）
 * @param[in] permR 行列成分の CSC 形式との置換ベクトル
 * @param[inout] A 係数行列
 * @param[inout] b 右辺ベクトル
 * @param[in] node_id 自由度番号
 * @param[in] n_dof_bc ブロック番号
 * @param[in] val 設定値
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
 * @brief スカラ値を疎行列に設定（メイン関数、複素数型）
 * @param[in] n_node 
 * @param[in] nz
 * @param[in] n_dof ブロック自由度
 * @param[in] index index 配列
 * @param[in] item item 配列
 * @param[inout] A 係数行列
 * @param[in] i 行番号
 * @param[in] j 列番号
 * @param[in] submat_i ブロック中の行番号
 * @param[in] submat_j ブロック中の列番号
 * @param[in] val 設定値
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
 * @brief スカラ値を疎行列に足込（メイン関数、複素数型）
 * @param[in] n_node 
 * @param[in] nz
 * @param[in] n_dof ブロック自由度
 * @param[in] index index 配列
 * @param[in] item item 配列
 * @param[inout] A 係数行列
 * @param[in] i 行番号
 * @param[in] j 列番号
 * @param[in] submat_i ブロック中の行番号
 * @param[in] submat_j ブロック中の列番号
 * @param[in] val 設定値
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
 * @brief スカラ値を疎行列から取得（メイン関数、複素数型）
 * @param[in] n_node 
 * @param[in] nz
 * @param[in] n_dof ブロック自由度
 * @param[in] index index 配列
 * @param[in] item item 配列
 * @param[in] A 係数行列
 * @param[in] i 行番号
 * @param[in] j 列番号
 * @param[in] submat_i ブロック中の行番号
 * @param[in] submat_j ブロック中の列番号
 * @param[out] val 設定値
 * @param[out] is_find 取得判定フラグ
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
 * @brief 行列値を疎行列に足込（メイン関数、複素数型）
 * @param[in] n_node 
 * @param[in] nz
 * @param[in] n_dof ブロック自由度
 * @param[in] index index 配列
 * @param[in] item item 配列
 * @param[inout] A 係数行列
 * @param[in] n_base1
 * @param[in] n_base2
 * @param[in] connectivity1
 * @param[in] connectivity2
 * @param[in] val 設定値
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
 * @brief 境界条件処理（複素数型、メイン関数）
 * @param[in] n_node 
 * @param[in] nz
 * @param[in] n_dof ブロック自由度
 * @param[in] index index 配列
 * @param[in] item item 配列
 * @param[in] indexR index 配列（CSC 形式）
 * @param[in] itemR item 配列（CSC 形式）
 * @param[in] permR 行列成分の CSC 形式との置換ベクトル
 * @param[inout] A 係数行列
 * @param[inout] b 右辺ベクトル
 * @param[in] node_id 自由度番号
 * @param[in] n_dof_bc ブロック番号
 * @param[in] val 設定値
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
