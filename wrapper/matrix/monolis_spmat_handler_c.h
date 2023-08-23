/* monolis_spmat_handler.h */
#ifndef MONOLIS_SPMAT_HANDLER_H
#define MONOLIS_SPMAT_HANDLER_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief スカラ値を疎行列に設定（実数型）
 * @param[inout] mat monolis 構造体
 * @param[in] i 行番号
 * @param[in] j 列番号
 * @param[in] submat_i ブロック中の行番号
 * @param[in] submat_j ブロック中の列番号
 * @param[in] val 行列値
 * @ingroup matrix
 */
void monolis_set_scalar_to_sparse_matrix_R(
  MONOLIS* mat,
  int      i,
  int      j,
  int      submat_i,
  int      submat_j,
  double   val);

/**
 * @brief スカラ値を疎行列に足込（実数型）
 * @param[inout] mat monolis 構造体
 * @param[in] i 行番号
 * @param[in] j 列番号
 * @param[in] submat_i ブロック中の行番号
 * @param[in] submat_j ブロック中の列番号
 * @param[in] val 行列値
 * @ingroup matrix
 */
void monolis_add_scalar_to_sparse_matrix_R(
  MONOLIS* mat,
  int      i,
  int      j,
  int      submat_i,
  int      submat_j,
  double   val);

/**
 * @brief スカラ値を疎行列から取得（実数型）
 * @param[inout] mat monolis 構造体
 * @param[in] i 行番号
 * @param[in] j 列番号
 * @param[in] submat_i ブロック中の行番号
 * @param[in] submat_j ブロック中の列番号
 * @param[out] val 行列値
 * @param[out] is_find 取得判定フラグ
 * @ingroup matrix
 */
void monolis_get_scalar_from_sparse_matrix_R(
  MONOLIS* mat,
  int      i,
  int      j,
  int      submat_i,
  int      submat_j,
  double*  val,
  bool*    is_find);

/**
 * @brief 行列を疎行列に足込（実数型）
 * @param[inout] mat monolis 構造体
 * @param[in] n_base 節点数
 * @param[in] connectivity 要素コネクティビティ
 * @param[in] val 行列値
 * @ingroup matrix
 */
void monolis_add_matrix_to_sparse_matrix_R(
  MONOLIS* mat,
  int      n_base,
  int*     connectivity,
  double** val);

/**
 * @brief 行列を疎行列の非対角部分に足込（実数型）
 * @param[inout] mat monolis 構造体
 * @param[in] n_base1 節点数
 * @param[in] n_base2 節点数
 * @param[in] connectivity1 要素コネクティビティ
 * @param[in] connectivity2 要素コネクティビティ
 * @param[in] val 行列値
 * @ingroup matrix
 */
void monolis_add_matrix_to_sparse_matrix_offdiag_R(
  MONOLIS* mat,
  int      n_base1,
  int      n_base2,
  int*     connectivity1,
  int*     connectivity2,
  double** val);

/**
 * @brief BCSR 形式の疎行列を直接設定（実数型）
 * @param[inout] mat monolis 構造体
 * @param[in] n 内部節点数
 * @param[in] np 全自由度数
 * @param[in] n_dof ブロックサイズ
 * @param[in] nz 非零要素数
 * @param[in] A 行列値配列
 * @param[in] index index 配列
 * @param[in] item item 配列
 * @ingroup matrix
 */
void monolis_set_matrix_BCSR_R(
  MONOLIS* mat,
  int      n,
  int      np,
  int      n_dof,
  int      nz,
  double*  A,
  int*     index,
  int*     item);

/**
 * @brief BCSR 形式の疎行列の行列値を直接設定（実数型）
 * @param[inout] mat monolis 構造体
 * @param[in] n_dof ブロックサイズ
 * @param[in] nz 非零要素数
 * @param[in] A 行列値配列
 * @ingroup matrix
 */
void monolis_set_matrix_BCSR_mat_val_R(
  MONOLIS* mat,
  int      n_dof,
  int      nz,
  double*  A);

/**
 * @brief 境界条件処理（実数型）
 * @param[inout] mat monolis 構造体
 * @param[inout] b 右辺ベクトル
 * @param[in] node_id 自由度番号
 * @param[in] n_dof_bc ブロック番号
 * @param[in] val 境界条件の設定値
 * @ingroup matrix
 */
void monolis_set_Dirichlet_bc_R(
  MONOLIS* mat,
  double*  b,
  int      node_id,
  int      n_dof_bc,
  double   val);

/**
 * @brief スカラ値を疎行列に設定（複素数型）
 * @param[inout] mat monolis 構造体
 * @param[in] i 行番号
 * @param[in] j 列番号
 * @param[in] submat_i ブロック中の行番号
 * @param[in] submat_j ブロック中の列番号
 * @param[out] val 行列値
 * @ingroup matrix
 */
void monolis_set_scalar_to_sparse_matrix_C(
  MONOLIS*       mat,
  int            i,
  int            j,
  int            submat_i,
  int            submat_j,
  double _Complex val);

/**
 * @brief スカラ値を疎行列に足込（複素数型）
 * @param[inout] mat monolis 構造体
 * @param[in] i 行番号
 * @param[in] j 列番号
 * @param[in] submat_i ブロック中の行番号
 * @param[in] submat_j ブロック中の列番号
 * @param[out] val 行列値
 * @ingroup matrix
 */
void monolis_add_scalar_to_sparse_matrix_C(
  MONOLIS*       mat,
  int            i,
  int            j,
  int            submat_i,
  int            submat_j,
  double _Complex val);

/**
 * @brief スカラ値を疎行列から取得（複素数型）
 * @param[inout] mat monolis 構造体
 * @param[in] i 行番号
 * @param[in] j 列番号
 * @param[in] submat_i ブロック中の行番号
 * @param[in] submat_j ブロック中の列番号
 * @param[out] val 行列値
 * @param[out] is_find 取得判定フラグ
 * @ingroup matrix
 */
void monolis_get_scalar_from_sparse_matrix_C(
  MONOLIS*        mat,
  int             i,
  int             j,
  int             submat_i,
  int             submat_j,
  double _Complex* val,
  bool*           is_find);

/**
 * @brief 行列を疎行列に足込（複素数型）
 * @param[inout] mat monolis 構造体
 * @param[in] n_base 節点数
 * @param[in] connectivity 要素コネクティビティ
 * @param[in] val 行列値
 * @ingroup matrix
 */
void monolis_add_matrix_to_sparse_matrix_C(
  MONOLIS*         mat,
  int              n_base,
  int*             connectivity,
  double _Complex** val);

/**
 * @brief 行列を疎行列の非対角部分に足込（複素数型）
 * @param[inout] mat monolis 構造体
 * @param[in] n_base1 節点数
 * @param[in] n_base2 節点数
 * @param[in] connectivity1 要素コネクティビティ
 * @param[in] connectivity2 要素コネクティビティ
 * @param[in] val 行列値
 * @ingroup matrix
 */
void monolis_add_matrix_to_sparse_matrix_offdiag_C(
  MONOLIS*         mat,
  int              n_base1,
  int              n_base2,
  int*             connectivity1,
  int*             connectivity2,
  double _Complex** val);

/**
 * @brief BCSR 形式の疎行列を直接設定（複素数型）
 * @param[inout] mat monolis 構造体
 * @param[in] n 内部節点数
 * @param[in] np 全自由度数
 * @param[in] n_dof ブロックサイズ
 * @param[in] nz 非零要素数
 * @param[in] A 行列値配列
 * @param[in] index index 配列
 * @param[in] item item 配列
 * @ingroup matrix
 */
void monolis_set_matrix_BCSR_C(
  MONOLIS*        mat,
  int             n,
  int             np,
  int             n_dof,
  int             nz,
  double _Complex* A,
  int*            index,
  int*            item);

/**
 * @brief BCSR 形式の疎行列の行列値を直接設定（複素数型）
 * @param[inout] mat monolis 構造体
 * @param[in] n_dof ブロックサイズ
 * @param[in] nz 非零要素数
 * @param[in] A 行列値配列
 * @ingroup matrix
 */
void monolis_set_matrix_BCSR_mat_val_C(
  MONOLIS*        mat,
  int             n_dof,
  int             nz,
  double _Complex* A);

/**
 * @brief 境界条件処理（実数型）
 * @param[inout] mat monolis 構造体
 * @param[inout] b 右辺ベクトル
 * @param[in] node_id 自由度番号
 * @param[in] n_dof_bc ブロック番号
 * @param[in] val 境界条件の設定値
 * @ingroup matrix
 */
void monolis_set_Dirichlet_bc_C(
  MONOLIS*        mat,
  double _Complex* b,
  int             node_id,
  int             n_dof_bc,
  double _Complex  val);

#ifdef __cplusplus
}
#endif

#endif
