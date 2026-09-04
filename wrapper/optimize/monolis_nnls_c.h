/* monolis_nnls_c.h */
#ifndef MONOLIS_NNLS_H
#define MONOLIS_NNLS_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Non-Negetive Least Squares 関数（Lawson-Hanson 法）
 * @details 係数行列 A は列ブロック分散（各ランクが m x n_loc の列ブロックを保持し、
 * 全体の列順序はランク番号順に連結した順序）とし、右辺ベクトル b は全ランクに複製する。
 * 逐次実行する場合はセルフコミュニケータを渡す。
 * @param[in] A 係数行列（m x n_loc、列ブロック分散）
 * @param[in] b 右辺ベクトル（m、全ランク複製）
 * @param[out] x 解ベクトル（n_loc、自ランク保持列に対応）
 * @param[in] m 行列の行数
 * @param[in] n_loc 自ランクが保持する列数
 * @param[in] max_iter 最大反復回数
 * @param[in] tol 収束判定閾値（KKT 条件の判定と微小解成分の除外に使用）
 * @param[out] residual 残差
 * @param[in] comm MPI コミュニケータ
 * @ingroup opt
 */
void monolis_optimize_nnls_R(
  double** A,
  double*  b,
  double*  x,
  int      m,
  int      n_loc,
  int      max_iter,
  double   tol,
  double*  residual,
  int      comm);

void monolis_optimize_nnls_R_c_main(
  double*  A,
  double*  b,
  double*  x,
  int      m,
  int      n_loc,
  int      max_iter,
  double   tol,
  double*  residual,
  int      comm);

/**
 * @brief Non-Negetive Least Squares 関数（疎な解ベクトル）
 * @details 係数行列 A は列ブロック分散（各ランクが m x n_loc の列ブロックを保持し、
 * 全体の列順序はランク番号順に連結した順序）とし、右辺ベクトル b は全ランクに複製する。
 * 逐次実行する場合はセルフコミュニケータを渡す。
 * @param[in] A 係数行列（m x n_loc、列ブロック分散）
 * @param[in] b 右辺ベクトル（m、全ランク複製）
 * @param[out] x 解ベクトル（n_loc、自ランク保持列に対応）
 * @param[in] m 行列の行数
 * @param[in] n_loc 自ランクが保持する列数
 * @param[in] max_iter 最大反復回数
 * @param[in] tol_outer 外側反復の収束判定閾値（相対残差に対する閾値）
 * @param[in] tol_inner 内部 NNLS の閾値（この値以下の解成分を active 集合から除外）
 * @param[out] residual 残差
 * @param[in] comm MPI コミュニケータ
 * @ingroup opt
 */
void monolis_optimize_nnls_R_with_sparse_solution(
  double** A,
  double*  b,
  double*  x,
  int      m,
  int      n_loc,
  int      max_iter,
  double   tol_outer,
  double   tol_inner,
  double*  residual,
  int      comm);

void monolis_optimize_nnls_R_with_sparse_solution_c_main(
  double*  A,
  double*  b,
  double*  x,
  int      m,
  int      n_loc,
  int      max_iter,
  double   tol_outer,
  double   tol_inner,
  double*  residual,
  int      comm);

#ifdef __cplusplus
}
#endif

#endif
