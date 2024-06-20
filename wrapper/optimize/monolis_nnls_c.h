/* monolis_nnls_c.h */
#ifndef MONOLIS_NNLS_H
#define MONOLIS_NNLS_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Non-Negetive Least Squares 関数
 * @param[in] A 係数行列（m x n）
 * @param[in] b 右辺ベクトル（m）
 * @param[out] x 解ベクトル（n）
 * @param[in] m 行列の行数
 * @param[in] n 行列の列数
 * @param[in] max_iter 最大反復回数
 * @param[in] tol 収束判定閾値
 * @param[out] residual 残差
 * @param[in] comm MPI コミュニケータ
 * @ingroup opt
 */
void monolis_optimize_nnls_R(
  double** A,
  double*  b,
  double*  x,
  int      m,
  int      n,
  int      max_iter,
  double   tol,
  double*  residual,
  int      comm);

void monolis_optimize_nnls_R_c_main(
  double*  A,
  double*  b,
  double*  x,
  int      m,
  int      n,
  int      max_iter,
  double   tol,
  double*  residual,
  int      comm);

/**
 * @brief Non-Negetive Least Squares 関数（疎な解ベクトル）
 * @param[in] A 係数行列（m x n）
 * @param[in] b 右辺ベクトル（m）
 * @param[out] x 解ベクトル（n）
 * @param[in] m 行列の行数
 * @param[in] n 行列の列数
 * @param[in] max_iter 最大反復回数
 * @param[in] tol 収束判定閾値
 * @param[out] residual 残差
 * @param[in] comm MPI コミュニケータ
 * @ingroup opt
 */
void monolis_optimize_nnls_R_with_sparse_solution(
  double** A,
  double*  b,
  double*  x,
  int      m,
  int      n,
  int      max_iter,
  double   tol,
  double*  residual,
  int      comm);

void monolis_optimize_nnls_R_with_sparse_solution_c_main(
  double*  A,
  double*  b,
  double*  x,
  int      m,
  int      n,
  int      max_iter,
  double   tol,
  double*  residual,
  int      comm);

#ifdef __cplusplus
}
#endif

#endif
