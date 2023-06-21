/* monolis_def_solver_util.h */
#ifndef MONOLIS_DEF_SOLVER_UTIL_H
#define MONOLIS_DEF_SOLVER_UTIL_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include "monolis_def_struc_c.h"

/**
 * @brief ソルバの設定
 * @param[inout] mat monolis 構造体
 * @param[in] prm パラメータ
 * @ingroup param
 */
void monolis_set_method(
  MONOLIS* mat,
  int      param);

/**
 * @brief 前処理の設定
 * @param[inout] mat monolis 構造体
 * @param[in] prm パラメータ
 * @ingroup param
 */
void monolis_set_precond(
  MONOLIS* mat,
  int      param);

/**
 * @brief 最大反復回数の設定
 * @param[inout] mat monolis 構造体
 * @param[in] prm パラメータ
 * @ingroup param
 */
void monolis_set_maxiter(
  MONOLIS* mat,
  int      param);

/**
 * @brief 現在の反復回数の取得
 * @param[in] mat monolis 構造体
 * @param[out] prm パラメータ
 * @ingroup param
 */
void monolis_get_converge_iter(
  MONOLIS* mat,
  int*     param);

/**
 * @brief エラー番号の取得
 * @param[in] mat monolis 構造体
 * @param[out] prm パラメータ
 * @ingroup param
 */
void monolis_get_error_tag(
  MONOLIS* mat,
  int*     param);

/**
 * @brief 解ベクトル初期化の有無の設定
 * @param[inout] mat monolis 構造体
 * @param[in] prm パラメータ
 * @ingroup param
 */
void monolis_set_init_x(
  MONOLIS* mat,
  bool     param);

/**
 * @brief 対称行列向け処理の有無の設定
 * @param[inout] mat monolis 構造体
 * @param[in] prm パラメータ
 * @ingroup param
 */
void monolis_set_sym_matrix(
  MONOLIS* mat,
  bool     param);

/**
 * @brief デバッグ出力の有無の設定
 * @param[inout] mat monolis 構造体
 * @param[in] prm パラメータ
 * @ingroup param
 */
void monolis_set_debug(
  MONOLIS* mat,
  bool     param);

/**
 * @brief 詳細な計算時間測定の有無の設定
 * @param[inout] mat monolis 構造体
 * @param[in] prm パラメータ
 * @ingroup param
 */
void monolis_set_performance_measurement(
  MONOLIS* mat,
  bool     param);

/**
 * @brief 行列対角成分確認の有無の設定
 * @param[inout] mat monolis 構造体
 * @param[in] prm パラメータ
 * @ingroup param
 */
void monolis_set_check_diag(
  MONOLIS* mat,
  bool     param);

/**
 * @brief 前処理情報保存の有無の設定
 * @param[inout] mat monolis 構造体
 * @param[in] prm パラメータ
 * @ingroup param
 */
void monolis_set_prec_stored(
  MONOLIS* mat,
  bool     param);

/**
 * @brief 反復回数と残差履歴の表示の設定
 * @param[inout] mat monolis 構造体
 * @param[in] prm パラメータ
 * @ingroup param
 */
void monolis_show_iterlog(
  MONOLIS* mat,
  bool     param);

/**
 * @brief 詳細な計算時間の表示の設定
 * @param[inout] mat monolis 構造体
 * @param[in] prm パラメータ
 * @ingroup param
 */
void monolis_show_timelog(
  MONOLIS* mat,
  bool     param);

/**
 * @brief ソルバ収束後のサマリの表示の設定
 * @param[inout] mat monolis 構造体
 * @param[in] prm パラメータ
 * @ingroup param
 */
void monolis_show_summary(
  MONOLIS* mat,
  bool     param);

/**
 * @brief 計算時間の統計的処理結果の表示の設定
 * @param[inout] mat monolis 構造体
 * @param[in] prm パラメータ
 * @ingroup param
 */
void monolis_show_timelog_statistics(
  MONOLIS* mat,
  bool     param);

/**
 * @brief 収束判定閾値の設定
 * @param[inout] mat monolis 構造体
 * @param[in] prm パラメータ
 * @ingroup param
 */
void monolis_set_tolerance(
  MONOLIS* mat,
  double   val);

/**
 * @brief 現在の残差の取得
 * @param[in] mat monolis 構造体
 * @param[inout] prm パラメータ
 * @ingroup param
 */
void monolis_get_converge_residual(
  MONOLIS* mat,
  double*  val);

/**
 * @brief ソルバの全計算時間の取得
 * @param[in] mat monolis 構造体
 * @param[inout] prm パラメータ
 * @ingroup param
 */
void monolis_get_time_solver(
  MONOLIS* mat,
  double*  val);

/**
 * @brief 前処理時間（生成時間）の取得
 * @param[in] mat monolis 構造体
 * @param[inout] prm パラメータ
 * @ingroup param
 */
void monolis_get_time_preparing(
  MONOLIS* mat,
  double*  val);

/**
 * @brief 疎行列ベクトル積時間の取得
 * @param[in] mat monolis 構造体
 * @param[inout] prm パラメータ
 * @ingroup param
 */
void monolis_get_time_spmv(
  MONOLIS* mat,
  double*  val);

/**
 * @brief ベクトル内積時間の取得
 * @param[in] mat monolis 構造体
 * @param[inout] prm パラメータ
 * @ingroup param
 */
void monolis_get_time_inner_product(
  MONOLIS* mat,
  double*  val);

/**
 * @brief 前処理時間（適用時間）の取得
 * @param[in] mat monolis 構造体
 * @param[inout] prm パラメータ
 * @ingroup param
 */
void monolis_get_time_precondition(
  MONOLIS* mat,
  double*  val);

/**
 * @brief ベクトル内積の通信時間の取得
 * @param[in] mat monolis 構造体
 * @param[inout] prm パラメータ
 * @ingroup param
 */
void monolis_get_time_comm_inner_product(
  MONOLIS* mat,
  double*  val);

/**
 * @brief 疎行列ベクトル積の通信時間の取得
 * @param[in] mat monolis 構造体
 * @param[inout] prm パラメータ
 * @ingroup param
 */
void monolis_get_time_comm_spmv(
  MONOLIS* mat,
  double*  val);

#ifdef __cplusplus
}
#endif

#endif
