/* monolis_opt_cvae_c.h */
#ifndef MONOLIS_OPT_CVAE_H
#define MONOLIS_OPT_CVAE_H

#include "monolis_opt_vae_c.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Conditional VAE モデル構造体 (C 側ハンドル)
 * @details Fortran 側の monolis_opt_cvae_t を指す不透明ハンドル net と、
 *          利便のための次元情報を保持する。データ配列はすべて単精度 (float)、
 *          Fortran と同じ列優先 (column-major) で格納する。学習オプションは
 *          MONOLIS_OPT_VAE_TRAIN_OPTS を共用する。
 * @ingroup opt
 */
typedef struct {
  /** Fortran monolis_opt_cvae_t への不透明ハンドル */
  void* net;
  /** 入力次元数 */
  int   D;
  /** 条件次元数 */
  int   C;
  /** 最終エンコーダ隠れ層次元 (互換用) */
  int   H;
  /** 潜在次元数 */
  int   Z;
} MONOLIS_OPT_CVAE;

/**
 * @brief CVAE モデルを初期化する (隠れ層 1 層 H)
 * @param[out] cvae CVAE モデル構造体
 * @param[in] D 入力次元数
 * @param[in] C 条件次元数
 * @param[in] H 隠れ層次元数
 * @param[in] Z 潜在次元数
 * @ingroup opt
 */
void monolis_opt_cvae_init(
  MONOLIS_OPT_CVAE* cvae,
  int               D,
  int               C,
  int               H,
  int               Z);

/**
 * @brief 任意層数の CVAE モデルを初期化する
 * @param[out] cvae CVAE モデル構造体
 * @param[in] D 入力次元数
 * @param[in] C 条件次元数
 * @param[in] H_enc エンコーダ隠れ層次元の配列 (Le)
 * @param[in] Le エンコーダ隠れ層数
 * @param[in] Z 潜在次元数
 * @param[in] H_dec デコーダ隠れ層次元の配列 (Ld)
 * @param[in] Ld デコーダ隠れ層数
 * @ingroup opt
 */
void monolis_opt_cvae_init_layers(
  MONOLIS_OPT_CVAE* cvae,
  int               D,
  int               C,
  int*              H_enc,
  int               Le,
  int               Z,
  int*              H_dec,
  int               Ld);

/**
 * @brief CVAE モデルが保持するメモリを解放する
 * @param[in,out] cvae CVAE モデル構造体
 * @ingroup opt
 */
void monolis_opt_cvae_finalize(
  MONOLIS_OPT_CVAE* cvae);

/**
 * @brief 1 ミニバッチに対する順伝播・逆伝播・Adam 更新を実行
 * @param[in,out] cvae CVAE モデル構造体
 * @param[in] X 入力ミニバッチ (D x B、列優先)
 * @param[in] Cnd 条件ミニバッチ (C x B、列優先)
 * @param[in] B バッチサイズ
 * @param[in] lr 学習率
 * @param[in] beta KL 重み
 * @param[out] loss バッチ平均損失
 * @param[out] recon バッチ平均再構成誤差
 * @param[out] kl バッチ平均 KL 損失
 * @ingroup opt
 */
void monolis_opt_cvae_train_step(
  MONOLIS_OPT_CVAE* cvae,
  float*            X,
  float*            Cnd,
  int               B,
  float             lr,
  float             beta,
  float*            loss,
  float*            recon,
  float*            kl);

/**
 * @brief ミニバッチ学習ループ (KL ウォームアップと EarlyStopping 付き)
 * @param[in,out] cvae CVAE モデル構造体
 * @param[in] X 学習データ (D x N、列優先)
 * @param[in] Cnd 条件データ (C x N、列優先)
 * @param[in] N 学習データ数
 * @param[in] opts 学習オプション
 * @ingroup opt
 */
void monolis_opt_cvae_fit(
  MONOLIS_OPT_CVAE*           cvae,
  float*                      X,
  float*                      Cnd,
  int                         N,
  MONOLIS_OPT_VAE_TRAIN_OPTS* opts);

/**
 * @brief 決定論的 (z = mu) な順伝播による再構成
 * @param[in] cvae CVAE モデル構造体
 * @param[in] X 入力データ (D x B、列優先)
 * @param[in] Cnd 条件 (C x B、列優先)
 * @param[in] B バッチサイズ
 * @param[out] xhat 再構成出力 (D x B、列優先)
 * @ingroup opt
 */
void monolis_opt_cvae_reconstruct(
  MONOLIS_OPT_CVAE* cvae,
  float*            X,
  float*            Cnd,
  int               B,
  float*            xhat);

/**
 * @brief 入力をエンコードし、再パラメータ化サンプル z を返す
 * @param[in] cvae CVAE モデル構造体
 * @param[in] X 入力データ (D x B、列優先)
 * @param[in] Cnd 条件 (C x B、列優先)
 * @param[in] B バッチサイズ
 * @param[out] zlat サンプリングされた潜在ベクトル (Z x B、列優先)
 * @ingroup opt
 */
void monolis_opt_cvae_encode(
  MONOLIS_OPT_CVAE* cvae,
  float*            X,
  float*            Cnd,
  int               B,
  float*            zlat);

/**
 * @brief 潜在ベクトルと条件をデコードして出力 (sigmoid) を得る
 * @param[in] cvae CVAE モデル構造体
 * @param[in] zlat 潜在ベクトル (Z x B、列優先)
 * @param[in] Cnd 条件 (C x B、列優先)
 * @param[in] B バッチサイズ
 * @param[out] xhat 出力 (D x B、列優先)
 * @ingroup opt
 */
void monolis_opt_cvae_decode(
  MONOLIS_OPT_CVAE* cvae,
  float*            zlat,
  float*            Cnd,
  int               B,
  float*            xhat);

/**
 * @brief 与えた条件ごとに標準正規分布から z をサンプリングしてデコードする
 * @param[in] cvae CVAE モデル構造体
 * @param[in] Cnd 条件 (C x n、列優先)
 * @param[in] n サンプル数
 * @param[out] xhat 出力 (D x n、列優先)
 * @ingroup opt
 */
void monolis_opt_cvae_sample_prior(
  MONOLIS_OPT_CVAE* cvae,
  float*            Cnd,
  int               n,
  float*            xhat);

/* Fortran bind(c) ブリッジ (内部利用) */
void* monolis_opt_cvae_init_c_main(
  int D, int C, int H, int Z);

void* monolis_opt_cvae_init_layers_c_main(
  int D, int C, int* H_enc, int Le, int Z, int* H_dec, int Ld);

void monolis_opt_cvae_finalize_c_main(
  void* handle);

void monolis_opt_cvae_train_step_c_main(
  void* handle, float* X, float* Cnd, int D, int C, int B, float lr, float beta,
  float* loss, float* recon, float* kl);

void monolis_opt_cvae_fit_c_main(
  void* handle, float* X, float* Cnd, int D, int C, int N,
  int batch_size, int epochs, float lr, float r_loss_factor,
  int kl_warmup_epochs, int early_stop_patience,
  int log_every_batches, int verbose);

void monolis_opt_cvae_reconstruct_c_main(
  void* handle, float* X, float* Cnd, int D, int C, int B, float* xhat);

void monolis_opt_cvae_encode_c_main(
  void* handle, float* X, float* Cnd, int D, int C, int B, int Z, float* zlat);

void monolis_opt_cvae_decode_c_main(
  void* handle, float* zlat, float* Cnd, int Z, int C, int B, int D, float* xhat);

void monolis_opt_cvae_sample_prior_c_main(
  void* handle, float* Cnd, int C, int n, int D, float* xhat);

#ifdef __cplusplus
}
#endif

#endif
