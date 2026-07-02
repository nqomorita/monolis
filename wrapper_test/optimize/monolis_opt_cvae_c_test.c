#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "monolis_solver.h"
#include "monolis_opt_cvae_c_test.h"

/* 配列の全要素が [0,1] に収まるか判定 */
static int monolis_opt_cvae_all_in_unit(const float* a, int n)
{
  int i;
  for (i = 0; i < n; ++i) {
    if (!(a[i] >= 0.0f && a[i] <= 1.0f)) return 0;
  }
  return 1;
}

/* スカラ値が有限か判定 */
static int monolis_opt_cvae_is_finite(float v)
{
  return (v == v) && (fabsf(v) < HUGE_VALF);
}

void monolis_opt_cvae_c_test()
{
  MONOLIS_OPT_CVAE cvae;
  MONOLIS_OPT_VAE_TRAIN_OPTS opts;
  const int D = 4, C = 3, H = 8, Z = 2, N = 32;
  float X[4 * 32];
  float Cnd[3 * 32];
  float xhat[4 * 32];
  float zlat[2 * 32];
  float samples[4 * 4];
  float Csmp[3 * 4];
  float loss0, loss, recon, kl;
  int H_enc[2] = {6, 5};
  int H_dec[2] = {5, 6};
  int i, j;

  monolis_std_log_string("monolis_opt_cvae_c_test");

  /* ダミーデータ (列優先) */
  for (j = 0; j < N; ++j) {
    for (i = 0; i < D; ++i) {
      X[i + j * D] = 0.5f + 0.4f * sinf((float)((i + 1) + (j + 1)));
    }
    for (i = 0; i < C; ++i) {
      Cnd[i + j * C] = 0.5f + 0.3f * cosf((float)((i + 1) + (j + 1)));
    }
  }
  for (j = 0; j < 4; ++j) {
    for (i = 0; i < C; ++i) {
      Csmp[i + j * C] = 0.5f;
    }
  }

  /* 初期化 (隠れ層 1 層) */
  monolis_opt_cvae_init(&cvae, D, C, H, Z);
  monolis_test_check_eq_I1("cvae_c_test dim D", cvae.D, D);
  monolis_test_check_eq_I1("cvae_c_test dim C", cvae.C, C);
  monolis_test_check_eq_I1("cvae_c_test dim H", cvae.H, H);
  monolis_test_check_eq_I1("cvae_c_test dim Z", cvae.Z, Z);

  /* 1 ステップ学習で損失が有限値であること */
  monolis_opt_cvae_train_step(&cvae, X, Cnd, N, 1.0e-3f, 1.0f, &loss0, &recon, &kl);
  monolis_test_check_eq_I1("cvae_c_test loss finite",
    monolis_opt_cvae_is_finite(loss0), 1);

  /* 短い fit 呼び出し */
  opts.batch_size          = 8;
  opts.epochs              = 5;
  opts.lr                  = 1.0e-2f;
  opts.r_loss_factor       = 1.0f;
  opts.kl_warmup_epochs    = 2;
  opts.early_stop_patience = 100;
  opts.log_every_batches   = 0;
  opts.verbose             = 0;
  monolis_opt_cvae_fit(&cvae, X, Cnd, N, &opts);

  /* 再構成と値域 [0,1] */
  monolis_opt_cvae_reconstruct(&cvae, X, Cnd, N, xhat);
  monolis_test_check_eq_I1("cvae_c_test reconstruct in [0,1]",
    monolis_opt_cvae_all_in_unit(xhat, D * N), 1);

  /* エンコード / デコード */
  monolis_opt_cvae_encode(&cvae, X, Cnd, N, zlat);
  monolis_opt_cvae_decode(&cvae, zlat, Cnd, N, xhat);
  monolis_test_check_eq_I1("cvae_c_test decode in [0,1]",
    monolis_opt_cvae_all_in_unit(xhat, D * N), 1);

  /* 条件付き事前分布からのサンプル */
  monolis_opt_cvae_sample_prior(&cvae, Csmp, 4, samples);
  monolis_test_check_eq_I1("cvae_c_test prior in [0,1]",
    monolis_opt_cvae_all_in_unit(samples, D * 4), 1);

  monolis_opt_cvae_finalize(&cvae);
  monolis_test_check_eq_I1("cvae_c_test finalize", cvae.D, 0);

  /* 任意層数の初期化も学習できること */
  monolis_opt_cvae_init_layers(&cvae, D, C, H_enc, 2, Z, H_dec, 2);
  monolis_test_check_eq_I1("cvae_c_test layers dim H", cvae.H, H_enc[1]);
  monolis_opt_cvae_train_step(&cvae, X, Cnd, N, 1.0e-3f, 1.0f, &loss, &recon, &kl);
  monolis_test_check_eq_I1("cvae_c_test layers loss finite",
    monolis_opt_cvae_is_finite(loss), 1);
  monolis_opt_cvae_finalize(&cvae);
}
