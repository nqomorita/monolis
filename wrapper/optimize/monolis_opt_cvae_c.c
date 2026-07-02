#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "monolis_utils.h"
#include "monolis_opt_cvae_c.h"

void monolis_opt_cvae_init(
  MONOLIS_OPT_CVAE* cvae,
  int               D,
  int               C,
  int               H,
  int               Z)
{
  cvae->net = monolis_opt_cvae_init_c_main(D, C, H, Z);
  cvae->D = D;
  cvae->C = C;
  cvae->H = H;
  cvae->Z = Z;
}

void monolis_opt_cvae_init_layers(
  MONOLIS_OPT_CVAE* cvae,
  int               D,
  int               C,
  int*              H_enc,
  int               Le,
  int               Z,
  int*              H_dec,
  int               Ld)
{
  cvae->net = monolis_opt_cvae_init_layers_c_main(D, C, H_enc, Le, Z, H_dec, Ld);
  cvae->D = D;
  cvae->C = C;
  cvae->H = H_enc[Le - 1];
  cvae->Z = Z;
}

void monolis_opt_cvae_finalize(
  MONOLIS_OPT_CVAE* cvae)
{
  monolis_opt_cvae_finalize_c_main(cvae->net);
  cvae->net = NULL;
  cvae->D = 0;
  cvae->C = 0;
  cvae->H = 0;
  cvae->Z = 0;
}

void monolis_opt_cvae_train_step(
  MONOLIS_OPT_CVAE* cvae,
  float*            X,
  float*            Cnd,
  int               B,
  float             lr,
  float             beta,
  float*            loss,
  float*            recon,
  float*            kl)
{
  monolis_opt_cvae_train_step_c_main(
    cvae->net, X, Cnd, cvae->D, cvae->C, B, lr, beta, loss, recon, kl);
}

void monolis_opt_cvae_fit(
  MONOLIS_OPT_CVAE*           cvae,
  float*                      X,
  float*                      Cnd,
  int                         N,
  MONOLIS_OPT_VAE_TRAIN_OPTS* opts)
{
  monolis_opt_cvae_fit_c_main(
    cvae->net, X, Cnd, cvae->D, cvae->C, N,
    opts->batch_size, opts->epochs, opts->lr, opts->r_loss_factor,
    opts->kl_warmup_epochs, opts->early_stop_patience,
    opts->log_every_batches, opts->verbose);
}

void monolis_opt_cvae_reconstruct(
  MONOLIS_OPT_CVAE* cvae,
  float*            X,
  float*            Cnd,
  int               B,
  float*            xhat)
{
  monolis_opt_cvae_reconstruct_c_main(
    cvae->net, X, Cnd, cvae->D, cvae->C, B, xhat);
}

void monolis_opt_cvae_encode(
  MONOLIS_OPT_CVAE* cvae,
  float*            X,
  float*            Cnd,
  int               B,
  float*            zlat)
{
  monolis_opt_cvae_encode_c_main(
    cvae->net, X, Cnd, cvae->D, cvae->C, B, cvae->Z, zlat);
}

void monolis_opt_cvae_decode(
  MONOLIS_OPT_CVAE* cvae,
  float*            zlat,
  float*            Cnd,
  int               B,
  float*            xhat)
{
  monolis_opt_cvae_decode_c_main(
    cvae->net, zlat, Cnd, cvae->Z, cvae->C, B, cvae->D, xhat);
}

void monolis_opt_cvae_sample_prior(
  MONOLIS_OPT_CVAE* cvae,
  float*            Cnd,
  int               n,
  float*            xhat)
{
  monolis_opt_cvae_sample_prior_c_main(
    cvae->net, Cnd, cvae->C, n, cvae->D, xhat);
}
