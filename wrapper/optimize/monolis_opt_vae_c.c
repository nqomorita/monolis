#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "monolis_utils.h"
#include "monolis_opt_vae_c.h"

void monolis_opt_vae_init(
  MONOLIS_OPT_VAE* vae,
  int              D,
  int              H,
  int              Z)
{
  vae->net = monolis_opt_vae_init_c_main(D, H, Z);
  vae->D = D;
  vae->H = H;
  vae->Z = Z;
}

void monolis_opt_vae_init_layers(
  MONOLIS_OPT_VAE* vae,
  int              D,
  int*             H_enc,
  int              Le,
  int              Z,
  int*             H_dec,
  int              Ld)
{
  vae->net = monolis_opt_vae_init_layers_c_main(D, H_enc, Le, Z, H_dec, Ld);
  vae->D = D;
  vae->H = H_enc[Le - 1];
  vae->Z = Z;
}

void monolis_opt_vae_finalize(
  MONOLIS_OPT_VAE* vae)
{
  monolis_opt_vae_finalize_c_main(vae->net);
  vae->net = NULL;
  vae->D = 0;
  vae->H = 0;
  vae->Z = 0;
}

void monolis_opt_vae_train_step(
  MONOLIS_OPT_VAE* vae,
  float*           X,
  int              B,
  float            lr,
  float            beta,
  float*           loss,
  float*           recon,
  float*           kl)
{
  monolis_opt_vae_train_step_c_main(
    vae->net, X, vae->D, B, lr, beta, loss, recon, kl);
}

void monolis_opt_vae_fit(
  MONOLIS_OPT_VAE*            vae,
  float*                     X,
  int                        N,
  MONOLIS_OPT_VAE_TRAIN_OPTS* opts)
{
  monolis_opt_vae_fit_c_main(
    vae->net, X, vae->D, N,
    opts->batch_size, opts->epochs, opts->lr, opts->r_loss_factor,
    opts->kl_warmup_epochs, opts->early_stop_patience,
    opts->log_every_batches, opts->verbose);
}

void monolis_opt_vae_reconstruct(
  MONOLIS_OPT_VAE* vae,
  float*           X,
  int              B,
  float*           xhat)
{
  monolis_opt_vae_reconstruct_c_main(vae->net, X, vae->D, B, xhat);
}

void monolis_opt_vae_encode(
  MONOLIS_OPT_VAE* vae,
  float*           X,
  int              B,
  float*           zlat)
{
  monolis_opt_vae_encode_c_main(vae->net, X, vae->D, B, vae->Z, zlat);
}

void monolis_opt_vae_decode(
  MONOLIS_OPT_VAE* vae,
  float*           zlat,
  int              B,
  float*           xhat)
{
  monolis_opt_vae_decode_c_main(vae->net, zlat, vae->Z, B, vae->D, xhat);
}

void monolis_opt_vae_sample_prior(
  MONOLIS_OPT_VAE* vae,
  int              n,
  float*           xhat)
{
  monolis_opt_vae_sample_prior_c_main(vae->net, n, vae->D, xhat);
}

void monolis_opt_vae_generate_spx(
  MONOLIS_OPT_VAE* vae,
  float*           Xtrain,
  int              Ntr,
  int              n,
  float*           xhat)
{
  monolis_opt_vae_generate_spx_c_main(vae->net, Xtrain, vae->D, Ntr, n, xhat);
}
