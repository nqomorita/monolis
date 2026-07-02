!> Conditional VAE C ラッパモジュール
!> @details C から Conditional VAE (mod_monolis_opt_cvae) を利用するための bind(c) ブリッジ。
!>          Fortran 派生型 monolis_opt_cvae_t をヒープに確保し、その C アドレス
!>          (type(c_ptr)) を不透明ハンドルとして C 側へ渡す。計算は既存の
!>          Fortran ルーチンをそのまま呼び出す (OpenACC のデータ管理も各ルーチン側に委譲)。
!>          機械学習に関わる実数は 32bit 浮動小数点 (c_float = kdouble_ml) で扱う。
module mod_monolis_cvae_wrapper
  use mod_monolis_utils
  use mod_monolis_def_opt
  use mod_monolis_opt_vae_util
  use mod_monolis_opt_cvae
  use iso_c_binding

  implicit none

contains

  !> CVAE モデルを初期化し不透明ハンドルを返す (隠れ層 1 層 H)
  function monolis_opt_cvae_init_c(D, C, H, Z) result(handle) &
    & bind(c, name = "monolis_opt_cvae_init_c_main")
    implicit none
    !> [in] 入力次元数
    integer(kint_c), value :: D
    !> [in] 条件次元数
    integer(kint_c), value :: C
    !> [in] 隠れ層次元数
    integer(kint_c), value :: H
    !> [in] 潜在次元数
    integer(kint_c), value :: Z
    !> [out] CVAE モデルの不透明ハンドル
    type(c_ptr) :: handle
    type(monolis_opt_cvae_t), pointer :: net

    allocate(net)
    call monolis_opt_cvae_init(net, D, C, H, Z)
    handle = c_loc(net)
  end function monolis_opt_cvae_init_c

  !> 任意層数の CVAE モデルを初期化し不透明ハンドルを返す
  function monolis_opt_cvae_init_layers_c(D, C, H_enc, Le, Z, H_dec, Ld) result(handle) &
    & bind(c, name = "monolis_opt_cvae_init_layers_c_main")
    implicit none
    !> [in] 入力次元数
    integer(kint_c), value :: D
    !> [in] 条件次元数
    integer(kint_c), value :: C
    !> [in] エンコーダ隠れ層数
    integer(kint_c), value :: Le
    !> [in] エンコーダ隠れ層次元の列 (Le)
    integer(kint_c) :: H_enc(Le)
    !> [in] 潜在次元数
    integer(kint_c), value :: Z
    !> [in] デコーダ隠れ層数
    integer(kint_c), value :: Ld
    !> [in] デコーダ隠れ層次元の列 (Ld)
    integer(kint_c) :: H_dec(Ld)
    !> [out] CVAE モデルの不透明ハンドル
    type(c_ptr) :: handle
    type(monolis_opt_cvae_t), pointer :: net

    allocate(net)
    call monolis_opt_cvae_init_layers(net, D, C, H_enc, Z, H_dec)
    handle = c_loc(net)
  end function monolis_opt_cvae_init_layers_c

  !> CVAE モデルが保持するメモリを解放しハンドルを破棄する
  subroutine monolis_opt_cvae_finalize_c(handle) &
    & bind(c, name = "monolis_opt_cvae_finalize_c_main")
    implicit none
    !> [in] CVAE モデルの不透明ハンドル
    type(c_ptr), value :: handle
    type(monolis_opt_cvae_t), pointer :: net

    call c_f_pointer(handle, net)
    call monolis_opt_cvae_finalize(net)
    deallocate(net)
  end subroutine monolis_opt_cvae_finalize_c

  !> 1 ミニバッチに対する順伝播・逆伝播・Adam 更新を実行
  subroutine monolis_opt_cvae_train_step_c(handle, X, Cnd, D, C, B, lr, beta, &
    & loss, recon, kl) &
    & bind(c, name = "monolis_opt_cvae_train_step_c_main")
    implicit none
    !> [in,out] CVAE モデルの不透明ハンドル
    type(c_ptr), value :: handle
    !> [in] 入力次元数
    integer(kint_c), value :: D
    !> [in] 条件次元数
    integer(kint_c), value :: C
    !> [in] バッチサイズ
    integer(kint_c), value :: B
    !> [in] 入力ミニバッチ (D x B)
    real(c_float) :: X(D, B)
    !> [in] 条件ミニバッチ (C x B)
    real(c_float) :: Cnd(C, B)
    !> [in] 学習率
    real(c_float), value :: lr
    !> [in] KL 重み (beta)
    real(c_float), value :: beta
    !> [out] バッチ平均損失
    real(c_float) :: loss
    !> [out] バッチ平均再構成誤差
    real(c_float) :: recon
    !> [out] バッチ平均 KL 損失
    real(c_float) :: kl
    type(monolis_opt_cvae_t), pointer :: net

    call c_f_pointer(handle, net)
    call monolis_opt_cvae_train_step(net, X, Cnd, lr, beta, loss, recon, kl)
  end subroutine monolis_opt_cvae_train_step_c

  !> ミニバッチ学習ループ (KL ウォームアップと EarlyStopping 付き)
  subroutine monolis_opt_cvae_fit_c(handle, X, Cnd, D, C, N, batch_size, epochs, lr, &
    & r_loss_factor, kl_warmup_epochs, early_stop_patience, log_every_batches, verbose) &
    & bind(c, name = "monolis_opt_cvae_fit_c_main")
    implicit none
    !> [in,out] CVAE モデルの不透明ハンドル
    type(c_ptr), value :: handle
    !> [in] 入力次元数
    integer(kint_c), value :: D
    !> [in] 条件次元数
    integer(kint_c), value :: C
    !> [in] 学習データ数
    integer(kint_c), value :: N
    !> [in] 学習データ (D x N)
    real(c_float) :: X(D, N)
    !> [in] 条件データ (C x N)
    real(c_float) :: Cnd(C, N)
    !> [in] ミニバッチサイズ
    integer(kint_c), value :: batch_size
    !> [in] エポック数
    integer(kint_c), value :: epochs
    !> [in] 学習率
    real(c_float), value :: lr
    !> [in] KL 損失の最終スケール (beta)
    real(c_float), value :: r_loss_factor
    !> [in] KL ウォームアップエポック数
    integer(kint_c), value :: kl_warmup_epochs
    !> [in] EarlyStopping の patience (負値で自動)
    integer(kint_c), value :: early_stop_patience
    !> [in] 何バッチごとに進捗をログするか
    integer(kint_c), value :: log_every_batches
    !> [in] エポックログを出すか (0 で抑制)
    integer(kint_c), value :: verbose
    type(monolis_opt_cvae_t), pointer :: net
    type(monolis_opt_vae_train_opts) :: opts

    call c_f_pointer(handle, net)
    opts%batch_size          = batch_size
    opts%epochs              = epochs
    opts%lr                  = lr
    opts%r_loss_factor       = r_loss_factor
    opts%kl_warmup_epochs    = kl_warmup_epochs
    opts%early_stop_patience = early_stop_patience
    opts%log_every_batches   = log_every_batches
    opts%verbose             = (verbose /= 0)
    call monolis_opt_cvae_fit(net, X, Cnd, opts)
  end subroutine monolis_opt_cvae_fit_c

  !> 決定論的 (z = mu) な順伝播による再構成
  subroutine monolis_opt_cvae_reconstruct_c(handle, X, Cnd, D, C, B, xhat) &
    & bind(c, name = "monolis_opt_cvae_reconstruct_c_main")
    implicit none
    !> [in] CVAE モデルの不透明ハンドル
    type(c_ptr), value :: handle
    !> [in] 入力次元数
    integer(kint_c), value :: D
    !> [in] 条件次元数
    integer(kint_c), value :: C
    !> [in] バッチサイズ
    integer(kint_c), value :: B
    !> [in] 入力データ (D x B)
    real(c_float) :: X(D, B)
    !> [in] 条件 (C x B)
    real(c_float) :: Cnd(C, B)
    !> [out] 再構成出力 (D x B)
    real(c_float) :: xhat(D, B)
    type(monolis_opt_cvae_t), pointer :: net

    call c_f_pointer(handle, net)
    call monolis_opt_cvae_reconstruct(net, X, Cnd, xhat)
  end subroutine monolis_opt_cvae_reconstruct_c

  !> 入力をエンコードし、再パラメータ化サンプル z を返す
  subroutine monolis_opt_cvae_encode_c(handle, X, Cnd, D, C, B, Z, zlat) &
    & bind(c, name = "monolis_opt_cvae_encode_c_main")
    implicit none
    !> [in] CVAE モデルの不透明ハンドル
    type(c_ptr), value :: handle
    !> [in] 入力次元数
    integer(kint_c), value :: D
    !> [in] 条件次元数
    integer(kint_c), value :: C
    !> [in] バッチサイズ
    integer(kint_c), value :: B
    !> [in] 潜在次元数
    integer(kint_c), value :: Z
    !> [in] 入力データ (D x B)
    real(c_float) :: X(D, B)
    !> [in] 条件 (C x B)
    real(c_float) :: Cnd(C, B)
    !> [out] サンプリングされた潜在ベクトル (Z x B)
    real(c_float) :: zlat(Z, B)
    type(monolis_opt_cvae_t), pointer :: net

    call c_f_pointer(handle, net)
    call monolis_opt_cvae_encode(net, X, Cnd, zlat)
  end subroutine monolis_opt_cvae_encode_c

  !> 潜在ベクトルと条件をデコードして出力 (sigmoid) を得る
  subroutine monolis_opt_cvae_decode_c(handle, zlat, Cnd, Z, C, B, D, xhat) &
    & bind(c, name = "monolis_opt_cvae_decode_c_main")
    implicit none
    !> [in] CVAE モデルの不透明ハンドル
    type(c_ptr), value :: handle
    !> [in] 潜在次元数
    integer(kint_c), value :: Z
    !> [in] 条件次元数
    integer(kint_c), value :: C
    !> [in] バッチサイズ
    integer(kint_c), value :: B
    !> [in] 入力次元数
    integer(kint_c), value :: D
    !> [in] 潜在ベクトル (Z x B)
    real(c_float) :: zlat(Z, B)
    !> [in] 条件 (C x B)
    real(c_float) :: Cnd(C, B)
    !> [out] 出力 (D x B)
    real(c_float) :: xhat(D, B)
    type(monolis_opt_cvae_t), pointer :: net

    call c_f_pointer(handle, net)
    call monolis_opt_cvae_decode(net, zlat, Cnd, xhat)
  end subroutine monolis_opt_cvae_decode_c

  !> 与えた条件ごとに標準正規分布から z をサンプリングしてデコードする
  subroutine monolis_opt_cvae_sample_prior_c(handle, Cnd, C, n, D, xhat) &
    & bind(c, name = "monolis_opt_cvae_sample_prior_c_main")
    implicit none
    !> [in] CVAE モデルの不透明ハンドル
    type(c_ptr), value :: handle
    !> [in] 条件次元数
    integer(kint_c), value :: C
    !> [in] サンプル数
    integer(kint_c), value :: n
    !> [in] 入力次元数
    integer(kint_c), value :: D
    !> [in] 条件 (C x n)
    real(c_float) :: Cnd(C, n)
    !> [out] 出力 (D x n)
    real(c_float) :: xhat(D, n)
    type(monolis_opt_cvae_t), pointer :: net

    call c_f_pointer(handle, net)
    call monolis_opt_cvae_sample_prior(net, Cnd, xhat)
  end subroutine monolis_opt_cvae_sample_prior_c
end module mod_monolis_cvae_wrapper
