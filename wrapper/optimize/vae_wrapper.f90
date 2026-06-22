!> VAE C ラッパモジュール
!> @details C から VAE (mod_monolis_opt_vae) を利用するための bind(c) ブリッジ。
!>          Fortran 派生型 monolis_opt_vae_t をヒープに確保し、その C アドレス
!>          (type(c_ptr)) を不透明ハンドルとして C 側へ渡す。計算は既存の
!>          Fortran ルーチンをそのまま呼び出す (OpenACC のデータ管理も各ルーチン側に委譲)。
!>          機械学習に関わる実数は 32bit 浮動小数点 (c_float = kdouble_ml) で扱う。
module mod_monolis_vae_wrapper
  use mod_monolis_utils
  use mod_monolis_def_opt
  use mod_monolis_opt_vae_util
  use mod_monolis_opt_vae
  use iso_c_binding

  implicit none

contains

  !> VAE モデルを初期化し不透明ハンドルを返す (隠れ層 1 層 H)
  function monolis_opt_vae_init_c(D, H, Z) result(handle) &
    & bind(c, name = "monolis_opt_vae_init_c_main")
    implicit none
    !> [in] 入力次元数
    integer(kint_c), value :: D
    !> [in] 隠れ層次元数
    integer(kint_c), value :: H
    !> [in] 潜在次元数
    integer(kint_c), value :: Z
    !> [out] VAE モデルの不透明ハンドル
    type(c_ptr) :: handle
    type(monolis_opt_vae_t), pointer :: net

    allocate(net)
    call monolis_opt_vae_init(net, D, H, Z)
    handle = c_loc(net)
  end function monolis_opt_vae_init_c

  !> 任意層数の VAE モデルを初期化し不透明ハンドルを返す
  function monolis_opt_vae_init_layers_c(D, H_enc, Le, Z, H_dec, Ld) result(handle) &
    & bind(c, name = "monolis_opt_vae_init_layers_c_main")
    implicit none
    !> [in] 入力次元数
    integer(kint_c), value :: D
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
    !> [out] VAE モデルの不透明ハンドル
    type(c_ptr) :: handle
    type(monolis_opt_vae_t), pointer :: net

    allocate(net)
    call monolis_opt_vae_init_layers(net, D, H_enc, Z, H_dec)
    handle = c_loc(net)
  end function monolis_opt_vae_init_layers_c

  !> VAE モデルが保持するメモリを解放しハンドルを破棄する
  subroutine monolis_opt_vae_finalize_c(handle) &
    & bind(c, name = "monolis_opt_vae_finalize_c_main")
    implicit none
    !> [in] VAE モデルの不透明ハンドル
    type(c_ptr), value :: handle
    type(monolis_opt_vae_t), pointer :: net

    call c_f_pointer(handle, net)
    call monolis_opt_vae_finalize(net)
    deallocate(net)
  end subroutine monolis_opt_vae_finalize_c

  !> 1 ミニバッチに対する順伝播・逆伝播・Adam 更新を実行
  subroutine monolis_opt_vae_train_step_c(handle, X, D, B, lr, beta, loss, recon, kl) &
    & bind(c, name = "monolis_opt_vae_train_step_c_main")
    implicit none
    !> [in,out] VAE モデルの不透明ハンドル
    type(c_ptr), value :: handle
    !> [in] 入力次元数
    integer(kint_c), value :: D
    !> [in] バッチサイズ
    integer(kint_c), value :: B
    !> [in] 入力ミニバッチ (D x B)
    real(c_float) :: X(D, B)
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
    type(monolis_opt_vae_t), pointer :: net

    call c_f_pointer(handle, net)
    call monolis_opt_vae_train_step(net, X, lr, beta, loss, recon, kl)
  end subroutine monolis_opt_vae_train_step_c

  !> ミニバッチ学習ループ (KL ウォームアップと EarlyStopping 付き)
  subroutine monolis_opt_vae_fit_c(handle, X, D, N, batch_size, epochs, lr, &
    & r_loss_factor, kl_warmup_epochs, early_stop_patience, log_every_batches, verbose) &
    & bind(c, name = "monolis_opt_vae_fit_c_main")
    implicit none
    !> [in,out] VAE モデルの不透明ハンドル
    type(c_ptr), value :: handle
    !> [in] 入力次元数
    integer(kint_c), value :: D
    !> [in] 学習データ数
    integer(kint_c), value :: N
    !> [in] 学習データ (D x N)
    real(c_float) :: X(D, N)
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
    type(monolis_opt_vae_t), pointer :: net
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
    call monolis_opt_vae_fit(net, X, opts)
  end subroutine monolis_opt_vae_fit_c

  !> 決定論的 (z = mu) な順伝播による再構成
  subroutine monolis_opt_vae_reconstruct_c(handle, X, D, B, xhat) &
    & bind(c, name = "monolis_opt_vae_reconstruct_c_main")
    implicit none
    !> [in] VAE モデルの不透明ハンドル
    type(c_ptr), value :: handle
    !> [in] 入力次元数
    integer(kint_c), value :: D
    !> [in] バッチサイズ
    integer(kint_c), value :: B
    !> [in] 入力データ (D x B)
    real(c_float) :: X(D, B)
    !> [out] 再構成出力 (D x B)
    real(c_float) :: xhat(D, B)
    type(monolis_opt_vae_t), pointer :: net

    call c_f_pointer(handle, net)
    call monolis_opt_vae_reconstruct(net, X, xhat)
  end subroutine monolis_opt_vae_reconstruct_c

  !> 入力をエンコードし、再パラメータ化サンプル z を返す
  subroutine monolis_opt_vae_encode_c(handle, X, D, B, Z, zlat) &
    & bind(c, name = "monolis_opt_vae_encode_c_main")
    implicit none
    !> [in] VAE モデルの不透明ハンドル
    type(c_ptr), value :: handle
    !> [in] 入力次元数
    integer(kint_c), value :: D
    !> [in] バッチサイズ
    integer(kint_c), value :: B
    !> [in] 潜在次元数
    integer(kint_c), value :: Z
    !> [in] 入力データ (D x B)
    real(c_float) :: X(D, B)
    !> [out] サンプリングされた潜在ベクトル (Z x B)
    real(c_float) :: zlat(Z, B)
    type(monolis_opt_vae_t), pointer :: net

    call c_f_pointer(handle, net)
    call monolis_opt_vae_encode(net, X, zlat)
  end subroutine monolis_opt_vae_encode_c

  !> 潜在ベクトルをデコードして出力 (sigmoid) を得る
  subroutine monolis_opt_vae_decode_c(handle, zlat, Z, B, D, xhat) &
    & bind(c, name = "monolis_opt_vae_decode_c_main")
    implicit none
    !> [in] VAE モデルの不透明ハンドル
    type(c_ptr), value :: handle
    !> [in] 潜在次元数
    integer(kint_c), value :: Z
    !> [in] バッチサイズ
    integer(kint_c), value :: B
    !> [in] 入力次元数
    integer(kint_c), value :: D
    !> [in] 潜在ベクトル (Z x B)
    real(c_float) :: zlat(Z, B)
    !> [out] 出力 (D x B)
    real(c_float) :: xhat(D, B)
    type(monolis_opt_vae_t), pointer :: net

    call c_f_pointer(handle, net)
    call monolis_opt_vae_decode(net, zlat, xhat)
  end subroutine monolis_opt_vae_decode_c

  !> 標準正規分布から n 個サンプリングしてデコードする
  subroutine monolis_opt_vae_sample_prior_c(handle, n, D, xhat) &
    & bind(c, name = "monolis_opt_vae_sample_prior_c_main")
    implicit none
    !> [in] VAE モデルの不透明ハンドル
    type(c_ptr), value :: handle
    !> [in] サンプル数
    integer(kint_c), value :: n
    !> [in] 入力次元数
    integer(kint_c), value :: D
    !> [out] 出力 (D x n)
    real(c_float) :: xhat(D, n)
    type(monolis_opt_vae_t), pointer :: net

    call c_f_pointer(handle, net)
    call monolis_opt_vae_sample_prior(net, n, xhat)
  end subroutine monolis_opt_vae_sample_prior_c

  !> 学習データを潜在空間にエンコードし SPX 交叉により n 個生成・デコード
  subroutine monolis_opt_vae_generate_spx_c(handle, Xtrain, D, Ntr, n, xhat) &
    & bind(c, name = "monolis_opt_vae_generate_spx_c_main")
    implicit none
    !> [in] VAE モデルの不透明ハンドル
    type(c_ptr), value :: handle
    !> [in] 入力次元数
    integer(kint_c), value :: D
    !> [in] 学習データ数
    integer(kint_c), value :: Ntr
    !> [in] 生成サンプル数
    integer(kint_c), value :: n
    !> [in] 学習データ (D x Ntr)
    real(c_float) :: Xtrain(D, Ntr)
    !> [out] 出力 (D x n)
    real(c_float) :: xhat(D, n)
    type(monolis_opt_vae_t), pointer :: net

    call c_f_pointer(handle, net)
    call monolis_opt_vae_generate_spx(net, Xtrain, n, xhat)
  end subroutine monolis_opt_vae_generate_spx_c
end module mod_monolis_vae_wrapper
