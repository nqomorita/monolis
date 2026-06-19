!> 変分オートエンコーダ (VAE) ライブラリ
!> @details 任意層数 MLP のエンコーダ・デコーダをサポート。
!>          encoder : D -> H_enc(1) -> ... -> H_enc(L) (ReLU)
!>                    -> head_mu / head_lv (linear, Z)
!>          z       = mu + exp(0.5*logvar) * eps,  eps ~ N(0,1)
!>          decoder : Z -> H_dec(1) -> ... -> H_dec(M) (ReLU) -> D (sigmoid)
!>          loss    = mean_batch_dim( (x - xhat)^2 )
!>                  + beta * mean_batch( -0.5 * sum (1 + lv - mu^2 - exp(lv)) )
!>          層構造・Adam・順逆伝播などの共通処理は mod_monolis_opt_vae_util /
!>          mod_monolis_opt_adam に切り出している。
!>          機械学習に関わる実数は 32bit 浮動小数点 (kdouble_ml) で計算する。
module mod_monolis_opt_vae
  use mod_monolis_utils
  use mod_monolis_def_opt
  use mod_monolis_opt_adam
  use mod_monolis_opt_vae_util
  implicit none

  private
  public :: monolis_opt_vae_t
  public :: monolis_opt_vae_init
  public :: monolis_opt_vae_init_layers
  public :: monolis_opt_vae_finalize
  public :: monolis_opt_vae_train_step
  public :: monolis_opt_vae_fit
  public :: monolis_opt_vae_reconstruct
  public :: monolis_opt_vae_encode
  public :: monolis_opt_vae_decode
  public :: monolis_opt_vae_sample_prior
  public :: monolis_opt_vae_generate_spx

  !> @ingroup optimize
  !> VAE モデル (任意層数のエンコーダ・デコーダを保持)
  type :: monolis_opt_vae_t
    !> [in] 入力次元数
    integer(kint) :: D = 0
    !> [in] 互換用: 最終エンコーダ隠れ層次元 (旧 API との互換)
    integer(kint) :: H = 0
    !> [in] 潜在次元数
    integer(kint) :: Z = 0
    !> [in,out] Adam のステップ数
    integer(kint) :: t = 0
    !> [in,out] エンコーダ層列 (D -> H_enc(1) -> ... -> H_enc(L))
    type(monolis_opt_vae_layer_t), allocatable :: enc(:)
    !> [in,out] mu 出力ヘッド (linear)
    type(monolis_opt_vae_layer_t) :: head_mu
    !> [in,out] logvar 出力ヘッド (linear)
    type(monolis_opt_vae_layer_t) :: head_lv
    !> [in,out] デコーダ層列 (Z -> H_dec(1) -> ... -> H_dec(M) -> D)
    type(monolis_opt_vae_layer_t), allocatable :: dec(:)
  end type monolis_opt_vae_t

contains

  !> @ingroup optimize
  !> VAE モデルを初期化する (旧 API。隠れ層 1 層 H で構成)
  subroutine monolis_opt_vae_init(net, D, H, Z)
    implicit none
    !> [out] 初期化対象の VAE モデル
    type(monolis_opt_vae_t), intent(out) :: net
    !> [in] 入力次元数
    integer(kint), intent(in) :: D
    !> [in] 隠れ層次元数
    integer(kint), intent(in) :: H
    !> [in] 潜在次元数
    integer(kint), intent(in) :: Z

    call monolis_opt_vae_init_layers(net, D, (/ H /), Z, (/ H /))
  end subroutine monolis_opt_vae_init

  !> @ingroup optimize
  !> 任意層数の VAE モデルを初期化する
  !> @details 重みは Glorot 一様分布、バイアスは 0 で初期化する。
  !>          乱数消費順 (旧 API 互換): enc(1)..enc(L), head_mu, head_lv, dec(1)..dec(M+1)
  subroutine monolis_opt_vae_init_layers(net, D, H_enc, Z, H_dec)
    implicit none
    !> [out] 初期化対象の VAE モデル
    type(monolis_opt_vae_t), intent(out) :: net
    !> [in] 入力次元数
    integer(kint), intent(in) :: D
    !> [in] エンコーダ隠れ層次元の列 (size>=1、ReLU 層)
    integer(kint), intent(in) :: H_enc(:)
    !> [in] 潜在次元数
    integer(kint), intent(in) :: Z
    !> [in] デコーダ隠れ層次元の列 (size>=1、ReLU 層)
    integer(kint), intent(in) :: H_dec(:)
    integer(kint) :: l, in_dim, Le, Ld

    Le = size(H_enc)
    Ld = size(H_dec)
    if(Le < 1 .or. Ld < 1)then
      call monolis_std_error_string("monolis_opt_vae_init_layers: hidden size must be >= 1")
      call monolis_std_error_stop()
    endif

    net%D = D
    net%H = H_enc(Le)
    net%Z = Z
    net%t = 0

    !> エンコーダ層列
    allocate(net%enc(Le))
    in_dim = D
    do l = 1, Le
      call monolis_opt_vae_layer_init(net%enc(l), in_dim, H_enc(l), monolis_opt_vae_act_relu)
      in_dim = H_enc(l)
    end do

    !> mu / lv ヘッド (線形)
    call monolis_opt_vae_layer_init(net%head_mu, in_dim, Z, monolis_opt_vae_act_linear)
    call monolis_opt_vae_layer_init(net%head_lv, in_dim, Z, monolis_opt_vae_act_linear)

    !> デコーダ層列 (最終層のみ sigmoid)
    allocate(net%dec(Ld + 1))
    in_dim = Z
    do l = 1, Ld
      call monolis_opt_vae_layer_init(net%dec(l), in_dim, H_dec(l), monolis_opt_vae_act_relu)
      in_dim = H_dec(l)
    end do
    call monolis_opt_vae_layer_init(net%dec(Ld + 1), in_dim, D, monolis_opt_vae_act_sigmoid)
  end subroutine monolis_opt_vae_init_layers

  !> @ingroup optimize
  !> VAE モデルが保持するメモリを解放する
  subroutine monolis_opt_vae_finalize(net)
    implicit none
    !> [in,out] 解放対象の VAE モデル
    type(monolis_opt_vae_t), intent(inout) :: net
    integer(kint) :: l

    if(allocated(net%enc))then
      do l = 1, size(net%enc)
        call monolis_opt_vae_layer_free(net%enc(l))
      end do
      deallocate(net%enc)
    endif
    call monolis_opt_vae_layer_free(net%head_mu)
    call monolis_opt_vae_layer_free(net%head_lv)
    if(allocated(net%dec))then
      do l = 1, size(net%dec)
        call monolis_opt_vae_layer_free(net%dec(l))
      end do
      deallocate(net%dec)
    endif
    net%D = 0; net%H = 0; net%Z = 0; net%t = 0
  end subroutine monolis_opt_vae_finalize

  !> @ingroup optimize
  !> 1 ミニバッチに対する順伝播・逆伝播・Adam 更新を実行
  subroutine monolis_opt_vae_train_step(net, X, lr, beta, loss_avg, recon_avg, kl_avg)
    implicit none
    !> [in,out] VAE モデル
    type(monolis_opt_vae_t), intent(inout) :: net
    !> [in] 入力ミニバッチ (D x B)
    real(kdouble_ml), intent(in) :: X(:,:)
    !> [in] 学習率
    real(kdouble_ml), intent(in) :: lr
    !> [in] KL 重み (beta)
    real(kdouble_ml), intent(in) :: beta
    !> [out] バッチ平均損失
    real(kdouble_ml), intent(out) :: loss_avg
    !> [out] バッチ平均再構成誤差
    real(kdouble_ml), intent(out) :: recon_avg
    !> [out] バッチ平均 KL 損失
    real(kdouble_ml), intent(out) :: kl_avg
    integer(kint) :: D, Z, B, j, l, Le, Ld
    type(monolis_opt_vae_cache_t), allocatable :: enc_cache(:), dec_cache(:)
    type(monolis_opt_vae_grad_t),  allocatable :: enc_grads(:), dec_grads(:)
    real(kdouble_ml), allocatable :: h_last(:,:)
    real(kdouble_ml), allocatable :: mu(:,:), logvar(:,:), eps(:,:), zlat(:,:), xhat(:,:)
    real(kdouble_ml), allocatable :: dxhat(:,:), dz_dec(:,:)
    real(kdouble_ml), allocatable :: dmu(:,:), dlv(:,:)
    real(kdouble_ml), allocatable :: gWmu(:,:), gWlv(:,:), gbmu(:), gblv(:)
    real(kdouble_ml), allocatable :: dh1_mu(:,:), dh1_lv(:,:), dh_last(:,:), dX_dummy(:,:)
    real(kdouble_ml) :: invB, invBD, kl_w

    D = net%D; Z = net%Z; B = size(X, 2)
    Le = size(net%enc); Ld = size(net%dec)
    invB  = 1.0_kdouble_ml / real(B, kdouble_ml)
    invBD = 1.0_kdouble_ml / real(B*D, kdouble_ml)
    kl_w  = beta * invB

    !> エンコーダ順伝播
    call monolis_opt_vae_forward_stack(net%enc, X, enc_cache)
    h_last = monolis_opt_vae_top_post(enc_cache)

    !> mu / logvar (線形ヘッド: out_dim x B)
    allocate(mu(Z, B),     source = 0.0_kdouble_ml)
    allocate(logvar(Z, B), source = 0.0_kdouble_ml)
    mu     = matmul(transpose(net%head_mu%W), h_last)
    logvar = matmul(transpose(net%head_lv%W), h_last)
    do j = 1, B
      mu(:,j)     = mu(:,j)     + net%head_mu%b
      logvar(:,j) = logvar(:,j) + net%head_lv%b
    end do
    where(logvar >  10.0_kdouble_ml) logvar =  10.0_kdouble_ml
    where(logvar < -10.0_kdouble_ml) logvar = -10.0_kdouble_ml

    !> 再パラメータ化
    allocate(eps(Z, B),  source = 0.0_kdouble_ml)
    allocate(zlat(Z, B), source = 0.0_kdouble_ml)
    call monolis_opt_vae_fill_randn(eps)
    zlat = mu + exp(0.5_kdouble_ml*logvar) * eps

    !> デコーダ順伝播
    call monolis_opt_vae_forward_stack(net%dec, zlat, dec_cache)
    xhat = monolis_opt_vae_top_post(dec_cache)

    !> 損失
    recon_avg = sum((X - xhat)**2) * invBD
    kl_avg    = -0.5_kdouble_ml * sum(1.0_kdouble_ml + logvar - mu*mu - exp(logvar)) * invB
    loss_avg  = recon_avg + beta * kl_avg

    !> デコーダ逆伝播: 出力 post=xhat に対する勾配 (sigmoid は層内で処理)
    allocate(dxhat(D, B), source = 0.0_kdouble_ml)
    dxhat = 2.0_kdouble_ml * (xhat - X) * invBD
    call monolis_opt_vae_backward_stack(net%dec, zlat, dec_cache, dxhat, dec_grads, dz_dec)

    !> dmu / dlv
    allocate(dmu(Z, B), source = 0.0_kdouble_ml)
    allocate(dlv(Z, B), source = 0.0_kdouble_ml)
    dmu = dz_dec + kl_w * mu
    dlv = dz_dec * eps * 0.5_kdouble_ml * exp(0.5_kdouble_ml*logvar) + kl_w * 0.5_kdouble_ml*(exp(logvar) - 1.0_kdouble_ml)

    !> ヘッド (線形) の勾配と h_last への逆伝播
    allocate(gWmu(net%head_mu%in_dim, Z), source = 0.0_kdouble_ml)
    allocate(gbmu(Z),                     source = 0.0_kdouble_ml)
    allocate(gWlv(net%head_lv%in_dim, Z), source = 0.0_kdouble_ml)
    allocate(gblv(Z),                     source = 0.0_kdouble_ml)
    gWmu = matmul(h_last, transpose(dmu))
    gbmu = sum(dmu, dim=2)
    gWlv = matmul(h_last, transpose(dlv))
    gblv = sum(dlv, dim=2)
    allocate(dh1_mu(net%head_mu%in_dim, B), source = 0.0_kdouble_ml)
    allocate(dh1_lv(net%head_lv%in_dim, B), source = 0.0_kdouble_ml)
    dh1_mu = matmul(net%head_mu%W, dmu)
    dh1_lv = matmul(net%head_lv%W, dlv)
    allocate(dh_last(net%head_mu%in_dim, B), source = 0.0_kdouble_ml)
    dh_last = dh1_mu + dh1_lv

    !> エンコーダ逆伝播
    call monolis_opt_vae_backward_stack(net%enc, X, enc_cache, dh_last, enc_grads, dX_dummy)

    !> Adam 更新
    net%t = net%t + 1
    do l = 1, Le
      call monolis_opt_adam_apply(net%enc(l)%W, net%enc(l)%b, &
        enc_grads(l)%gW, enc_grads(l)%gb, net%enc(l)%a, net%t, lr)
    end do
    call monolis_opt_adam_apply(net%head_mu%W, net%head_mu%b, gWmu, gbmu, &
      net%head_mu%a, net%t, lr)
    call monolis_opt_adam_apply(net%head_lv%W, net%head_lv%b, gWlv, gblv, &
      net%head_lv%a, net%t, lr)
    do l = 1, Ld
      call monolis_opt_adam_apply(net%dec(l)%W, net%dec(l)%b, &
        dec_grads(l)%gW, dec_grads(l)%gb, net%dec(l)%a, net%t, lr)
    end do

    !> 一時メモリ解放
    call monolis_opt_vae_cache_free(enc_cache)
    call monolis_opt_vae_cache_free(dec_cache)
    call monolis_opt_vae_grads_free(enc_grads)
    call monolis_opt_vae_grads_free(dec_grads)
  end subroutine monolis_opt_vae_train_step

  !> @ingroup optimize
  !> ミニバッチ学習ループ (KL ウォームアップと EarlyStopping 付き)
  subroutine monolis_opt_vae_fit(net, X, opts)
    implicit none
    !> [in,out] 学習対象の VAE モデル
    type(monolis_opt_vae_t), intent(inout) :: net
    !> [in] 学習データ (D x N)
    real(kdouble_ml), intent(in) :: X(:,:)
    !> [in] 学習オプション
    type(monolis_opt_vae_train_opts), intent(in) :: opts
    integer(kint) :: D, N, e, b, k, no_improve, nbatches, patience
    integer(kint), allocatable :: perm(:)
    real(kdouble_ml), allocatable :: Xb(:,:)
    real(kdouble_ml) :: loss, recon, kl, sum_loss, sum_recon, sum_kl
    real(kdouble_ml) :: best_loss, epoch_loss, beta

    D = size(X, 1)
    N = size(X, 2)
    if(D /= net%D)then
      call monolis_std_error_string("monolis_opt_vae_fit: data dim mismatch")
      call monolis_std_error_stop()
    endif

    nbatches = N / opts%batch_size
    if(nbatches < 1)then
      call monolis_std_error_string("monolis_opt_vae_fit: not enough samples for one batch")
      call monolis_std_error_stop()
    endif
    patience = opts%early_stop_patience
    if(patience < 0) patience = max(1, opts%epochs / 10)

    call monolis_alloc_I_1d(perm, N)
    allocate(Xb(D, opts%batch_size), source = 0.0_kdouble_ml)
    do k = 1, N
      perm(k) = k
    end do

    best_loss  = huge(0.0_kdouble_ml)
    no_improve = 0

    do e = 1, opts%epochs
      if(opts%kl_warmup_epochs <= 0)then
        beta = opts%r_loss_factor
      else
        beta = opts%r_loss_factor * &
               min(1.0_kdouble_ml, real(e, kdouble_ml) / real(opts%kl_warmup_epochs, kdouble_ml))
      endif

      call monolis_opt_vae_shuffle(perm)
      sum_loss = 0.0_kdouble_ml; sum_recon = 0.0_kdouble_ml; sum_kl = 0.0_kdouble_ml
      do b = 1, nbatches
        do k = 1, opts%batch_size
          Xb(:, k) = X(:, perm((b-1)*opts%batch_size + k))
        end do
        call monolis_opt_vae_train_step(net, Xb, opts%lr, beta, loss, recon, kl)
        sum_loss  = sum_loss  + loss
        sum_recon = sum_recon + recon
        sum_kl    = sum_kl    + kl
        if(opts%verbose .and. opts%log_every_batches > 0)then
          if(mod(b, opts%log_every_batches) == 0)then
            write(*,'(a,i0,a,i0,a,i0,a,es12.4,a,es12.4,a,es12.4)') &
              ' epoch ', e, '  batch ', b, '/', nbatches, &
              '  loss=', loss, '  recon=', recon, '  kl=', kl
          endif
        endif
      end do
      epoch_loss = sum_loss / nbatches
      if(opts%verbose)then
        write(*,'(a,i0,a,f6.3,a,es12.4,a,es12.4,a,es12.4)') &
          '== epoch ', e, '  beta=', beta, &
          '  avg loss=', epoch_loss, &
          '  recon=', sum_recon/nbatches, '  kl=', sum_kl/nbatches
      endif

      if(e <= opts%kl_warmup_epochs)then
        best_loss  = epoch_loss
        no_improve = 0
      else if(epoch_loss < best_loss)then
        best_loss  = epoch_loss
        no_improve = 0
      else
        no_improve = no_improve + 1
        if(no_improve >= patience)then
          if(opts%verbose)then
            write(*,'(a,i0,a,i0,a)') &
              'EarlyStopping: no improvement for ', no_improve, &
              ' epochs (patience=', patience, ')'
          endif
          exit
        endif
      endif
    end do
  end subroutine monolis_opt_vae_fit

  !> @ingroup optimize
  !> 共通の順伝播: X -> (enc + heads) -> mu/logvar/h_last
  subroutine monolis_opt_vae_encode_mu_lv(net, X, mu, logvar, enc_cache, h_last)
    implicit none
    !> [in] VAE モデル
    type(monolis_opt_vae_t), intent(in) :: net
    !> [in] 入力 (D x B)
    real(kdouble_ml), intent(in) :: X(:,:)
    !> [out] mu (Z x B)
    real(kdouble_ml), allocatable, intent(out) :: mu(:,:)
    !> [out] logvar (Z x B、[-10,10] でクリップ)
    real(kdouble_ml), allocatable, intent(out) :: logvar(:,:)
    !> [out] エンコーダキャッシュ (不要なら呼び出し後に解放)
    type(monolis_opt_vae_cache_t), allocatable, intent(out) :: enc_cache(:)
    !> [out] 最終エンコーダ隠れ層出力
    real(kdouble_ml), allocatable, intent(out) :: h_last(:,:)
    integer(kint) :: j, B

    B = size(X, 2)
    call monolis_opt_vae_forward_stack(net%enc, X, enc_cache)
    h_last = monolis_opt_vae_top_post(enc_cache)
    allocate(mu(net%Z, B),     source = 0.0_kdouble_ml)
    allocate(logvar(net%Z, B), source = 0.0_kdouble_ml)
    mu     = matmul(transpose(net%head_mu%W), h_last)
    logvar = matmul(transpose(net%head_lv%W), h_last)
    do j = 1, B
      mu(:,j)     = mu(:,j)     + net%head_mu%b
      logvar(:,j) = logvar(:,j) + net%head_lv%b
    end do
    where(logvar >  10.0_kdouble_ml) logvar =  10.0_kdouble_ml
    where(logvar < -10.0_kdouble_ml) logvar = -10.0_kdouble_ml
  end subroutine monolis_opt_vae_encode_mu_lv

  !> @ingroup optimize
  !> 決定論的 (z = mu) な順伝播による再構成
  subroutine monolis_opt_vae_reconstruct(net, X, xhat)
    implicit none
    !> [in] VAE モデル
    type(monolis_opt_vae_t), intent(in) :: net
    !> [in] 入力データ (D x B)
    real(kdouble_ml), intent(in) :: X(:,:)
    !> [out] 再構成出力 (D x B)
    real(kdouble_ml), intent(out) :: xhat(:,:)
    real(kdouble_ml), allocatable :: mu(:,:), logvar(:,:), h_last(:,:), xhat_loc(:,:)
    type(monolis_opt_vae_cache_t), allocatable :: enc_cache(:)

    call monolis_opt_vae_encode_mu_lv(net, X, mu, logvar, enc_cache, h_last)
    call monolis_opt_vae_decode_alloc(net, mu, xhat_loc)
    xhat = xhat_loc
    call monolis_opt_vae_cache_free(enc_cache)
  end subroutine monolis_opt_vae_reconstruct

  !> @ingroup optimize
  !> 入力をエンコードし、再パラメータ化サンプル z を返す
  subroutine monolis_opt_vae_encode(net, X, zlat)
    implicit none
    !> [in] VAE モデル
    type(monolis_opt_vae_t), intent(in) :: net
    !> [in] 入力データ (D x B)
    real(kdouble_ml), intent(in) :: X(:,:)
    !> [out] サンプリングされた潜在ベクトル (Z x B)
    real(kdouble_ml), intent(out) :: zlat(:,:)
    real(kdouble_ml), allocatable :: mu(:,:), logvar(:,:), eps(:,:), h_last(:,:)
    type(monolis_opt_vae_cache_t), allocatable :: enc_cache(:)

    call monolis_opt_vae_encode_mu_lv(net, X, mu, logvar, enc_cache, h_last)
    allocate(eps(net%Z, size(X, 2)), source = 0.0_kdouble_ml)
    call monolis_opt_vae_fill_randn(eps)
    zlat = mu + exp(0.5_kdouble_ml*logvar) * eps
    call monolis_opt_vae_cache_free(enc_cache)
  end subroutine monolis_opt_vae_encode

  !> @ingroup optimize
  !> デコード本体 (戻り値は allocatable)
  subroutine monolis_opt_vae_decode_alloc(net, zlat, xhat)
    implicit none
    !> [in] VAE モデル
    type(monolis_opt_vae_t), intent(in) :: net
    !> [in] 潜在ベクトル (Z x B)
    real(kdouble_ml), intent(in) :: zlat(:,:)
    !> [out] 出力 (D x B)
    real(kdouble_ml), allocatable, intent(out) :: xhat(:,:)
    type(monolis_opt_vae_cache_t), allocatable :: dec_cache(:)

    call monolis_opt_vae_forward_stack(net%dec, zlat, dec_cache)
    xhat = monolis_opt_vae_top_post(dec_cache)
    call monolis_opt_vae_cache_free(dec_cache)
  end subroutine monolis_opt_vae_decode_alloc

  !> @ingroup optimize
  !> 潜在ベクトルをデコードして出力 (sigmoid) を得る
  subroutine monolis_opt_vae_decode(net, zlat, xhat)
    implicit none
    !> [in] VAE モデル
    type(monolis_opt_vae_t), intent(in) :: net
    !> [in] 潜在ベクトル (Z x B)
    real(kdouble_ml), intent(in) :: zlat(:,:)
    !> [out] 出力 (D x B)
    real(kdouble_ml), intent(out) :: xhat(:,:)
    real(kdouble_ml), allocatable :: xhat_loc(:,:)

    call monolis_opt_vae_decode_alloc(net, zlat, xhat_loc)
    xhat = xhat_loc
  end subroutine monolis_opt_vae_decode

  !> @ingroup optimize
  !> 標準正規分布から n 個サンプリングしてデコードする
  subroutine monolis_opt_vae_sample_prior(net, n, xhat)
    implicit none
    !> [in] VAE モデル
    type(monolis_opt_vae_t), intent(in) :: net
    !> [in] サンプル数
    integer(kint), intent(in) :: n
    !> [out] 出力 (D x n)
    real(kdouble_ml), intent(out) :: xhat(:,:)
    real(kdouble_ml), allocatable :: zlat(:,:)

    allocate(zlat(net%Z, n), source = 0.0_kdouble_ml)
    call monolis_opt_vae_fill_randn(zlat)
    call monolis_opt_vae_decode(net, zlat, xhat)
  end subroutine monolis_opt_vae_sample_prior

  !> @ingroup optimize
  !> 学習データを潜在空間にエンコードしたうえで SPX 交叉により n 個生成・デコード
  subroutine monolis_opt_vae_generate_spx(net, Xtrain, n, xhat)
    implicit none
    !> [in] VAE モデル
    type(monolis_opt_vae_t), intent(in) :: net
    !> [in] 学習データ (D x Ntr)
    real(kdouble_ml), intent(in) :: Xtrain(:,:)
    !> [in] 生成サンプル数
    integer(kint), intent(in) :: n
    !> [out] 出力 (D x n)
    real(kdouble_ml), intent(out) :: xhat(:,:)
    integer(kint) :: Ntr, Np, i, k, idx
    real(kdouble_ml) :: r
    real(kdouble_ml), allocatable :: zall(:,:), parents(:,:), zchild(:,:)

    Ntr = size(Xtrain, 2)
    Np  = net%Z + 1
    allocate(zall(net%Z, Ntr),   source = 0.0_kdouble_ml)
    call monolis_opt_vae_encode(net, Xtrain, zall)
    allocate(parents(net%Z, Np), source = 0.0_kdouble_ml)
    allocate(zchild(net%Z, n),   source = 0.0_kdouble_ml)
    do i = 1, n
      do k = 1, Np
        call random_number(r)
        idx = 1 + int(r * Ntr, kint)
        if(idx > Ntr) idx = Ntr
        parents(:, k) = zall(:, idx)
      end do
      call monolis_opt_vae_spx_child(parents, net%Z, zchild(:, i))
    end do
    call monolis_opt_vae_decode(net, zchild, xhat)
  end subroutine monolis_opt_vae_generate_spx

end module mod_monolis_opt_vae
