!> 条件付き変分オートエンコーダ (Conditional VAE) ライブラリ
!> @details 条件ベクトル c を入力・潜在変数に連結する CVAE。
!>          encoder : [x; c] (D+C) -> H_enc(1) -> ... -> H_enc(L) (ReLU)
!>                    -> head_mu / head_lv (linear, Z)
!>          z       = mu + exp(0.5*logvar) * eps,  eps ~ N(0,1)
!>          decoder : [z; c] (Z+C) -> H_dec(1) -> ... -> H_dec(M) (ReLU) -> D (sigmoid)
!>          loss    = mean_batch_dim( (x - xhat)^2 )
!>                  + beta * mean_batch( -0.5 * sum (1 + lv - mu^2 - exp(lv)) )
!>          層構造・Adam・順逆伝播などの共通処理は mod_monolis_opt_vae_util /
!>          mod_monolis_opt_adam に切り出している。
module mod_monolis_opt_cvae
  use mod_monolis_utils
  use mod_monolis_opt_adam
  use mod_monolis_opt_vae_util
  implicit none

  private
  public :: monolis_opt_cvae_t
  public :: monolis_opt_cvae_init
  public :: monolis_opt_cvae_init_layers
  public :: monolis_opt_cvae_finalize
  public :: monolis_opt_cvae_train_step
  public :: monolis_opt_cvae_fit
  public :: monolis_opt_cvae_reconstruct
  public :: monolis_opt_cvae_encode
  public :: monolis_opt_cvae_decode
  public :: monolis_opt_cvae_sample_prior

  !> @ingroup optimize
  !> Conditional VAE モデル (任意層数のエンコーダ・デコーダを保持)
  type :: monolis_opt_cvae_t
    !> [in] 入力次元数
    integer(kint) :: D = 0
    !> [in] 条件次元数
    integer(kint) :: C = 0
    !> [in] 互換用: 最終エンコーダ隠れ層次元
    integer(kint) :: H = 0
    !> [in] 潜在次元数
    integer(kint) :: Z = 0
    !> [in,out] Adam のステップ数
    integer(kint) :: t = 0
    !> [in,out] エンコーダ層列 ([x; c] (D+C) -> H_enc(1) -> ... -> H_enc(L))
    type(monolis_opt_vae_layer_t), allocatable :: enc(:)
    !> [in,out] mu 出力ヘッド (linear)
    type(monolis_opt_vae_layer_t) :: head_mu
    !> [in,out] logvar 出力ヘッド (linear)
    type(monolis_opt_vae_layer_t) :: head_lv
    !> [in,out] デコーダ層列 ([z; c] (Z+C) -> H_dec(1) -> ... -> H_dec(M) -> D)
    type(monolis_opt_vae_layer_t), allocatable :: dec(:)
  end type monolis_opt_cvae_t

contains

  !> @ingroup optimize
  !> 条件次元を連結して 2 次元配列を縦方向に結合する ([A; B])
  subroutine monolis_opt_cvae_concat_rows(A, Bm, AB)
    implicit none
    !> [in] 上段配列 (nA x ncol)
    real(kdouble), intent(in) :: A(:,:)
    !> [in] 下段配列 (nB x ncol)
    real(kdouble), intent(in) :: Bm(:,:)
    !> [out] 結合配列 ((nA+nB) x ncol)
    real(kdouble), allocatable, intent(out) :: AB(:,:)
    integer(kint) :: nA, nB, ncol

    nA = size(A, 1); nB = size(Bm, 1); ncol = size(A, 2)
    call monolis_alloc_R_2d(AB, nA + nB, ncol)
    AB(1:nA, :)            = A
    AB(nA+1:nA+nB, :)      = Bm
  end subroutine monolis_opt_cvae_concat_rows

  !> @ingroup optimize
  !> CVAE モデルを初期化する (隠れ層 1 層 H で構成)
  subroutine monolis_opt_cvae_init(net, D, C, H, Z)
    implicit none
    !> [out] 初期化対象の CVAE モデル
    type(monolis_opt_cvae_t), intent(out) :: net
    !> [in] 入力次元数
    integer(kint), intent(in) :: D
    !> [in] 条件次元数
    integer(kint), intent(in) :: C
    !> [in] 隠れ層次元数
    integer(kint), intent(in) :: H
    !> [in] 潜在次元数
    integer(kint), intent(in) :: Z

    call monolis_opt_cvae_init_layers(net, D, C, (/ H /), Z, (/ H /))
  end subroutine monolis_opt_cvae_init

  !> @ingroup optimize
  !> 任意層数の CVAE モデルを初期化する
  !> @details 重みは Glorot 一様分布、バイアスは 0 で初期化する。
  !>          エンコーダ入力は D+C、デコーダ入力は Z+C 次元。
  subroutine monolis_opt_cvae_init_layers(net, D, C, H_enc, Z, H_dec)
    implicit none
    !> [out] 初期化対象の CVAE モデル
    type(monolis_opt_cvae_t), intent(out) :: net
    !> [in] 入力次元数
    integer(kint), intent(in) :: D
    !> [in] 条件次元数
    integer(kint), intent(in) :: C
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
      call monolis_std_error_string("monolis_opt_cvae_init_layers: hidden size must be >= 1")
      call monolis_std_error_stop()
    endif

    net%D = D
    net%C = C
    net%H = H_enc(Le)
    net%Z = Z
    net%t = 0

    !> エンコーダ層列 (入力 D+C)
    allocate(net%enc(Le))
    in_dim = D + C
    do l = 1, Le
      call monolis_opt_vae_layer_init(net%enc(l), in_dim, H_enc(l), monolis_opt_vae_act_relu)
      in_dim = H_enc(l)
    end do

    !> mu / lv ヘッド (線形)
    call monolis_opt_vae_layer_init(net%head_mu, in_dim, Z, monolis_opt_vae_act_linear)
    call monolis_opt_vae_layer_init(net%head_lv, in_dim, Z, monolis_opt_vae_act_linear)

    !> デコーダ層列 (入力 Z+C、最終層のみ sigmoid)
    allocate(net%dec(Ld + 1))
    in_dim = Z + C
    do l = 1, Ld
      call monolis_opt_vae_layer_init(net%dec(l), in_dim, H_dec(l), monolis_opt_vae_act_relu)
      in_dim = H_dec(l)
    end do
    call monolis_opt_vae_layer_init(net%dec(Ld + 1), in_dim, D, monolis_opt_vae_act_sigmoid)
  end subroutine monolis_opt_cvae_init_layers

  !> @ingroup optimize
  !> CVAE モデルが保持するメモリを解放する
  subroutine monolis_opt_cvae_finalize(net)
    implicit none
    !> [in,out] 解放対象の CVAE モデル
    type(monolis_opt_cvae_t), intent(inout) :: net
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
    net%D = 0; net%C = 0; net%H = 0; net%Z = 0; net%t = 0
  end subroutine monolis_opt_cvae_finalize

  !> @ingroup optimize
  !> 1 ミニバッチに対する順伝播・逆伝播・Adam 更新を実行
  subroutine monolis_opt_cvae_train_step(net, X, Cnd, lr, beta, loss_avg, recon_avg, kl_avg)
    implicit none
    !> [in,out] CVAE モデル
    type(monolis_opt_cvae_t), intent(inout) :: net
    !> [in] 入力ミニバッチ (D x B)
    real(kdouble), intent(in) :: X(:,:)
    !> [in] 条件ミニバッチ (C x B)
    real(kdouble), intent(in) :: Cnd(:,:)
    !> [in] 学習率
    real(kdouble), intent(in) :: lr
    !> [in] KL 重み (beta)
    real(kdouble), intent(in) :: beta
    !> [out] バッチ平均損失
    real(kdouble), intent(out) :: loss_avg
    !> [out] バッチ平均再構成誤差
    real(kdouble), intent(out) :: recon_avg
    !> [out] バッチ平均 KL 損失
    real(kdouble), intent(out) :: kl_avg
    integer(kint) :: D, C, Z, B, j, l, Le, Ld
    type(monolis_opt_vae_cache_t), allocatable :: enc_cache(:), dec_cache(:)
    type(monolis_opt_vae_grad_t),  allocatable :: enc_grads(:), dec_grads(:)
    real(kdouble), allocatable :: enc_in(:,:), dec_in(:,:), h_last(:,:)
    real(kdouble), allocatable :: mu(:,:), logvar(:,:), eps(:,:), zlat(:,:), xhat(:,:)
    real(kdouble), allocatable :: dxhat(:,:), dzc_dec(:,:), dz_dec(:,:)
    real(kdouble), allocatable :: dmu(:,:), dlv(:,:)
    real(kdouble), allocatable :: gWmu(:,:), gWlv(:,:), gbmu(:), gblv(:)
    real(kdouble), allocatable :: dh1_mu(:,:), dh1_lv(:,:), dh_last(:,:), dX_dummy(:,:)
    real(kdouble) :: invB, invBD, kl_w

    D = net%D; C = net%C; Z = net%Z; B = size(X, 2)
    Le = size(net%enc); Ld = size(net%dec)
    invB  = 1.0d0 / real(B, kdouble)
    invBD = 1.0d0 / real(B*D, kdouble)
    kl_w  = beta * invB

    !> エンコーダ順伝播 (入力 [x; c])
    call monolis_opt_cvae_concat_rows(X, Cnd, enc_in)
    call monolis_opt_vae_forward_stack(net%enc, enc_in, enc_cache)
    h_last = monolis_opt_vae_top_post(enc_cache)

    !> mu / logvar (線形ヘッド: out_dim x B)
    call monolis_alloc_R_2d(mu,     Z, B)
    call monolis_alloc_R_2d(logvar, Z, B)
    mu     = matmul(transpose(net%head_mu%W), h_last)
    logvar = matmul(transpose(net%head_lv%W), h_last)
    do j = 1, B
      mu(:,j)     = mu(:,j)     + net%head_mu%b
      logvar(:,j) = logvar(:,j) + net%head_lv%b
    end do
    where(logvar >  10.0d0) logvar =  10.0d0
    where(logvar < -10.0d0) logvar = -10.0d0

    !> 再パラメータ化
    call monolis_alloc_R_2d(eps,  Z, B)
    call monolis_alloc_R_2d(zlat, Z, B)
    call monolis_opt_vae_fill_randn(eps)
    zlat = mu + exp(0.5d0*logvar) * eps

    !> デコーダ順伝播 (入力 [z; c])
    call monolis_opt_cvae_concat_rows(zlat, Cnd, dec_in)
    call monolis_opt_vae_forward_stack(net%dec, dec_in, dec_cache)
    xhat = monolis_opt_vae_top_post(dec_cache)

    !> 損失
    recon_avg = sum((X - xhat)**2) * invBD
    kl_avg    = -0.5d0 * sum(1.0d0 + logvar - mu*mu - exp(logvar)) * invB
    loss_avg  = recon_avg + beta * kl_avg

    !> デコーダ逆伝播: 出力 post=xhat に対する勾配 (sigmoid は層内で処理)
    call monolis_alloc_R_2d(dxhat, D, B)
    dxhat = 2.0d0 * (xhat - X) * invBD
    call monolis_opt_vae_backward_stack(net%dec, dec_in, dec_cache, dxhat, dec_grads, dzc_dec)

    !> デコーダ入力勾配のうち潜在変数 z に対応する部分のみ使用 (条件 c は入力)
    call monolis_alloc_R_2d(dz_dec, Z, B)
    dz_dec = dzc_dec(1:Z, :)

    !> dmu / dlv
    call monolis_alloc_R_2d(dmu, Z, B)
    call monolis_alloc_R_2d(dlv, Z, B)
    dmu = dz_dec + kl_w * mu
    dlv = dz_dec * eps * 0.5d0 * exp(0.5d0*logvar) + kl_w * 0.5d0*(exp(logvar) - 1.0d0)

    !> ヘッド (線形) の勾配と h_last への逆伝播
    call monolis_alloc_R_2d(gWmu, net%head_mu%in_dim, Z)
    call monolis_alloc_R_1d(gbmu, Z)
    call monolis_alloc_R_2d(gWlv, net%head_lv%in_dim, Z)
    call monolis_alloc_R_1d(gblv, Z)
    gWmu = matmul(h_last, transpose(dmu))
    gbmu = sum(dmu, dim=2)
    gWlv = matmul(h_last, transpose(dlv))
    gblv = sum(dlv, dim=2)
    call monolis_alloc_R_2d(dh1_mu, net%head_mu%in_dim, B)
    call monolis_alloc_R_2d(dh1_lv, net%head_lv%in_dim, B)
    dh1_mu = matmul(net%head_mu%W, dmu)
    dh1_lv = matmul(net%head_lv%W, dlv)
    call monolis_alloc_R_2d(dh_last, net%head_mu%in_dim, B)
    dh_last = dh1_mu + dh1_lv

    !> エンコーダ逆伝播
    call monolis_opt_vae_backward_stack(net%enc, enc_in, enc_cache, dh_last, enc_grads, dX_dummy)

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
  end subroutine monolis_opt_cvae_train_step

  !> @ingroup optimize
  !> ミニバッチ学習ループ (KL ウォームアップと EarlyStopping 付き)
  subroutine monolis_opt_cvae_fit(net, X, Cnd, opts)
    implicit none
    !> [in,out] 学習対象の CVAE モデル
    type(monolis_opt_cvae_t), intent(inout) :: net
    !> [in] 学習データ (D x N)
    real(kdouble), intent(in) :: X(:,:)
    !> [in] 条件データ (C x N)
    real(kdouble), intent(in) :: Cnd(:,:)
    !> [in] 学習オプション
    type(monolis_opt_vae_train_opts), intent(in) :: opts
    integer(kint) :: D, C, N, e, b, k, no_improve, nbatches, patience
    integer(kint), allocatable :: perm(:)
    real(kdouble), allocatable :: Xb(:,:), Cb(:,:)
    real(kdouble) :: loss, recon, kl, sum_loss, sum_recon, sum_kl
    real(kdouble) :: best_loss, epoch_loss, beta

    D = size(X, 1)
    C = size(Cnd, 1)
    N = size(X, 2)
    if(D /= net%D)then
      call monolis_std_error_string("monolis_opt_cvae_fit: data dim mismatch")
      call monolis_std_error_stop()
    endif
    if(C /= net%C)then
      call monolis_std_error_string("monolis_opt_cvae_fit: condition dim mismatch")
      call monolis_std_error_stop()
    endif

    nbatches = N / opts%batch_size
    if(nbatches < 1)then
      call monolis_std_error_string("monolis_opt_cvae_fit: not enough samples for one batch")
      call monolis_std_error_stop()
    endif
    patience = opts%early_stop_patience
    if(patience < 0) patience = max(1, opts%epochs / 10)

    call monolis_alloc_I_1d(perm, N)
    call monolis_alloc_R_2d(Xb, D, opts%batch_size)
    call monolis_alloc_R_2d(Cb, C, opts%batch_size)
    do k = 1, N
      perm(k) = k
    end do

    best_loss  = huge(1.0d0)
    no_improve = 0

    do e = 1, opts%epochs
      if(opts%kl_warmup_epochs <= 0)then
        beta = opts%r_loss_factor
      else
        beta = opts%r_loss_factor * &
               min(1.0d0, real(e, kdouble) / real(opts%kl_warmup_epochs, kdouble))
      endif

      call monolis_opt_vae_shuffle(perm)
      sum_loss = 0.0d0; sum_recon = 0.0d0; sum_kl = 0.0d0
      do b = 1, nbatches
        do k = 1, opts%batch_size
          Xb(:, k) = X(:,   perm((b-1)*opts%batch_size + k))
          Cb(:, k) = Cnd(:, perm((b-1)*opts%batch_size + k))
        end do
        call monolis_opt_cvae_train_step(net, Xb, Cb, opts%lr, beta, loss, recon, kl)
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
  end subroutine monolis_opt_cvae_fit

  !> @ingroup optimize
  !> 共通の順伝播: [x; c] -> (enc + heads) -> mu/logvar/h_last
  subroutine monolis_opt_cvae_encode_mu_lv(net, X, Cnd, mu, logvar, enc_cache, h_last)
    implicit none
    !> [in] CVAE モデル
    type(monolis_opt_cvae_t), intent(in) :: net
    !> [in] 入力 (D x B)
    real(kdouble), intent(in) :: X(:,:)
    !> [in] 条件 (C x B)
    real(kdouble), intent(in) :: Cnd(:,:)
    !> [out] mu (Z x B)
    real(kdouble), allocatable, intent(out) :: mu(:,:)
    !> [out] logvar (Z x B、[-10,10] でクリップ)
    real(kdouble), allocatable, intent(out) :: logvar(:,:)
    !> [out] エンコーダキャッシュ (不要なら呼び出し後に解放)
    type(monolis_opt_vae_cache_t), allocatable, intent(out) :: enc_cache(:)
    !> [out] 最終エンコーダ隠れ層出力
    real(kdouble), allocatable, intent(out) :: h_last(:,:)
    real(kdouble), allocatable :: enc_in(:,:)
    integer(kint) :: j, B

    B = size(X, 2)
    call monolis_opt_cvae_concat_rows(X, Cnd, enc_in)
    call monolis_opt_vae_forward_stack(net%enc, enc_in, enc_cache)
    h_last = monolis_opt_vae_top_post(enc_cache)
    call monolis_alloc_R_2d(mu,     net%Z, B)
    call monolis_alloc_R_2d(logvar, net%Z, B)
    mu     = matmul(transpose(net%head_mu%W), h_last)
    logvar = matmul(transpose(net%head_lv%W), h_last)
    do j = 1, B
      mu(:,j)     = mu(:,j)     + net%head_mu%b
      logvar(:,j) = logvar(:,j) + net%head_lv%b
    end do
    where(logvar >  10.0d0) logvar =  10.0d0
    where(logvar < -10.0d0) logvar = -10.0d0
  end subroutine monolis_opt_cvae_encode_mu_lv

  !> @ingroup optimize
  !> デコード本体 (戻り値は allocatable、入力 [z; c])
  subroutine monolis_opt_cvae_decode_alloc(net, zlat, Cnd, xhat)
    implicit none
    !> [in] CVAE モデル
    type(monolis_opt_cvae_t), intent(in) :: net
    !> [in] 潜在ベクトル (Z x B)
    real(kdouble), intent(in) :: zlat(:,:)
    !> [in] 条件 (C x B)
    real(kdouble), intent(in) :: Cnd(:,:)
    !> [out] 出力 (D x B)
    real(kdouble), allocatable, intent(out) :: xhat(:,:)
    real(kdouble), allocatable :: dec_in(:,:)
    type(monolis_opt_vae_cache_t), allocatable :: dec_cache(:)

    call monolis_opt_cvae_concat_rows(zlat, Cnd, dec_in)
    call monolis_opt_vae_forward_stack(net%dec, dec_in, dec_cache)
    xhat = monolis_opt_vae_top_post(dec_cache)
    call monolis_opt_vae_cache_free(dec_cache)
  end subroutine monolis_opt_cvae_decode_alloc

  !> @ingroup optimize
  !> 決定論的 (z = mu) な順伝播による再構成
  subroutine monolis_opt_cvae_reconstruct(net, X, Cnd, xhat)
    implicit none
    !> [in] CVAE モデル
    type(monolis_opt_cvae_t), intent(in) :: net
    !> [in] 入力データ (D x B)
    real(kdouble), intent(in) :: X(:,:)
    !> [in] 条件 (C x B)
    real(kdouble), intent(in) :: Cnd(:,:)
    !> [out] 再構成出力 (D x B)
    real(kdouble), intent(out) :: xhat(:,:)
    real(kdouble), allocatable :: mu(:,:), logvar(:,:), h_last(:,:), xhat_loc(:,:)
    type(monolis_opt_vae_cache_t), allocatable :: enc_cache(:)

    call monolis_opt_cvae_encode_mu_lv(net, X, Cnd, mu, logvar, enc_cache, h_last)
    call monolis_opt_cvae_decode_alloc(net, mu, Cnd, xhat_loc)
    xhat = xhat_loc
    call monolis_opt_vae_cache_free(enc_cache)
  end subroutine monolis_opt_cvae_reconstruct

  !> @ingroup optimize
  !> 入力をエンコードし、再パラメータ化サンプル z を返す
  subroutine monolis_opt_cvae_encode(net, X, Cnd, zlat)
    implicit none
    !> [in] CVAE モデル
    type(monolis_opt_cvae_t), intent(in) :: net
    !> [in] 入力データ (D x B)
    real(kdouble), intent(in) :: X(:,:)
    !> [in] 条件 (C x B)
    real(kdouble), intent(in) :: Cnd(:,:)
    !> [out] サンプリングされた潜在ベクトル (Z x B)
    real(kdouble), intent(out) :: zlat(:,:)
    real(kdouble), allocatable :: mu(:,:), logvar(:,:), eps(:,:), h_last(:,:)
    type(monolis_opt_vae_cache_t), allocatable :: enc_cache(:)

    call monolis_opt_cvae_encode_mu_lv(net, X, Cnd, mu, logvar, enc_cache, h_last)
    call monolis_alloc_R_2d(eps, net%Z, size(X, 2))
    call monolis_opt_vae_fill_randn(eps)
    zlat = mu + exp(0.5d0*logvar) * eps
    call monolis_opt_vae_cache_free(enc_cache)
  end subroutine monolis_opt_cvae_encode

  !> @ingroup optimize
  !> 潜在ベクトルと条件をデコードして出力 (sigmoid) を得る
  subroutine monolis_opt_cvae_decode(net, zlat, Cnd, xhat)
    implicit none
    !> [in] CVAE モデル
    type(monolis_opt_cvae_t), intent(in) :: net
    !> [in] 潜在ベクトル (Z x B)
    real(kdouble), intent(in) :: zlat(:,:)
    !> [in] 条件 (C x B)
    real(kdouble), intent(in) :: Cnd(:,:)
    !> [out] 出力 (D x B)
    real(kdouble), intent(out) :: xhat(:,:)
    real(kdouble), allocatable :: xhat_loc(:,:)

    call monolis_opt_cvae_decode_alloc(net, zlat, Cnd, xhat_loc)
    xhat = xhat_loc
  end subroutine monolis_opt_cvae_decode

  !> @ingroup optimize
  !> 与えた条件ごとに標準正規分布から z をサンプリングしてデコードする
  subroutine monolis_opt_cvae_sample_prior(net, Cnd, xhat)
    implicit none
    !> [in] CVAE モデル
    type(monolis_opt_cvae_t), intent(in) :: net
    !> [in] 条件 (C x n)
    real(kdouble), intent(in) :: Cnd(:,:)
    !> [out] 出力 (D x n)
    real(kdouble), intent(out) :: xhat(:,:)
    real(kdouble), allocatable :: zlat(:,:)
    integer(kint) :: n

    n = size(Cnd, 2)
    call monolis_alloc_R_2d(zlat, net%Z, n)
    call monolis_opt_vae_fill_randn(zlat)
    call monolis_opt_cvae_decode(net, zlat, Cnd, xhat)
  end subroutine monolis_opt_cvae_sample_prior

end module mod_monolis_opt_cvae
