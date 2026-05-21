!> 変分オートエンコーダ (VAE) ライブラリ
!> @details モデル構造体・順伝播/逆伝播・Adam 更新・学習ループ・潜在空間
!>          サンプリング (標準正規および SPX 交叉) を 1 ファイルに集約。
!>          encoder : D -> H (ReLU) -> mu(Z), logvar(Z)
!>          z       = mu + exp(0.5*logvar) * eps,  eps ~ N(0,1)
!>          decoder : Z -> H (ReLU) -> D (sigmoid)
!>          loss    = mean_batch( sum (x - xhat)^2 )
!>                  + beta * mean_batch( -0.5 * sum (1 + lv - mu^2 - exp(lv)) )
module mod_monolis_opt_vae
  use mod_monolis_utils
  implicit none

  private
  public :: monolis_opt_vae_adam_state
  public :: monolis_opt_vae_t
  public :: monolis_opt_vae_train_opts
  public :: monolis_opt_vae_init
  public :: monolis_opt_vae_finalize
  public :: monolis_opt_vae_train_step
  public :: monolis_opt_vae_fit
  public :: monolis_opt_vae_reconstruct
  public :: monolis_opt_vae_encode
  public :: monolis_opt_vae_decode
  public :: monolis_opt_vae_sample_prior
  public :: monolis_opt_vae_spx_child
  public :: monolis_opt_vae_generate_spx

  !> @ingroup optimize
  !> Adam オプティマイザの状態 (1 つの重み行列 + バイアスに対応)
  type :: monolis_opt_vae_adam_state
    !> [in,out] 重み行列の 1 次モーメント
    real(kdouble), allocatable :: m(:,:)
    !> [in,out] 重み行列の 2 次モーメント
    real(kdouble), allocatable :: v(:,:)
    !> [in,out] バイアスの 1 次モーメント
    real(kdouble), allocatable :: mb(:)
    !> [in,out] バイアスの 2 次モーメント
    real(kdouble), allocatable :: vb(:)
  end type monolis_opt_vae_adam_state

  !> @ingroup optimize
  !> VAE モデル (重み・バイアス・Adam 状態を保持)
  type :: monolis_opt_vae_t
    !> [in] 入力次元数
    integer(kint) :: D = 0
    !> [in] 隠れ層次元数
    integer(kint) :: H = 0
    !> [in] 潜在次元数
    integer(kint) :: Z = 0
    !> [in,out] エンコーダ重み (D x H)
    real(kdouble), allocatable :: W1(:,:)
    !> [in,out] エンコーダバイアス (H)
    real(kdouble), allocatable :: b1(:)
    !> [in,out] mu 用重み (H x Z)
    real(kdouble), allocatable :: Wmu(:,:)
    !> [in,out] mu 用バイアス (Z)
    real(kdouble), allocatable :: bmu(:)
    !> [in,out] logvar 用重み (H x Z)
    real(kdouble), allocatable :: Wlv(:,:)
    !> [in,out] logvar 用バイアス (Z)
    real(kdouble), allocatable :: blv(:)
    !> [in,out] デコーダ第 1 層重み (Z x H)
    real(kdouble), allocatable :: W2(:,:)
    !> [in,out] デコーダ第 1 層バイアス (H)
    real(kdouble), allocatable :: b2(:)
    !> [in,out] デコーダ第 2 層重み (H x D)
    real(kdouble), allocatable :: W3(:,:)
    !> [in,out] デコーダ第 2 層バイアス (D)
    real(kdouble), allocatable :: b3(:)
    !> [in,out] 各重みに対応する Adam 状態
    type(monolis_opt_vae_adam_state) :: a_W1, a_Wmu, a_Wlv, a_W2, a_W3
    !> [in,out] Adam のステップ数
    integer(kint) :: t = 0
  end type monolis_opt_vae_t

  !> @ingroup optimize
  !> 学習ハイパーパラメータ
  type :: monolis_opt_vae_train_opts
    !> [in] ミニバッチサイズ
    integer(kint) :: batch_size = 128
    !> [in] エポック数
    integer(kint) :: epochs = 10
    !> [in] 学習率
    real(kdouble) :: lr = 1.0d-3
    !> [in] KL 損失の最終スケール (beta)
    real(kdouble) :: r_loss_factor = 1.0d0
    !> [in] KL ウォームアップエポック数 (0 -> r_loss_factor へ線形、0 で無効)
    integer(kint) :: kl_warmup_epochs = 0
    !> [in] EarlyStopping の patience。負値で max(1, epochs/10)
    integer(kint) :: early_stop_patience = -1
    !> [in] 何バッチごとに進捗をログするか (0 で抑制)
    integer(kint) :: log_every_batches = 0
    !> [in] エポックログを出すか
    logical :: verbose = .false.
  end type monolis_opt_vae_train_opts

contains

  !> @ingroup optimize
  !> VAE モデルを初期化する (Glorot 一様分布で重みを初期化、バイアスは 0)
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

    net%D = D
    net%H = H
    net%Z = Z
    call monolis_alloc_R_2d(net%W1,  D, H)
    call monolis_alloc_R_1d(net%b1,  H)
    call monolis_alloc_R_2d(net%Wmu, H, Z)
    call monolis_alloc_R_1d(net%bmu, Z)
    call monolis_alloc_R_2d(net%Wlv, H, Z)
    call monolis_alloc_R_1d(net%blv, Z)
    call monolis_alloc_R_2d(net%W2,  Z, H)
    call monolis_alloc_R_1d(net%b2,  H)
    call monolis_alloc_R_2d(net%W3,  H, D)
    call monolis_alloc_R_1d(net%b3,  D)

    call monolis_opt_vae_glorot_uniform(net%W1)
    call monolis_opt_vae_glorot_uniform(net%Wmu)
    call monolis_opt_vae_glorot_uniform(net%Wlv)
    call monolis_opt_vae_glorot_uniform(net%W2)
    call monolis_opt_vae_glorot_uniform(net%W3)

    call monolis_opt_vae_adam_init(net%a_W1,  D, H)
    call monolis_opt_vae_adam_init(net%a_Wmu, H, Z)
    call monolis_opt_vae_adam_init(net%a_Wlv, H, Z)
    call monolis_opt_vae_adam_init(net%a_W2,  Z, H)
    call monolis_opt_vae_adam_init(net%a_W3,  H, D)

    net%t = 0
  end subroutine monolis_opt_vae_init

  !> @ingroup optimize
  !> VAE モデルが保持するメモリを解放する
  subroutine monolis_opt_vae_finalize(net)
    implicit none
    !> [in,out] 解放対象の VAE モデル
    type(monolis_opt_vae_t), intent(inout) :: net

    if(allocated(net%W1))  deallocate(net%W1)
    if(allocated(net%b1))  deallocate(net%b1)
    if(allocated(net%Wmu)) deallocate(net%Wmu)
    if(allocated(net%bmu)) deallocate(net%bmu)
    if(allocated(net%Wlv)) deallocate(net%Wlv)
    if(allocated(net%blv)) deallocate(net%blv)
    if(allocated(net%W2))  deallocate(net%W2)
    if(allocated(net%b2))  deallocate(net%b2)
    if(allocated(net%W3))  deallocate(net%W3)
    if(allocated(net%b3))  deallocate(net%b3)
    call monolis_opt_vae_adam_free(net%a_W1)
    call monolis_opt_vae_adam_free(net%a_Wmu)
    call monolis_opt_vae_adam_free(net%a_Wlv)
    call monolis_opt_vae_adam_free(net%a_W2)
    call monolis_opt_vae_adam_free(net%a_W3)
    net%D = 0; net%H = 0; net%Z = 0; net%t = 0
  end subroutine monolis_opt_vae_finalize

  !> @ingroup optimize
  !> Adam 状態を初期化 (ゼロ埋め)
  subroutine monolis_opt_vae_adam_init(st, n1, n2)
    implicit none
    !> [out] 初期化対象の Adam 状態
    type(monolis_opt_vae_adam_state), intent(out) :: st
    !> [in] 重み行列の行数
    integer(kint), intent(in) :: n1
    !> [in] 重み行列の列数 (バイアス次元)
    integer(kint), intent(in) :: n2

    call monolis_alloc_R_2d(st%m,  n1, n2)
    call monolis_alloc_R_2d(st%v,  n1, n2)
    call monolis_alloc_R_1d(st%mb, n2)
    call monolis_alloc_R_1d(st%vb, n2)
  end subroutine monolis_opt_vae_adam_init

  !> @ingroup optimize
  !> Adam 状態のメモリを解放
  subroutine monolis_opt_vae_adam_free(st)
    implicit none
    !> [in,out] 解放対象の Adam 状態
    type(monolis_opt_vae_adam_state), intent(inout) :: st

    if(allocated(st%m))  deallocate(st%m)
    if(allocated(st%v))  deallocate(st%v)
    if(allocated(st%mb)) deallocate(st%mb)
    if(allocated(st%vb)) deallocate(st%vb)
  end subroutine monolis_opt_vae_adam_free

  !> @ingroup optimize
  !> Adam による 1 ステップ更新
  subroutine monolis_opt_vae_adam_apply(W, b, gW, gb, st, t, lr)
    implicit none
    !> [in,out] 重み行列
    real(kdouble), intent(inout) :: W(:,:)
    !> [in,out] バイアス
    real(kdouble), intent(inout) :: b(:)
    !> [in] 重み勾配
    real(kdouble), intent(in) :: gW(:,:)
    !> [in] バイアス勾配
    real(kdouble), intent(in) :: gb(:)
    !> [in,out] Adam 状態
    type(monolis_opt_vae_adam_state), intent(inout) :: st
    !> [in] グローバルステップ数 (1 始まり)
    integer(kint), intent(in) :: t
    !> [in] 学習率
    real(kdouble), intent(in) :: lr
    real(kdouble), parameter :: b1 = 0.9d0, b2 = 0.999d0, eps = 1.0d-7
    real(kdouble) :: bc1, bc2

    st%m  = b1*st%m  + (1.0d0 - b1)*gW
    st%v  = b2*st%v  + (1.0d0 - b2)*gW*gW
    st%mb = b1*st%mb + (1.0d0 - b1)*gb
    st%vb = b2*st%vb + (1.0d0 - b2)*gb*gb
    bc1 = 1.0d0 - b1**t
    bc2 = 1.0d0 - b2**t
    W = W - lr * (st%m / bc1) / (sqrt(st%v / bc2) + eps)
    b = b - lr * (st%mb/ bc1) / (sqrt(st%vb/ bc2) + eps)
  end subroutine monolis_opt_vae_adam_apply

  !> @ingroup optimize
  !> Glorot/Xavier 一様初期化 : U(-sqrt(6/(in+out)), +sqrt(6/(in+out)))
  subroutine monolis_opt_vae_glorot_uniform(W)
    implicit none
    !> [out] 初期化対象の重み行列 (fan_in x fan_out)
    real(kdouble), intent(out) :: W(:,:)
    real(kdouble) :: limit, r
    integer(kint) :: i, j, fan_in, fan_out

    fan_in  = size(W, 1)
    fan_out = size(W, 2)
    limit = sqrt(6.0d0 / real(fan_in + fan_out, kdouble))
    do j = 1, size(W, 2)
      do i = 1, size(W, 1)
        call random_number(r)
        W(i,j) = (2.0d0*r - 1.0d0) * limit
      end do
    end do
  end subroutine monolis_opt_vae_glorot_uniform

  !> @ingroup optimize
  !> Box-Muller による標準正規乱数で配列を埋める
  subroutine monolis_opt_vae_fill_randn(a)
    implicit none
    !> [out] 充填対象の 2 次元配列
    real(kdouble), intent(out) :: a(:,:)
    real(kdouble), parameter :: pi = 3.141592653589793d0
    real(kdouble) :: u1, u2
    integer(kint) :: i, j

    do j = 1, size(a, 2)
      do i = 1, size(a, 1)
        call random_number(u1)
        call random_number(u2)
        if(u1 < 1.0d-12) u1 = 1.0d-12
        a(i,j) = sqrt(-2.0d0*log(u1)) * cos(2.0d0*pi*u2)
      end do
    end do
  end subroutine monolis_opt_vae_fill_randn

  !> @ingroup optimize
  !> Fisher-Yates シャッフル
  subroutine monolis_opt_vae_shuffle(a)
    implicit none
    !> [in,out] シャッフル対象の整数配列
    integer(kint), intent(inout) :: a(:)
    integer(kint) :: i, j, tmp, n
    real(kdouble) :: r

    n = size(a)
    do i = n, 2, -1
      call random_number(r)
      j = 1 + int(r * i, kint)
      if(j > i) j = i
      tmp = a(i); a(i) = a(j); a(j) = tmp
    end do
  end subroutine monolis_opt_vae_shuffle

  !> @ingroup optimize
  !> 1 ミニバッチに対する順伝播・逆伝播・Adam 更新を実行
  subroutine monolis_opt_vae_train_step(net, X, lr, beta, loss_avg, recon_avg, kl_avg)
    implicit none
    !> [in,out] VAE モデル
    type(monolis_opt_vae_t), intent(inout) :: net
    !> [in] 入力ミニバッチ (D x B)
    real(kdouble), intent(in) :: X(:,:)
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
    integer(kint) :: D, H, Z, B, j
    real(kdouble), allocatable :: h1pre(:,:), h1(:,:)
    real(kdouble), allocatable :: mu(:,:), logvar(:,:), eps(:,:), zlat(:,:)
    real(kdouble), allocatable :: h2pre(:,:), h2(:,:), logits(:,:), xhat(:,:)
    real(kdouble), allocatable :: dlogits(:,:), dxhat(:,:), dh2(:,:), dh2pre(:,:)
    real(kdouble), allocatable :: dz(:,:), dmu(:,:), dlv(:,:)
    real(kdouble), allocatable :: dh1(:,:), dh1pre(:,:)
    real(kdouble), allocatable :: gW1(:,:), gWmu(:,:), gWlv(:,:), gW2(:,:), gW3(:,:)
    real(kdouble), allocatable :: gb1(:),  gbmu(:),  gblv(:),  gb2(:),  gb3(:)
    real(kdouble) :: invB, invBD, kl_w

    D = net%D; H = net%H; Z = net%Z; B = size(X, 2)
    invB  = 1.0d0 / real(B, kdouble)
    invBD = 1.0d0 / real(B*D, kdouble)
    kl_w  = beta * invB

    !> エンコーダ順伝播
    call monolis_alloc_R_2d(h1pre, H, B)
    call monolis_alloc_R_2d(h1,    H, B)
    h1pre = matmul(transpose(net%W1), X)
    do j = 1, B
      h1pre(:,j) = h1pre(:,j) + net%b1
    end do
    h1 = max(0.0d0, h1pre)

    call monolis_alloc_R_2d(mu,     Z, B)
    call monolis_alloc_R_2d(logvar, Z, B)
    mu     = matmul(transpose(net%Wmu), h1)
    logvar = matmul(transpose(net%Wlv), h1)
    do j = 1, B
      mu(:,j)     = mu(:,j)     + net%bmu
      logvar(:,j) = logvar(:,j) + net%blv
    end do
    where(logvar >  10.0d0) logvar =  10.0d0
    where(logvar < -10.0d0) logvar = -10.0d0

    !> 再パラメータ化
    call monolis_alloc_R_2d(eps,  Z, B)
    call monolis_alloc_R_2d(zlat, Z, B)
    call monolis_opt_vae_fill_randn(eps)
    zlat = mu + exp(0.5d0*logvar) * eps

    !> デコーダ順伝播
    call monolis_alloc_R_2d(h2pre,  H, B)
    call monolis_alloc_R_2d(h2,     H, B)
    call monolis_alloc_R_2d(logits, D, B)
    call monolis_alloc_R_2d(xhat,   D, B)
    h2pre = matmul(transpose(net%W2), zlat)
    do j = 1, B
      h2pre(:,j) = h2pre(:,j) + net%b2
    end do
    h2 = max(0.0d0, h2pre)
    logits = matmul(transpose(net%W3), h2)
    do j = 1, B
      logits(:,j) = logits(:,j) + net%b3
    end do
    xhat = 1.0d0 / (1.0d0 + exp(-logits))

    !> 損失 (再構成は特徴次元 D 方向にも平均: Keras mse 互換)
    recon_avg = sum((X - xhat)**2) * invBD
    kl_avg    = -0.5d0 * sum(1.0d0 + logvar - mu*mu - exp(logvar)) * invB
    loss_avg  = recon_avg + beta * kl_avg

    !> 逆伝播
    call monolis_alloc_R_2d(dxhat,   D, B)
    call monolis_alloc_R_2d(dlogits, D, B)
    dxhat   = 2.0d0 * (xhat - X) * invBD
    dlogits = dxhat * xhat * (1.0d0 - xhat)

    call monolis_alloc_R_2d(gW3, H, D)
    call monolis_alloc_R_1d(gb3, D)
    gW3 = matmul(h2, transpose(dlogits))
    gb3 = sum(dlogits, dim=2)

    call monolis_alloc_R_2d(dh2,    H, B)
    call monolis_alloc_R_2d(dh2pre, H, B)
    dh2 = matmul(net%W3, dlogits)
    where(h2pre > 0.0d0)
      dh2pre = dh2
    elsewhere
      dh2pre = 0.0d0
    end where

    call monolis_alloc_R_2d(gW2, Z, H)
    call monolis_alloc_R_1d(gb2, H)
    gW2 = matmul(zlat, transpose(dh2pre))
    gb2 = sum(dh2pre, dim=2)

    call monolis_alloc_R_2d(dz, Z, B)
    dz = matmul(net%W2, dh2pre)

    call monolis_alloc_R_2d(dmu, Z, B)
    call monolis_alloc_R_2d(dlv, Z, B)
    dmu = dz + kl_w * mu
    dlv = dz * eps * 0.5d0 * exp(0.5d0*logvar) + kl_w * 0.5d0*(exp(logvar) - 1.0d0)

    call monolis_alloc_R_2d(gWmu, H, Z)
    call monolis_alloc_R_1d(gbmu, Z)
    call monolis_alloc_R_2d(gWlv, H, Z)
    call monolis_alloc_R_1d(gblv, Z)
    gWmu = matmul(h1, transpose(dmu))
    gbmu = sum(dmu, dim=2)
    gWlv = matmul(h1, transpose(dlv))
    gblv = sum(dlv, dim=2)

    call monolis_alloc_R_2d(dh1,    H, B)
    call monolis_alloc_R_2d(dh1pre, H, B)
    dh1 = matmul(net%Wmu, dmu) + matmul(net%Wlv, dlv)
    where(h1pre > 0.0d0)
      dh1pre = dh1
    elsewhere
      dh1pre = 0.0d0
    end where

    call monolis_alloc_R_2d(gW1, D, H)
    call monolis_alloc_R_1d(gb1, H)
    gW1 = matmul(X, transpose(dh1pre))
    gb1 = sum(dh1pre, dim=2)

    !> Adam 更新
    net%t = net%t + 1
    call monolis_opt_vae_adam_apply(net%W1,  net%b1,  gW1,  gb1,  net%a_W1,  net%t, lr)
    call monolis_opt_vae_adam_apply(net%Wmu, net%bmu, gWmu, gbmu, net%a_Wmu, net%t, lr)
    call monolis_opt_vae_adam_apply(net%Wlv, net%blv, gWlv, gblv, net%a_Wlv, net%t, lr)
    call monolis_opt_vae_adam_apply(net%W2,  net%b2,  gW2,  gb2,  net%a_W2,  net%t, lr)
    call monolis_opt_vae_adam_apply(net%W3,  net%b3,  gW3,  gb3,  net%a_W3,  net%t, lr)
  end subroutine monolis_opt_vae_train_step

  !> @ingroup optimize
  !> ミニバッチ学習ループ (KL ウォームアップと EarlyStopping 付き)
  subroutine monolis_opt_vae_fit(net, X, opts)
    implicit none
    !> [in,out] 学習対象の VAE モデル
    type(monolis_opt_vae_t), intent(inout) :: net
    !> [in] 学習データ (D x N)
    real(kdouble), intent(in) :: X(:,:)
    !> [in] 学習オプション
    type(monolis_opt_vae_train_opts), intent(in) :: opts
    integer(kint) :: D, N, e, b, k, no_improve, nbatches, patience
    integer(kint), allocatable :: perm(:)
    real(kdouble), allocatable :: Xb(:,:)
    real(kdouble) :: loss, recon, kl, sum_loss, sum_recon, sum_kl
    real(kdouble) :: best_loss, epoch_loss, beta

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
    call monolis_alloc_R_2d(Xb, D, opts%batch_size)
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
  !> 決定論的 (z = mu) な順伝播による再構成
  subroutine monolis_opt_vae_reconstruct(net, X, xhat)
    implicit none
    !> [in] VAE モデル
    type(monolis_opt_vae_t), intent(in) :: net
    !> [in] 入力データ (D x B)
    real(kdouble), intent(in) :: X(:,:)
    !> [out] 再構成出力 (D x B)
    real(kdouble), intent(out) :: xhat(:,:)
    integer(kint) :: B, j
    real(kdouble), allocatable :: h1(:,:), mu(:,:), logvar(:,:), zlat(:,:)
    real(kdouble), allocatable :: h2(:,:), logits(:,:)

    B = size(X, 2)
    call monolis_alloc_R_2d(h1,     net%H, B)
    call monolis_alloc_R_2d(mu,     net%Z, B)
    call monolis_alloc_R_2d(logvar, net%Z, B)
    call monolis_alloc_R_2d(zlat,   net%Z, B)
    call monolis_alloc_R_2d(h2,     net%H, B)
    call monolis_alloc_R_2d(logits, net%D, B)

    h1 = matmul(transpose(net%W1), X)
    do j = 1, B
      h1(:,j) = h1(:,j) + net%b1
    end do
    h1 = max(0.0d0, h1)
    mu     = matmul(transpose(net%Wmu), h1)
    logvar = matmul(transpose(net%Wlv), h1)
    do j = 1, B
      mu(:,j)     = mu(:,j)     + net%bmu
      logvar(:,j) = logvar(:,j) + net%blv
    end do
    zlat = mu
    h2 = matmul(transpose(net%W2), zlat)
    do j = 1, B
      h2(:,j) = h2(:,j) + net%b2
    end do
    h2 = max(0.0d0, h2)
    logits = matmul(transpose(net%W3), h2)
    do j = 1, B
      logits(:,j) = logits(:,j) + net%b3
    end do
    xhat = 1.0d0 / (1.0d0 + exp(-logits))
  end subroutine monolis_opt_vae_reconstruct

  !> @ingroup optimize
  !> 入力をエンコードし、再パラメータ化サンプル z を返す
  subroutine monolis_opt_vae_encode(net, X, zlat)
    implicit none
    !> [in] VAE モデル
    type(monolis_opt_vae_t), intent(in) :: net
    !> [in] 入力データ (D x B)
    real(kdouble), intent(in) :: X(:,:)
    !> [out] サンプリングされた潜在ベクトル (Z x B)
    real(kdouble), intent(out) :: zlat(:,:)
    integer(kint) :: B, j
    real(kdouble), allocatable :: h1(:,:), mu(:,:), logvar(:,:), eps(:,:)

    B = size(X, 2)
    call monolis_alloc_R_2d(h1,     net%H, B)
    call monolis_alloc_R_2d(mu,     net%Z, B)
    call monolis_alloc_R_2d(logvar, net%Z, B)
    call monolis_alloc_R_2d(eps,    net%Z, B)

    h1 = matmul(transpose(net%W1), X)
    do j = 1, B
      h1(:,j) = h1(:,j) + net%b1
    end do
    h1 = max(0.0d0, h1)
    mu     = matmul(transpose(net%Wmu), h1)
    logvar = matmul(transpose(net%Wlv), h1)
    do j = 1, B
      mu(:,j)     = mu(:,j)     + net%bmu
      logvar(:,j) = logvar(:,j) + net%blv
    end do
    where(logvar >  10.0d0) logvar =  10.0d0
    where(logvar < -10.0d0) logvar = -10.0d0
    call monolis_opt_vae_fill_randn(eps)
    zlat = mu + exp(0.5d0*logvar) * eps
  end subroutine monolis_opt_vae_encode

  !> @ingroup optimize
  !> 潜在ベクトルをデコードして出力 (sigmoid) を得る
  subroutine monolis_opt_vae_decode(net, zlat, xhat)
    implicit none
    !> [in] VAE モデル
    type(monolis_opt_vae_t), intent(in) :: net
    !> [in] 潜在ベクトル (Z x B)
    real(kdouble), intent(in) :: zlat(:,:)
    !> [out] 出力 (D x B)
    real(kdouble), intent(out) :: xhat(:,:)
    integer(kint) :: B, j
    real(kdouble), allocatable :: h2(:,:), logits(:,:)

    B = size(zlat, 2)
    call monolis_alloc_R_2d(h2,     net%H, B)
    call monolis_alloc_R_2d(logits, net%D, B)

    h2 = matmul(transpose(net%W2), zlat)
    do j = 1, B
      h2(:,j) = h2(:,j) + net%b2
    end do
    h2 = max(0.0d0, h2)
    logits = matmul(transpose(net%W3), h2)
    do j = 1, B
      logits(:,j) = logits(:,j) + net%b3
    end do
    xhat = 1.0d0 / (1.0d0 + exp(-logits))
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
    real(kdouble), intent(out) :: xhat(:,:)
    real(kdouble), allocatable :: zlat(:,:)

    call monolis_alloc_R_2d(zlat, net%Z, n)
    call monolis_opt_vae_fill_randn(zlat)
    call monolis_opt_vae_decode(net, zlat, xhat)
  end subroutine monolis_opt_vae_sample_prior

  !> @ingroup optimize
  !> SPX (Simplex Crossover) で (Zdim+1) 個の親から 1 つの子を生成
  subroutine monolis_opt_vae_spx_child(parents, Zdim, child)
    implicit none
    !> [in] 親ベクトル (Zdim x (Zdim+1))
    real(kdouble), intent(in) :: parents(:,:)
    !> [in] 潜在次元数
    integer(kint), intent(in) :: Zdim
    !> [out] 子ベクトル (Zdim)
    real(kdouble), intent(out) :: child(:)
    integer(kint) :: i, Np
    real(kdouble) :: epsilon_, r
    real(kdouble), allocatable :: center(:), x(:,:), c(:,:)

    Np = Zdim + 1
    epsilon_ = sqrt(real(Zdim + 2, kdouble))
    call monolis_alloc_R_1d(center, Zdim)
    call monolis_alloc_R_2d(x,      Zdim, Np)
    call monolis_alloc_R_2d(c,      Zdim, Np)
    center = sum(parents, dim=2) / real(Np, kdouble)
    do i = 1, Np
      x(:, i) = center + epsilon_ * (parents(:, i) - center)
    end do
    c(:, 1) = 0.0d0
    do i = 2, Np
      call random_number(r)
      if(r < 1.0d-300) r = 1.0d-300
      r = r ** (1.0d0 / real(i - 1, kdouble))
      c(:, i) = r * (x(:, i-1) - x(:, i) + c(:, i-1))
    end do
    child = x(:, Np) + c(:, Np)
  end subroutine monolis_opt_vae_spx_child

  !> @ingroup optimize
  !> 学習データを潜在空間にエンコードしたうえで SPX 交叉により n 個生成・デコード
  subroutine monolis_opt_vae_generate_spx(net, Xtrain, n, xhat)
    implicit none
    !> [in] VAE モデル
    type(monolis_opt_vae_t), intent(in) :: net
    !> [in] 学習データ (D x Ntr)
    real(kdouble), intent(in) :: Xtrain(:,:)
    !> [in] 生成サンプル数
    integer(kint), intent(in) :: n
    !> [out] 出力 (D x n)
    real(kdouble), intent(out) :: xhat(:,:)
    integer(kint) :: Ntr, Np, i, k, idx
    real(kdouble) :: r
    real(kdouble), allocatable :: zall(:,:), parents(:,:), zchild(:,:)

    Ntr = size(Xtrain, 2)
    Np  = net%Z + 1
    call monolis_alloc_R_2d(zall, net%Z, Ntr)
    call monolis_opt_vae_encode(net, Xtrain, zall)
    call monolis_alloc_R_2d(parents, net%Z, Np)
    call monolis_alloc_R_2d(zchild,  net%Z, n)
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
