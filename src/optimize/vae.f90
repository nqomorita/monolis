!> 変分オートエンコーダ (VAE) ライブラリ
!> @details 任意層数 MLP のエンコーダ・デコーダをサポート。
!>          encoder : D -> H_enc(1) -> ... -> H_enc(L) (ReLU)
!>                    -> head_mu / head_lv (linear, Z)
!>          z       = mu + exp(0.5*logvar) * eps,  eps ~ N(0,1)
!>          decoder : Z -> H_dec(1) -> ... -> H_dec(M) (ReLU) -> D (sigmoid)
!>          loss    = mean_batch_dim( (x - xhat)^2 )
!>                  + beta * mean_batch( -0.5 * sum (1 + lv - mu^2 - exp(lv)) )
module mod_monolis_opt_vae
  use mod_monolis_utils
  implicit none

  private
  public :: monolis_opt_vae_adam_state
  public :: monolis_opt_vae_layer_t
  public :: monolis_opt_vae_t
  public :: monolis_opt_vae_train_opts
  public :: monolis_opt_vae_init
  public :: monolis_opt_vae_init_layers
  public :: monolis_opt_vae_finalize
  public :: monolis_opt_vae_train_step
  public :: monolis_opt_vae_fit
  public :: monolis_opt_vae_reconstruct
  public :: monolis_opt_vae_encode
  public :: monolis_opt_vae_decode
  public :: monolis_opt_vae_sample_prior
  public :: monolis_opt_vae_spx_child
  public :: monolis_opt_vae_generate_spx

  !> 活性関数種別
  integer(kint), parameter, public :: monolis_opt_vae_act_linear  = 0
  integer(kint), parameter, public :: monolis_opt_vae_act_relu    = 1
  integer(kint), parameter, public :: monolis_opt_vae_act_sigmoid = 2

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
  !> MLP の 1 層 (重み・バイアス・活性関数・Adam 状態)
  type :: monolis_opt_vae_layer_t
    !> [in] 入力次元
    integer(kint) :: in_dim  = 0
    !> [in] 出力次元
    integer(kint) :: out_dim = 0
    !> [in] 活性関数種別 (monolis_opt_vae_act_*)
    integer(kint) :: act     = 0
    !> [in,out] 重み行列 (in_dim x out_dim)
    real(kdouble), allocatable :: W(:,:)
    !> [in,out] バイアス (out_dim)
    real(kdouble), allocatable :: b(:)
    !> [in,out] Adam 状態
    type(monolis_opt_vae_adam_state) :: a
  end type monolis_opt_vae_layer_t

  !> 順伝播時のキャッシュ (層ごと)
  type :: monolis_opt_vae_cache_t
    !> 活性化前 (out_dim x B)
    real(kdouble), allocatable :: pre(:,:)
    !> 活性化後 (out_dim x B)
    real(kdouble), allocatable :: post(:,:)
  end type monolis_opt_vae_cache_t

  !> 1 層の勾配 (逆伝播時に動的に確保)
  type :: monolis_opt_vae_grad_t
    !> 重み勾配 (in_dim x out_dim)
    real(kdouble), allocatable :: gW(:,:)
    !> バイアス勾配 (out_dim)
    real(kdouble), allocatable :: gb(:)
  end type monolis_opt_vae_grad_t

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
  !> 1 層の確保と初期化 (Glorot 一様、バイアス 0、Adam 状態 0)
  subroutine monolis_opt_vae_layer_init(layer, in_dim, out_dim, act)
    implicit none
    !> [out] 初期化対象の層
    type(monolis_opt_vae_layer_t), intent(out) :: layer
    !> [in] 入力次元
    integer(kint), intent(in) :: in_dim
    !> [in] 出力次元
    integer(kint), intent(in) :: out_dim
    !> [in] 活性関数種別
    integer(kint), intent(in) :: act

    layer%in_dim  = in_dim
    layer%out_dim = out_dim
    layer%act     = act
    call monolis_alloc_R_2d(layer%W, in_dim, out_dim)
    call monolis_alloc_R_1d(layer%b, out_dim)
    call monolis_opt_vae_glorot_uniform(layer%W)
    call monolis_opt_vae_adam_init(layer%a, in_dim, out_dim)
  end subroutine monolis_opt_vae_layer_init

  !> @ingroup optimize
  !> 1 層のメモリ解放
  subroutine monolis_opt_vae_layer_free(layer)
    implicit none
    !> [in,out] 解放対象の層
    type(monolis_opt_vae_layer_t), intent(inout) :: layer

    if(allocated(layer%W)) deallocate(layer%W)
    if(allocated(layer%b)) deallocate(layer%b)
    call monolis_opt_vae_adam_free(layer%a)
    layer%in_dim = 0; layer%out_dim = 0; layer%act = 0
  end subroutine monolis_opt_vae_layer_free

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
  !> 1 層の順伝播 (post = act(W^T a_in + b)) とキャッシュ生成
  subroutine monolis_opt_vae_layer_forward(layer, a_in, cache)
    implicit none
    !> [in] 対象の層
    type(monolis_opt_vae_layer_t), intent(in) :: layer
    !> [in] 入力活性 (in_dim x B)
    real(kdouble), intent(in) :: a_in(:,:)
    !> [out] 活性化前/後をキャッシュ
    type(monolis_opt_vae_cache_t), intent(out) :: cache
    integer(kint) :: j, B

    B = size(a_in, 2)
    call monolis_alloc_R_2d(cache%pre,  layer%out_dim, B)
    call monolis_alloc_R_2d(cache%post, layer%out_dim, B)
    cache%pre = matmul(transpose(layer%W), a_in)
    do j = 1, B
      cache%pre(:,j) = cache%pre(:,j) + layer%b
    end do
    select case(layer%act)
    case(monolis_opt_vae_act_relu)
      cache%post = max(0.0d0, cache%pre)
    case(monolis_opt_vae_act_sigmoid)
      cache%post = 1.0d0 / (1.0d0 + exp(-cache%pre))
    case default
      cache%post = cache%pre
    end select
  end subroutine monolis_opt_vae_layer_forward

  !> @ingroup optimize
  !> 層列に対する順伝播
  subroutine monolis_opt_vae_forward_stack(layers, X, cache)
    implicit none
    !> [in] 層列
    type(monolis_opt_vae_layer_t), intent(in) :: layers(:)
    !> [in] 入力 (in_dim x B)
    real(kdouble), intent(in) :: X(:,:)
    !> [out] 層ごとのキャッシュ (size(layers))
    type(monolis_opt_vae_cache_t), allocatable, intent(out) :: cache(:)
    integer(kint) :: l, nL

    nL = size(layers)
    allocate(cache(nL))
    call monolis_opt_vae_layer_forward(layers(1), X, cache(1))
    do l = 2, nL
      call monolis_opt_vae_layer_forward(layers(l), cache(l-1)%post, cache(l))
    end do
  end subroutine monolis_opt_vae_forward_stack

  !> @ingroup optimize
  !> 層列の最終キャッシュ (post) を返すヘルパ
  function monolis_opt_vae_top_post(cache) result(p)
    implicit none
    !> [in] 層列キャッシュ
    type(monolis_opt_vae_cache_t), intent(in), target :: cache(:)
    !> 最終層の post への配列コピー
    real(kdouble), allocatable :: p(:,:)
    integer(kint) :: nL

    nL = size(cache)
    allocate(p(size(cache(nL)%post,1), size(cache(nL)%post,2)))
    p = cache(nL)%post
  end function monolis_opt_vae_top_post

  !> @ingroup optimize
  !> 1 層の逆伝播 (post 側勾配 dY -> pre 側勾配と重み勾配、入力側勾配)
  subroutine monolis_opt_vae_layer_backward(layer, a_in, cache, dY, gW, gb, dX_in)
    implicit none
    !> [in] 対象の層
    type(monolis_opt_vae_layer_t), intent(in) :: layer
    !> [in] 入力活性 (in_dim x B)
    real(kdouble), intent(in) :: a_in(:,:)
    !> [in] 順伝播キャッシュ
    type(monolis_opt_vae_cache_t), intent(in) :: cache
    !> [in] 上流からの勾配 (out_dim x B、post に対する勾配)
    real(kdouble), intent(in) :: dY(:,:)
    !> [out] 重み勾配 (in_dim x out_dim)
    real(kdouble), allocatable, intent(out) :: gW(:,:)
    !> [out] バイアス勾配 (out_dim)
    real(kdouble), allocatable, intent(out) :: gb(:)
    !> [out] 入力側勾配 (in_dim x B)
    real(kdouble), allocatable, intent(out) :: dX_in(:,:)
    real(kdouble), allocatable :: dpre(:,:)
    integer(kint) :: B

    B = size(a_in, 2)
    call monolis_alloc_R_2d(dpre, layer%out_dim, B)
    select case(layer%act)
    case(monolis_opt_vae_act_relu)
      where(cache%pre > 0.0d0)
        dpre = dY
      elsewhere
        dpre = 0.0d0
      end where
    case(monolis_opt_vae_act_sigmoid)
      dpre = dY * cache%post * (1.0d0 - cache%post)
    case default
      dpre = dY
    end select
    call monolis_alloc_R_2d(gW, layer%in_dim, layer%out_dim)
    call monolis_alloc_R_1d(gb, layer%out_dim)
    gW = matmul(a_in, transpose(dpre))
    gb = sum(dpre, dim=2)
    call monolis_alloc_R_2d(dX_in, layer%in_dim, B)
    dX_in = matmul(layer%W, dpre)
    deallocate(dpre)
  end subroutine monolis_opt_vae_layer_backward

  !> @ingroup optimize
  !> 層列に対する逆伝播 (最終層側 dY_top -> 入力側 dX_in、層ごとの勾配)
  subroutine monolis_opt_vae_backward_stack(layers, X, cache, dY_top, grads, dX_in)
    implicit none
    !> [in] 層列
    type(monolis_opt_vae_layer_t), intent(in) :: layers(:)
    !> [in] 入力 (in_dim x B)
    real(kdouble), intent(in) :: X(:,:)
    !> [in] 順伝播キャッシュ
    type(monolis_opt_vae_cache_t), intent(in) :: cache(:)
    !> [in] 最終層 post への勾配
    real(kdouble), intent(in) :: dY_top(:,:)
    !> [out] 層ごとの勾配
    type(monolis_opt_vae_grad_t), allocatable, intent(out) :: grads(:)
    !> [out] 入力側勾配
    real(kdouble), allocatable, intent(out) :: dX_in(:,:)
    real(kdouble), allocatable :: dY(:,:), dX(:,:)
    integer(kint) :: l, nL

    nL = size(layers)
    allocate(grads(nL))
    allocate(dY(size(dY_top,1), size(dY_top,2)))
    dY = dY_top
    do l = nL, 1, -1
      if(l == 1)then
        call monolis_opt_vae_layer_backward(layers(l), X, cache(l), dY, &
          grads(l)%gW, grads(l)%gb, dX)
      else
        call monolis_opt_vae_layer_backward(layers(l), cache(l-1)%post, cache(l), dY, &
          grads(l)%gW, grads(l)%gb, dX)
      endif
      deallocate(dY)
      call move_alloc(dX, dY)
    end do
    call move_alloc(dY, dX_in)
  end subroutine monolis_opt_vae_backward_stack

  !> @ingroup optimize
  !> キャッシュ列を解放する
  subroutine monolis_opt_vae_cache_free(cache)
    implicit none
    !> [in,out] 解放対象
    type(monolis_opt_vae_cache_t), allocatable, intent(inout) :: cache(:)
    integer(kint) :: l

    if(.not. allocated(cache)) return
    do l = 1, size(cache)
      if(allocated(cache(l)%pre))  deallocate(cache(l)%pre)
      if(allocated(cache(l)%post)) deallocate(cache(l)%post)
    end do
    deallocate(cache)
  end subroutine monolis_opt_vae_cache_free

  !> @ingroup optimize
  !> 勾配列を解放する
  subroutine monolis_opt_vae_grads_free(grads)
    implicit none
    !> [in,out] 解放対象
    type(monolis_opt_vae_grad_t), allocatable, intent(inout) :: grads(:)
    integer(kint) :: l

    if(.not. allocated(grads)) return
    do l = 1, size(grads)
      if(allocated(grads(l)%gW)) deallocate(grads(l)%gW)
      if(allocated(grads(l)%gb)) deallocate(grads(l)%gb)
    end do
    deallocate(grads)
  end subroutine monolis_opt_vae_grads_free

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
    integer(kint) :: D, Z, B, j, l, Le, Ld
    type(monolis_opt_vae_cache_t), allocatable :: enc_cache(:), dec_cache(:)
    type(monolis_opt_vae_grad_t),  allocatable :: enc_grads(:), dec_grads(:)
    real(kdouble), allocatable :: h_last(:,:)
    real(kdouble), allocatable :: mu(:,:), logvar(:,:), eps(:,:), zlat(:,:), xhat(:,:)
    real(kdouble), allocatable :: dxhat(:,:), dz_dec(:,:)
    real(kdouble), allocatable :: dmu(:,:), dlv(:,:)
    real(kdouble), allocatable :: gWmu(:,:), gWlv(:,:), gbmu(:), gblv(:)
    real(kdouble), allocatable :: dh1_mu(:,:), dh1_lv(:,:), dh_last(:,:), dX_dummy(:,:)
    real(kdouble) :: invB, invBD, kl_w

    D = net%D; Z = net%Z; B = size(X, 2)
    Le = size(net%enc); Ld = size(net%dec)
    invB  = 1.0d0 / real(B, kdouble)
    invBD = 1.0d0 / real(B*D, kdouble)
    kl_w  = beta * invB

    !> エンコーダ順伝播
    call monolis_opt_vae_forward_stack(net%enc, X, enc_cache)
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

    !> デコーダ順伝播
    call monolis_opt_vae_forward_stack(net%dec, zlat, dec_cache)
    xhat = monolis_opt_vae_top_post(dec_cache)

    !> 損失
    recon_avg = sum((X - xhat)**2) * invBD
    kl_avg    = -0.5d0 * sum(1.0d0 + logvar - mu*mu - exp(logvar)) * invB
    loss_avg  = recon_avg + beta * kl_avg

    !> デコーダ逆伝播: 出力 post=xhat に対する勾配 (sigmoid は層内で処理)
    call monolis_alloc_R_2d(dxhat, D, B)
    dxhat = 2.0d0 * (xhat - X) * invBD
    call monolis_opt_vae_backward_stack(net%dec, zlat, dec_cache, dxhat, dec_grads, dz_dec)

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
    call monolis_opt_vae_backward_stack(net%enc, X, enc_cache, dh_last, enc_grads, dX_dummy)

    !> Adam 更新
    net%t = net%t + 1
    do l = 1, Le
      call monolis_opt_vae_adam_apply(net%enc(l)%W, net%enc(l)%b, &
        enc_grads(l)%gW, enc_grads(l)%gb, net%enc(l)%a, net%t, lr)
    end do
    call monolis_opt_vae_adam_apply(net%head_mu%W, net%head_mu%b, gWmu, gbmu, &
      net%head_mu%a, net%t, lr)
    call monolis_opt_vae_adam_apply(net%head_lv%W, net%head_lv%b, gWlv, gblv, &
      net%head_lv%a, net%t, lr)
    do l = 1, Ld
      call monolis_opt_vae_adam_apply(net%dec(l)%W, net%dec(l)%b, &
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
  !> 共通の順伝播: X -> (enc + heads) -> mu/logvar/h_last
  subroutine monolis_opt_vae_encode_mu_lv(net, X, mu, logvar, enc_cache, h_last)
    implicit none
    !> [in] VAE モデル
    type(monolis_opt_vae_t), intent(in) :: net
    !> [in] 入力 (D x B)
    real(kdouble), intent(in) :: X(:,:)
    !> [out] mu (Z x B)
    real(kdouble), allocatable, intent(out) :: mu(:,:)
    !> [out] logvar (Z x B、[-10,10] でクリップ)
    real(kdouble), allocatable, intent(out) :: logvar(:,:)
    !> [out] エンコーダキャッシュ (不要なら呼び出し後に解放)
    type(monolis_opt_vae_cache_t), allocatable, intent(out) :: enc_cache(:)
    !> [out] 最終エンコーダ隠れ層出力
    real(kdouble), allocatable, intent(out) :: h_last(:,:)
    integer(kint) :: j, B

    B = size(X, 2)
    call monolis_opt_vae_forward_stack(net%enc, X, enc_cache)
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
  end subroutine monolis_opt_vae_encode_mu_lv

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
    real(kdouble), allocatable :: mu(:,:), logvar(:,:), h_last(:,:), xhat_loc(:,:)
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
    real(kdouble), intent(in) :: X(:,:)
    !> [out] サンプリングされた潜在ベクトル (Z x B)
    real(kdouble), intent(out) :: zlat(:,:)
    real(kdouble), allocatable :: mu(:,:), logvar(:,:), eps(:,:), h_last(:,:)
    type(monolis_opt_vae_cache_t), allocatable :: enc_cache(:)

    call monolis_opt_vae_encode_mu_lv(net, X, mu, logvar, enc_cache, h_last)
    call monolis_alloc_R_2d(eps, net%Z, size(X, 2))
    call monolis_opt_vae_fill_randn(eps)
    zlat = mu + exp(0.5d0*logvar) * eps
    call monolis_opt_vae_cache_free(enc_cache)
  end subroutine monolis_opt_vae_encode

  !> @ingroup optimize
  !> デコード本体 (戻り値は allocatable)
  subroutine monolis_opt_vae_decode_alloc(net, zlat, xhat)
    implicit none
    !> [in] VAE モデル
    type(monolis_opt_vae_t), intent(in) :: net
    !> [in] 潜在ベクトル (Z x B)
    real(kdouble), intent(in) :: zlat(:,:)
    !> [out] 出力 (D x B)
    real(kdouble), allocatable, intent(out) :: xhat(:,:)
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
    real(kdouble), intent(in) :: zlat(:,:)
    !> [out] 出力 (D x B)
    real(kdouble), intent(out) :: xhat(:,:)
    real(kdouble), allocatable :: xhat_loc(:,:)

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
