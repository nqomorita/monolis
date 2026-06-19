!> VAE / Conditional VAE 共通ユーティリティライブラリ
!> @details MLP の層構造・順伝播・逆伝播・初期化・乱数など、
!>          VAE と Conditional VAE で共有するデータ構造と処理をまとめる。
!>          機械学習に関わる実数は 32bit 浮動小数点 (kdouble_ml) で計算する。
module mod_monolis_opt_vae_util
  use mod_monolis_utils
  use mod_monolis_def_opt
  use mod_monolis_opt_adam
  implicit none

  private
  public :: monolis_opt_vae_layer_t
  public :: monolis_opt_vae_cache_t
  public :: monolis_opt_vae_grad_t
  public :: monolis_opt_vae_train_opts
  public :: monolis_opt_vae_act_linear
  public :: monolis_opt_vae_act_relu
  public :: monolis_opt_vae_act_sigmoid
  public :: monolis_opt_vae_layer_init
  public :: monolis_opt_vae_layer_free
  public :: monolis_opt_vae_fill_randn
  public :: monolis_opt_vae_shuffle
  public :: monolis_opt_vae_forward_stack
  public :: monolis_opt_vae_top_post
  public :: monolis_opt_vae_backward_stack
  public :: monolis_opt_vae_cache_free
  public :: monolis_opt_vae_grads_free
  public :: monolis_opt_vae_spx_child

  !> 活性関数種別
  integer(kint), parameter :: monolis_opt_vae_act_linear  = 0
  integer(kint), parameter :: monolis_opt_vae_act_relu    = 1
  integer(kint), parameter :: monolis_opt_vae_act_sigmoid = 2

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
    real(kdouble_ml), allocatable :: W(:,:)
    !> [in,out] バイアス (out_dim)
    real(kdouble_ml), allocatable :: b(:)
    !> [in,out] Adam 状態
    type(monolis_opt_adam_state) :: a
  end type monolis_opt_vae_layer_t

  !> 順伝播時のキャッシュ (層ごと)
  type :: monolis_opt_vae_cache_t
    !> 活性化前 (out_dim x B)
    real(kdouble_ml), allocatable :: pre(:,:)
    !> 活性化後 (out_dim x B)
    real(kdouble_ml), allocatable :: post(:,:)
  end type monolis_opt_vae_cache_t

  !> 1 層の勾配 (逆伝播時に動的に確保)
  type :: monolis_opt_vae_grad_t
    !> 重み勾配 (in_dim x out_dim)
    real(kdouble_ml), allocatable :: gW(:,:)
    !> バイアス勾配 (out_dim)
    real(kdouble_ml), allocatable :: gb(:)
  end type monolis_opt_vae_grad_t

  !> @ingroup optimize
  !> 学習ハイパーパラメータ
  type :: monolis_opt_vae_train_opts
    !> [in] ミニバッチサイズ
    integer(kint) :: batch_size = 128
    !> [in] エポック数
    integer(kint) :: epochs = 10
    !> [in] 学習率
    real(kdouble_ml) :: lr = 1.0e-3_kdouble_ml
    !> [in] KL 損失の最終スケール (beta)
    real(kdouble_ml) :: r_loss_factor = 1.0_kdouble_ml
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
    allocate(layer%W(in_dim, out_dim), source = 0.0_kdouble_ml)
    allocate(layer%b(out_dim),         source = 0.0_kdouble_ml)
    call monolis_opt_vae_glorot_uniform(layer%W)
    call monolis_opt_adam_init(layer%a, in_dim, out_dim)
  end subroutine monolis_opt_vae_layer_init

  !> @ingroup optimize
  !> 1 層のメモリ解放
  subroutine monolis_opt_vae_layer_free(layer)
    implicit none
    !> [in,out] 解放対象の層
    type(monolis_opt_vae_layer_t), intent(inout) :: layer

    if(allocated(layer%W)) deallocate(layer%W)
    if(allocated(layer%b)) deallocate(layer%b)
    call monolis_opt_adam_free(layer%a)
    layer%in_dim = 0; layer%out_dim = 0; layer%act = 0
  end subroutine monolis_opt_vae_layer_free

  !> @ingroup optimize
  !> Glorot/Xavier 一様初期化 : U(-sqrt(6/(in+out)), +sqrt(6/(in+out)))
  subroutine monolis_opt_vae_glorot_uniform(W)
    implicit none
    !> [out] 初期化対象の重み行列 (fan_in x fan_out)
    real(kdouble_ml), intent(out) :: W(:,:)
    real(kdouble_ml) :: limit, r
    integer(kint) :: i, j, fan_in, fan_out

    fan_in  = size(W, 1)
    fan_out = size(W, 2)
    limit = sqrt(6.0_kdouble_ml / real(fan_in + fan_out, kdouble_ml))
    do j = 1, size(W, 2)
      do i = 1, size(W, 1)
        call random_number(r)
        W(i,j) = (2.0_kdouble_ml*r - 1.0_kdouble_ml) * limit
      end do
    end do
  end subroutine monolis_opt_vae_glorot_uniform

  !> @ingroup optimize
  !> Box-Muller による標準正規乱数で配列を埋める
  subroutine monolis_opt_vae_fill_randn(a)
    implicit none
    !> [out] 充填対象の 2 次元配列
    real(kdouble_ml), intent(out) :: a(:,:)
    real(kdouble_ml), parameter :: pi = 3.141592653589793_kdouble_ml
    real(kdouble_ml) :: u1, u2
    integer(kint) :: i, j

    do j = 1, size(a, 2)
      do i = 1, size(a, 1)
        call random_number(u1)
        call random_number(u2)
        if(u1 < 1.0e-12_kdouble_ml) u1 = 1.0e-12_kdouble_ml
        a(i,j) = sqrt(-2.0_kdouble_ml*log(u1)) * cos(2.0_kdouble_ml*pi*u2)
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
    real(kdouble_ml) :: r

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
    real(kdouble_ml), intent(in) :: a_in(:,:)
    !> [out] 活性化前/後をキャッシュ
    type(monolis_opt_vae_cache_t), intent(out) :: cache
    integer(kint) :: j, B

    B = size(a_in, 2)
    allocate(cache%pre(layer%out_dim, B),  source = 0.0_kdouble_ml)
    allocate(cache%post(layer%out_dim, B), source = 0.0_kdouble_ml)
    cache%pre = matmul(transpose(layer%W), a_in)
    do j = 1, B
      cache%pre(:,j) = cache%pre(:,j) + layer%b
    end do
    select case(layer%act)
    case(monolis_opt_vae_act_relu)
      cache%post = max(0.0_kdouble_ml, cache%pre)
    case(monolis_opt_vae_act_sigmoid)
      cache%post = 1.0_kdouble_ml / (1.0_kdouble_ml + exp(-cache%pre))
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
    real(kdouble_ml), intent(in) :: X(:,:)
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
    real(kdouble_ml), allocatable :: p(:,:)
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
    real(kdouble_ml), intent(in) :: a_in(:,:)
    !> [in] 順伝播キャッシュ
    type(monolis_opt_vae_cache_t), intent(in) :: cache
    !> [in] 上流からの勾配 (out_dim x B、post に対する勾配)
    real(kdouble_ml), intent(in) :: dY(:,:)
    !> [out] 重み勾配 (in_dim x out_dim)
    real(kdouble_ml), allocatable, intent(out) :: gW(:,:)
    !> [out] バイアス勾配 (out_dim)
    real(kdouble_ml), allocatable, intent(out) :: gb(:)
    !> [out] 入力側勾配 (in_dim x B)
    real(kdouble_ml), allocatable, intent(out) :: dX_in(:,:)
    real(kdouble_ml), allocatable :: dpre(:,:)
    integer(kint) :: B

    B = size(a_in, 2)
    allocate(dpre(layer%out_dim, B), source = 0.0_kdouble_ml)
    select case(layer%act)
    case(monolis_opt_vae_act_relu)
      where(cache%pre > 0.0_kdouble_ml)
        dpre = dY
      elsewhere
        dpre = 0.0_kdouble_ml
      end where
    case(monolis_opt_vae_act_sigmoid)
      dpre = dY * cache%post * (1.0_kdouble_ml - cache%post)
    case default
      dpre = dY
    end select
    allocate(gW(layer%in_dim, layer%out_dim), source = 0.0_kdouble_ml)
    allocate(gb(layer%out_dim),               source = 0.0_kdouble_ml)
    gW = matmul(a_in, transpose(dpre))
    gb = sum(dpre, dim=2)
    allocate(dX_in(layer%in_dim, B), source = 0.0_kdouble_ml)
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
    real(kdouble_ml), intent(in) :: X(:,:)
    !> [in] 順伝播キャッシュ
    type(monolis_opt_vae_cache_t), intent(in) :: cache(:)
    !> [in] 最終層 post への勾配
    real(kdouble_ml), intent(in) :: dY_top(:,:)
    !> [out] 層ごとの勾配
    type(monolis_opt_vae_grad_t), allocatable, intent(out) :: grads(:)
    !> [out] 入力側勾配
    real(kdouble_ml), allocatable, intent(out) :: dX_in(:,:)
    real(kdouble_ml), allocatable :: dY(:,:), dX(:,:)
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
  !> SPX (Simplex Crossover) で (Zdim+1) 個の親から 1 つの子を生成
  subroutine monolis_opt_vae_spx_child(parents, Zdim, child)
    implicit none
    !> [in] 親ベクトル (Zdim x (Zdim+1))
    real(kdouble_ml), intent(in) :: parents(:,:)
    !> [in] 潜在次元数
    integer(kint), intent(in) :: Zdim
    !> [out] 子ベクトル (Zdim)
    real(kdouble_ml), intent(out) :: child(:)
    integer(kint) :: i, Np
    real(kdouble_ml) :: epsilon_, r
    real(kdouble_ml), allocatable :: center(:), x(:,:), c(:,:)

    Np = Zdim + 1
    epsilon_ = sqrt(real(Zdim + 2, kdouble_ml))
    allocate(center(Zdim),  source = 0.0_kdouble_ml)
    allocate(x(Zdim, Np),   source = 0.0_kdouble_ml)
    allocate(c(Zdim, Np),   source = 0.0_kdouble_ml)
    center = sum(parents, dim=2) / real(Np, kdouble_ml)
    do i = 1, Np
      x(:, i) = center + epsilon_ * (parents(:, i) - center)
    end do
    c(:, 1) = 0.0_kdouble_ml
    do i = 2, Np
      call random_number(r)
      if(r < 1.0e-30_kdouble_ml) r = 1.0e-30_kdouble_ml
      r = r ** (1.0_kdouble_ml / real(i - 1, kdouble_ml))
      c(:, i) = r * (x(:, i-1) - x(:, i) + c(:, i-1))
    end do
    child = x(:, Np) + c(:, Np)
  end subroutine monolis_opt_vae_spx_child

end module mod_monolis_opt_vae_util
