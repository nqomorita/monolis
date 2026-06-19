!> Adam オプティマイザライブラリ
!> @details VAE / Conditional VAE など MLP 学習で共有する Adam の状態と更新処理。
!>          1 つの重み行列 (in_dim x out_dim) とバイアス (out_dim) を対象とする。
module mod_monolis_opt_adam
  use mod_monolis_utils
  implicit none

  private
  public :: monolis_opt_adam_state
  public :: monolis_opt_adam_init
  public :: monolis_opt_adam_free
  public :: monolis_opt_adam_apply

  !> @ingroup optimize
  !> Adam オプティマイザの状態 (1 つの重み行列 + バイアスに対応)
  type :: monolis_opt_adam_state
    !> [in,out] 重み行列の 1 次モーメント
    real(kdouble), allocatable :: m(:,:)
    !> [in,out] 重み行列の 2 次モーメント
    real(kdouble), allocatable :: v(:,:)
    !> [in,out] バイアスの 1 次モーメント
    real(kdouble), allocatable :: mb(:)
    !> [in,out] バイアスの 2 次モーメント
    real(kdouble), allocatable :: vb(:)
  end type monolis_opt_adam_state

contains

  !> @ingroup optimize
  !> Adam 状態を初期化 (ゼロ埋め)
  subroutine monolis_opt_adam_init(st, n1, n2)
    implicit none
    !> [out] 初期化対象の Adam 状態
    type(monolis_opt_adam_state), intent(out) :: st
    !> [in] 重み行列の行数
    integer(kint), intent(in) :: n1
    !> [in] 重み行列の列数 (バイアス次元)
    integer(kint), intent(in) :: n2

    call monolis_alloc_R_2d(st%m,  n1, n2)
    call monolis_alloc_R_2d(st%v,  n1, n2)
    call monolis_alloc_R_1d(st%mb, n2)
    call monolis_alloc_R_1d(st%vb, n2)
  end subroutine monolis_opt_adam_init

  !> @ingroup optimize
  !> Adam 状態のメモリを解放
  subroutine monolis_opt_adam_free(st)
    implicit none
    !> [in,out] 解放対象の Adam 状態
    type(monolis_opt_adam_state), intent(inout) :: st

    if(allocated(st%m))  deallocate(st%m)
    if(allocated(st%v))  deallocate(st%v)
    if(allocated(st%mb)) deallocate(st%mb)
    if(allocated(st%vb)) deallocate(st%vb)
  end subroutine monolis_opt_adam_free

  !> @ingroup optimize
  !> Adam による 1 ステップ更新
  subroutine monolis_opt_adam_apply(W, b, gW, gb, st, t, lr)
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
    type(monolis_opt_adam_state), intent(inout) :: st
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
  end subroutine monolis_opt_adam_apply

end module mod_monolis_opt_adam
