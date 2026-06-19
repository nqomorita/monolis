!> Adam オプティマイザライブラリ
!> @details VAE / Conditional VAE など MLP 学習で共有する Adam の状態と更新処理。
!>          1 つの重み行列 (in_dim x out_dim) とバイアス (out_dim) を対象とする。
!>          機械学習に関わる実数は 32bit 浮動小数点 (kdouble_ml) で計算する。
module mod_monolis_opt_adam
  use mod_monolis_utils
  use mod_monolis_def_opt
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
    real(kdouble_ml), allocatable :: m(:,:)
    !> [in,out] 重み行列の 2 次モーメント
    real(kdouble_ml), allocatable :: v(:,:)
    !> [in,out] バイアスの 1 次モーメント
    real(kdouble_ml), allocatable :: mb(:)
    !> [in,out] バイアスの 2 次モーメント
    real(kdouble_ml), allocatable :: vb(:)
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

    call monolis_alloc_F_2d(st%m,  n1, n2)
    call monolis_alloc_F_2d(st%v,  n1, n2)
    call monolis_alloc_F_1d(st%mb, n2)
    call monolis_alloc_F_1d(st%vb, n2)
  end subroutine monolis_opt_adam_init

  !> @ingroup optimize
  !> Adam 状態のメモリを解放
  subroutine monolis_opt_adam_free(st)
    implicit none
    !> [in,out] 解放対象の Adam 状態
    type(monolis_opt_adam_state), intent(inout) :: st

    call monolis_dealloc_F_2d(st%m)
    call monolis_dealloc_F_2d(st%v)
    call monolis_dealloc_F_1d(st%mb)
    call monolis_dealloc_F_1d(st%vb)
  end subroutine monolis_opt_adam_free

  !> @ingroup optimize
  !> Adam による 1 ステップ更新
  subroutine monolis_opt_adam_apply(W, b, gW, gb, st, t, lr)
    implicit none
    !> [in,out] 重み行列
    real(kdouble_ml), intent(inout) :: W(:,:)
    !> [in,out] バイアス
    real(kdouble_ml), intent(inout) :: b(:)
    !> [in] 重み勾配
    real(kdouble_ml), intent(in) :: gW(:,:)
    !> [in] バイアス勾配
    real(kdouble_ml), intent(in) :: gb(:)
    !> [in,out] Adam 状態
    type(monolis_opt_adam_state), intent(inout) :: st
    !> [in] グローバルステップ数 (1 始まり)
    integer(kint), intent(in) :: t
    !> [in] 学習率
    real(kdouble_ml), intent(in) :: lr
    real(kdouble_ml), parameter :: b1 = 0.9_kdouble_ml, b2 = 0.999_kdouble_ml, eps = 1.0e-7_kdouble_ml
    real(kdouble_ml) :: bc1, bc2

    st%m  = b1*st%m  + (1.0_kdouble_ml - b1)*gW
    st%v  = b2*st%v  + (1.0_kdouble_ml - b2)*gW*gW
    st%mb = b1*st%mb + (1.0_kdouble_ml - b1)*gb
    st%vb = b2*st%vb + (1.0_kdouble_ml - b2)*gb*gb
    bc1 = 1.0_kdouble_ml - b1**t
    bc2 = 1.0_kdouble_ml - b2**t
    W = W - lr * (st%m / bc1) / (sqrt(st%v / bc2) + eps)
    b = b - lr * (st%mb/ bc1) / (sqrt(st%vb/ bc2) + eps)
  end subroutine monolis_opt_adam_apply

end module mod_monolis_opt_adam
