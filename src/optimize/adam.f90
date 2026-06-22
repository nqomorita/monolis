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
    !> Adam 状態をデバイスに常駐させる (学習ループ中の転送を避ける)
    !$acc enter data copyin(st%m, st%v, st%mb, st%vb)
  end subroutine monolis_opt_adam_init

  !> @ingroup optimize
  !> Adam 状態のメモリを解放
  subroutine monolis_opt_adam_free(st)
    implicit none
    !> [in,out] 解放対象の Adam 状態
    type(monolis_opt_adam_state), intent(inout) :: st

    if(allocated(st%m))then
      !$acc exit data delete(st%m, st%v, st%mb, st%vb)
    endif
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
    real(kdouble_ml), parameter :: b1 = 0.9_kdouble_ml, b2 = 0.999_kdouble_ml
    real(kdouble_ml) :: bc1, bc2

    bc1 = 1.0_kdouble_ml - b1**t
    bc2 = 1.0_kdouble_ml - b2**t
    call monolis_opt_adam_apply_kernel(W, b, gW, gb, st%m, st%v, st%mb, st%vb, &
      size(W, 1), size(W, 2), lr, bc1, bc2)
  end subroutine monolis_opt_adam_apply

  !> @ingroup optimize
  !> Adam 更新カーネル (OpenACC オフロード本体)
  !> @details 重み W,b と Adam 状態 m,v,mb,vb は present-or-copy (copy 句) で受ける。
  !>          デバイス常駐時は転送せずデバイス上で更新し、非常駐時 (スタンドアロン
  !>          呼び出し) はホストと同期する。
  subroutine monolis_opt_adam_apply_kernel(W, b, gW, gb, m, v, mb, vb, n1, n2, lr, bc1, bc2)
    implicit none
    !> [in] 重み行列の行数
    integer(kint), intent(in) :: n1
    !> [in] 重み行列の列数 (バイアス次元)
    integer(kint), intent(in) :: n2
    !> [in,out] 重み行列
    real(kdouble_ml), intent(inout) :: W(n1, n2)
    !> [in,out] バイアス
    real(kdouble_ml), intent(inout) :: b(n2)
    !> [in] 重み勾配
    real(kdouble_ml), intent(in) :: gW(n1, n2)
    !> [in] バイアス勾配
    real(kdouble_ml), intent(in) :: gb(n2)
    !> [in,out] 重み 1 次モーメント
    real(kdouble_ml), intent(inout) :: m(n1, n2)
    !> [in,out] 重み 2 次モーメント
    real(kdouble_ml), intent(inout) :: v(n1, n2)
    !> [in,out] バイアス 1 次モーメント
    real(kdouble_ml), intent(inout) :: mb(n2)
    !> [in,out] バイアス 2 次モーメント
    real(kdouble_ml), intent(inout) :: vb(n2)
    !> [in] 学習率
    real(kdouble_ml), intent(in) :: lr
    !> [in] 1 次モーメントのバイアス補正項
    real(kdouble_ml), intent(in) :: bc1
    !> [in] 2 次モーメントのバイアス補正項
    real(kdouble_ml), intent(in) :: bc2
    real(kdouble_ml), parameter :: b1 = 0.9_kdouble_ml, b2 = 0.999_kdouble_ml, eps = 1.0e-7_kdouble_ml
    integer(kint) :: i, k

    !$acc parallel loop collapse(2) present_or_copy(W, m, v) present_or_copyin(gW)
    do k = 1, n2
      do i = 1, n1
        m(i,k) = b1*m(i,k) + (1.0_kdouble_ml - b1)*gW(i,k)
        v(i,k) = b2*v(i,k) + (1.0_kdouble_ml - b2)*gW(i,k)*gW(i,k)
        W(i,k) = W(i,k) - lr*(m(i,k)/bc1)/(sqrt(v(i,k)/bc2) + eps)
      end do
    end do
    !$acc end parallel loop
    !$acc parallel loop present_or_copy(b, mb, vb) present_or_copyin(gb)
    do k = 1, n2
      mb(k) = b1*mb(k) + (1.0_kdouble_ml - b1)*gb(k)
      vb(k) = b2*vb(k) + (1.0_kdouble_ml - b2)*gb(k)*gb(k)
      b(k) = b(k) - lr*(mb(k)/bc1)/(sqrt(vb(k)/bc2) + eps)
    end do
    !$acc end parallel loop
  end subroutine monolis_opt_adam_apply_kernel

end module mod_monolis_opt_adam
