!> SOR 前処理関連モジュール
module mod_monolis_precond_sor
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_precond_sor_33
  use mod_monolis_precond_sor_nn

  implicit none

contains

  !> @ingroup prec
  !> 前処理生成：SOR 前処理（実数型）
  subroutine monolis_precond_sor_setup_R(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC

    call monolis_std_debug_log_header("monolis_precond_sor_setup_R")

    if(monoMAT%NDOF == 3)then
      call monolis_precond_sor_33_setup_R(monoMAT, monoPREC)
    else
      call monolis_precond_sor_nn_setup_R(monoMAT, monoPREC)
    endif
  end subroutine monolis_precond_sor_setup_R

  !> @ingroup prec
  !> 前処理生成：SOR 前処理（複素数型）
  subroutine monolis_precond_sor_setup_C(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC

    call monolis_std_debug_log_header("monolis_precond_sor_setup_C")

    call monolis_precond_sor_nn_setup_C(monoMAT, monoPREC)
  end subroutine monolis_precond_sor_setup_C

  !> @ingroup prec
  !> 前処理適用：SOR 前処理（実数型）
  subroutine monolis_precond_sor_apply_R(monoPRM, monoCOM, monoMAT, monoPREC, X, Y)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC
    real(kdouble) :: X(:), Y(:)

    call monolis_std_debug_log_header("monolis_precond_sor_apply_R")

    if(monoMAT%NDOF == 3)then
      call monolis_precond_sor_33_apply_R(monoMAT, monoPREC, X, Y)
    else
      call monolis_precond_sor_nn_apply_R(monoMAT, monoPREC, X, Y)
    endif
  end subroutine monolis_precond_sor_apply_R

  !> @ingroup prec
  !> 前処理適用：SOR 前処理（複素数型）
  subroutine monolis_precond_sor_apply_C(monoPRM, monoCOM, monoMAT, monoPREC, X, Y)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC
    complex(kdouble) :: X(:), Y(:)

    call monolis_std_debug_log_header("monolis_precond_sor_apply_C")

    call monolis_precond_sor_nn_apply_C(monoMAT, monoPREC, X, Y)
  end subroutine monolis_precond_sor_apply_C

  !> @ingroup prec
  !> 前処理初期化：SOR 前処理（実数型）
  subroutine monolis_precond_sor_clear_R(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC

    call monolis_std_debug_log_header("monolis_precond_sor_clear_R")

    if(monoMAT%NDOF == 3)then
      call monolis_precond_sor_33_clear_R(monoPREC)
    else
      call monolis_precond_sor_nn_clear_R(monoPREC)
    endif
  end subroutine monolis_precond_sor_clear_R

  !> @ingroup prec
  !> 前処理初期化：SOR 前処理（複素数型）
  subroutine monolis_precond_sor_clear_C(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC

    call monolis_std_debug_log_header("monolis_precond_sor_clear_C")

    call monolis_precond_sor_nn_clear_R(monoPREC)
  end subroutine monolis_precond_sor_clear_C
end module mod_monolis_precond_sor