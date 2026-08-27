!> Cholesky 前処理関連モジュール
module mod_monolis_precond_Cholesky
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_fact_Cholesky_nn

  implicit none

contains

  !> @ingroup prec
  !> 前処理生成：Cholesky 前処理（実数型）
  subroutine monolis_precond_Cholesky_setup_R(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC

    call monolis_std_debug_log_header("monolis_precond_Cholesky_setup_R")

    !> 保存済みの数値分解は、反復中の次の線形ソルバ呼び出しで再利用する
    if(monoPRM%Iarray(monolis_prm_I_is_prec_stored) == monolis_I_true .and. &
      & monoPREC%LU%factorized) return

    call monolis_fact_Cholesky_nn_setup_R(monoMAT, monoPREC)
  end subroutine monolis_precond_Cholesky_setup_R

  !> @ingroup prec
  !> 前処理生成：Cholesky 前処理（複素数型）
  subroutine monolis_precond_Cholesky_setup_C(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC

    call monolis_std_debug_log_header("monolis_precond_Cholesky_setup_C")

    !call monolis_fact_Cholesky_nn_setup_C(monoMAT, monoPREC)
  end subroutine monolis_precond_Cholesky_setup_C

  !> @ingroup prec
  !> 前処理適用：Cholesky 前処理（実数型）
  subroutine monolis_precond_Cholesky_apply_R(monoPRM, monoCOM, monoMAT, monoPREC, X, Y)
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

    call monolis_std_debug_log_header("monolis_precond_Cholesky_apply_R")

    call monolis_fact_Cholesky_nn_apply_R(monoMAT, monoPREC, X, Y)
  end subroutine monolis_precond_Cholesky_apply_R

  !> @ingroup prec
  !> 前処理適用：Cholesky 前処理（複素数型）
  subroutine monolis_precond_Cholesky_apply_C(monoPRM, monoCOM, monoMAT, monoPREC, X, Y)
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

    call monolis_std_debug_log_header("monolis_precond_Cholesky_apply_C")

    !call monolis_fact_Cholesky_nn_apply_C(monoMAT, monoPREC, X, Y)
  end subroutine monolis_precond_Cholesky_apply_C

  !> @ingroup prec
  !> 前処理初期化：Cholesky 前処理（実数型）
  subroutine monolis_precond_Cholesky_clear_R(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC

    call monolis_std_debug_log_header("monolis_precond_Cholesky_clear_R")

    !> 保存指定中は、線形ソルバ終了時にも分解データを保持する
    if(monoPRM%Iarray(monolis_prm_I_is_prec_stored) == monolis_I_true) return

    call monolis_fact_Cholesky_nn_clear_R(monoPREC)
  end subroutine monolis_precond_Cholesky_clear_R

  !> @ingroup prec
  !> 前処理初期化：Cholesky 前処理（複素数型）
  subroutine monolis_precond_Cholesky_clear_C(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC

    call monolis_std_debug_log_header("monolis_precond_Cholesky_clear_C")

    !call monolis_fact_Cholesky_nn_clear_C(monoPREC)
  end subroutine monolis_precond_Cholesky_clear_C
end module mod_monolis_precond_Cholesky
