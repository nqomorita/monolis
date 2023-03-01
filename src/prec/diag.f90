!> 対角スケーリング前処理関連モジュール
module mod_monolis_precond_diag
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_precond_diag_33
  use mod_monolis_precond_diag_nn

  implicit none

contains

  !> 前処理生成：対角スケーリング前処理
  subroutine monolis_precond_diag_setup(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> パラメータ構造体
    type(monolis_prm) :: monoPRM
    !> 通信テーブル構造体
    type(monolis_com) :: monoCOM
    !> 行列構造体
    type(monolis_mat) :: monoMAT
    !> 前処理構造体
    type(monolis_mat) :: monoPREC

    if(monoMAT%CSR%NDOF == 3)then
      call monolis_precond_diag_33_setup(monoMAT, monoPREC)
    else
      call monolis_precond_diag_nn_setup(monoMAT, monoPREC)
    endif
  end subroutine monolis_precond_diag_setup

  !> 前処理適用：対角スケーリング前処理
  subroutine monolis_precond_diag_apply(monoPRM, monoCOM, monoMAT, monoPREC, X, Y)
    implicit none
    !> パラメータ構造体
    type(monolis_prm) :: monoPRM
    !> 通信テーブル構造体
    type(monolis_com) :: monoCOM
    !> 行列構造体
    type(monolis_mat) :: monoMAT
    !> 前処理構造体
    type(monolis_mat) :: monoPREC
    real(kdouble) :: X(:), Y(:)

    if(monoMAT%CSR%NDOF == 3)then
      call monolis_precond_diag_33_apply(monoMAT, monoPREC, X, Y)
    else
      call monolis_precond_diag_nn_apply(monoMAT, monoPREC, X, Y)
    endif
  end subroutine monolis_precond_diag_apply

  !> 前処理初期化：対角スケーリング前処理
  subroutine monolis_precond_diag_clear(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> パラメータ構造体
    type(monolis_prm) :: monoPRM
    !> 通信テーブル構造体
    type(monolis_com) :: monoCOM
    !> 行列構造体
    type(monolis_mat) :: monoMAT
    !> 前処理構造体
    type(monolis_mat) :: monoPREC

    if(monoMAT%CSR%NDOF == 3)then
      call monolis_precond_diag_33_clear(monoPREC)
    else
      call monolis_precond_diag_nn_clear(monoPREC)
    endif
  end subroutine monolis_precond_diag_clear
end module mod_monolis_precond_diag