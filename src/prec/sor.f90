!> SOR 前処理関連モジュール
module mod_monolis_precond_sor
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_precond_sor_33
  use mod_monolis_precond_sor_nn

  implicit none

contains

  !> 前処理生成：SOR 前処理
  subroutine monolis_precond_sor_setup(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> パラメータ構造体
    type(monolis_prm) :: monoPRM
    !> 通信テーブル構造体
    type(monolis_com) :: monoCOM
    !> 行列構造体
    type(monolis_mat) :: monoMAT
    !> 前処理構造体
    type(monolis_mat) :: monoPREC

    if(monoMAT%NDOF == 3)then
      call monolis_precond_sor_33_setup(monoMAT, monoPREC)
    else
      call monolis_precond_sor_nn_setup(monoMAT, monoPREC)
    endif
  end subroutine monolis_precond_sor_setup

  !> 前処理適用：SOR 前処理
  subroutine monolis_precond_sor_apply(monoPRM, monoCOM, monoMAT, monoPREC, X, Y)
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

    if(monoMAT%NDOF == 3)then
      call monolis_precond_sor_33_apply(monoMAT, monoPREC, X, Y)
    else
      call monolis_precond_sor_nn_apply(monoMAT, monoPREC, X, Y)
    endif
  end subroutine monolis_precond_sor_apply

  !> 前処理初期化：SOR 前処理
  subroutine monolis_precond_sor_clear(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> パラメータ構造体
    type(monolis_prm) :: monoPRM
    !> 通信テーブル構造体
    type(monolis_com) :: monoCOM
    !> 行列構造体
    type(monolis_mat) :: monoMAT
    !> 前処理構造体
    type(monolis_mat) :: monoPREC

    if(monoMAT%NDOF == 3)then
      call monolis_precond_sor_33_clear(monoPREC)
    else
      call monolis_precond_sor_nn_clear(monoPREC)
    endif
  end subroutine monolis_precond_sor_clear
end module mod_monolis_precond_sor