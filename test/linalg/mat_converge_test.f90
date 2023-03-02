!> 収束判定テストモジュール
module mod_monolis_converge_test
  use mod_monolis
  use mod_monolis_converge

  implicit none

contains

  subroutine monolis_converge_test()
    implicit none
    !> パラメータ構造体
    type(monolis_prm) :: monoPRM
    !> 通信テーブル構造体
    type(monolis_com) :: monoCOM
    !> 行列構造体
    type(monolis_mat) :: monoMAT
    !> 残差ベクトル
    real(kdouble) :: R(4)
    !> L2 相対誤差の分母値
    real(kdouble) :: B2
    !> 反復回数
    integer(kint) :: iter
    !> 収束判定フラグ
    logical :: is_converge
    !> 内積時間
    real(kdouble) :: tdotp
    !> 通信時間
    real(kdouble) :: tcomm

    call monolis_std_log_string("monolis_converge_test")

    monoMAT%N = 2
    monoMAT%NDOF = 2

    B2 = 4.0d0

    R(1) = 2.0d0
    R(2) = 2.0d0
    R(3) = 2.0d0
    R(4) = 2.0d0

    iter = 10

    call monolis_check_converge_R(monoPRM, monoCOM, monoMAT, R, B2, iter, is_converge, tdotp, tcomm)

    call monolis_test_check_eq_I1("monolis_converge_test 1", monoPRM%Iarray(monolis_prm_I_cur_iter), 10)
    call monolis_test_check_eq_R1("monolis_converge_test 2", monoPRM%Rarray(monolis_prm_R_cur_resid), 2.0d0)
    call monolis_test_check_eq_L1("monolis_converge_test 3", is_converge, .false.)
  end subroutine monolis_converge_test

end module mod_monolis_converge_test
