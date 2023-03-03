!> 収束判定テストモジュール
module mod_monolis_converge_test
  use mod_monolis
  use mod_monolis_converge

  implicit none

contains

  subroutine monolis_converge_test()
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    real(kdouble) :: R(4)
    real(kdouble) :: B2
    integer(kint) :: iter
    logical :: is_converge
    real(kdouble) :: tdotp
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
