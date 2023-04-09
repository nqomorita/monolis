!> LAPACK ラッパーモジュール
module mod_monolis_lapack_test
  use mod_monolis_utils
  use mod_monolis_lapack

  implicit none

contains

  subroutine monolis_lapack_test()
    implicit none
    !> 行列の大きさ
    integer(kint) :: n
    !> 行列の対角成分
    real(kdouble) :: D(5)
    !> 行列の副対角成分
    real(kdouble) :: S(5)
    !> 固有値
    real(kdouble) :: eig_val(5)
    !> 固有ベクトル
    real(kdouble) :: eig_mode(5,5)
    real(kdouble) :: r_ans(5)

    call monolis_std_global_log_string("monolis_lapack_dstev")

    n = 5

    D(1) = 2.0d0
    D(2) = 2.0d0
    D(3) = 2.0d0
    D(4) = 2.0d0
    D(5) = 2.0d0

    S(1) = 1.0d0
    S(2) = 1.0d0
    S(3) = 1.0d0
    S(4) = 1.0d0
    S(5) = 0.0d0

    call monolis_lapack_dstev(n, D, S, eig_val, eig_mode)

    call monolis_test_check_eq_R1("monolis_lapack_test 1", eig_val(1), 0.267949192431122d0)
    call monolis_test_check_eq_R1("monolis_lapack_test 2", eig_val(2), 1.0d0)
    call monolis_test_check_eq_R1("monolis_lapack_test 3", eig_val(3), 2.0d0)
    call monolis_test_check_eq_R1("monolis_lapack_test 4", eig_val(4), 3.0d0)
    call monolis_test_check_eq_R1("monolis_lapack_test 5", eig_val(5), 3.732050807568877d0)

    r_ans(1) = 0.28867513459481281d0
    r_ans(2) =-0.5d0
    r_ans(3) = 0.57735026918962640d0
    r_ans(4) =-0.5d0
    r_ans(5) = 0.28867513459481270d0

    call monolis_test_check_eq_R("monolis_lapack_test 1 b", dabs(r_ans), dabs(eig_mode(:,1)))

    r_ans(1) =-0.5d0
    r_ans(2) = 0.5d0
    r_ans(3) = 0.0d0
    r_ans(4) =-0.5d0
    r_ans(5) = 0.5d0

    call monolis_test_check_eq_R("monolis_lapack_test 2 b", dabs(r_ans), dabs(eig_mode(:,2)))

    r_ans(1) = 0.5773502691896257d0
    r_ans(2) = 0.0d0
    r_ans(3) =-0.5773502691896257d0
    r_ans(4) = 0.0d0
    r_ans(5) = 0.5773502691896257d0

    call monolis_test_check_eq_R("monolis_lapack_test 3 b", dabs(r_ans), dabs(eig_mode(:,3)))

    r_ans(1) = 0.5d0
    r_ans(2) = 0.5d0
    r_ans(3) = 0.0d0
    r_ans(4) =-0.5d0
    r_ans(5) =-0.5d0

    call monolis_test_check_eq_R("monolis_lapack_test 4 b", dabs(r_ans), dabs(eig_mode(:,4)))

    r_ans(1) =-0.28867513459481281d0
    r_ans(2) =-0.5d0
    r_ans(3) =-0.57735026918962640d0
    r_ans(4) =-0.5d0
    r_ans(5) =-0.28867513459481270d0

    call monolis_test_check_eq_R("monolis_lapack_test 5 b", dabs(r_ans), dabs(eig_mode(:,5)))
  end subroutine monolis_lapack_test

end module mod_monolis_lapack_test
