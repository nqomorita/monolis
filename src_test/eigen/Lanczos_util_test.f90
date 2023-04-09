!> Lanczos 法 util テストモジュール
module mod_monolis_eigen_lanczos_util_test
  use mod_monolis
  use mod_monolis_eigen_lanczos_util

  implicit none

contains

  subroutine monolis_eigen_lanczos_util_test()
    implicit none

    call monolis_get_eigen_pair_from_tridiag_test()
  end subroutine monolis_eigen_lanczos_util_test

  subroutine monolis_get_eigen_pair_from_tridiag_test()
    implicit none
    integer(kint) :: iter
    integer(kint) :: n_get_eigen
    real(kdouble) :: alpha(5)
    real(kdouble) :: beta(5)
    real(kdouble) :: q(5,5)
    real(kdouble) :: eig_val(5)
    real(kdouble) :: eig_mode(5,5)
    real(kdouble) :: norm, r_ans(5)

    call monolis_std_global_log_string("monolis_get_inverted_eigen_pair_from_tridiag")

    iter = 5

    n_get_eigen = 5

    alpha(1) = 2.0d0
    alpha(2) = 2.0d0
    alpha(3) = 2.0d0
    alpha(4) = 2.0d0
    alpha(5) = 2.0d0

    beta(1) = 1.0d0
    beta(2) = 1.0d0
    beta(3) = 1.0d0
    beta(4) = 1.0d0
    beta(5) = 0.0d0

    q = 0.0d0
    q(1,1) = 1.0d0
    q(2,2) = 1.0d0
    q(3,3) = 1.0d0
    q(4,4) = 1.0d0
    q(5,5) = 1.0d0

    call monolis_get_inverted_eigen_pair_from_tridiag(iter, n_get_eigen, &
      & alpha, beta, q, eig_val, eig_mode, norm)

    call monolis_test_check_eq_R1("monolis_get_eigen_pair_from_tridiag_test 1", eig_val(1), 0.267949192431122d0)
    call monolis_test_check_eq_R1("monolis_get_eigen_pair_from_tridiag_test 2", eig_val(2), 0.333333333333333d0)
    call monolis_test_check_eq_R1("monolis_get_eigen_pair_from_tridiag_test 3", eig_val(3), 0.5d0)
    call monolis_test_check_eq_R1("monolis_get_eigen_pair_from_tridiag_test 4", eig_val(4), 1.0d0)
    call monolis_test_check_eq_R1("monolis_get_eigen_pair_from_tridiag_test 5", eig_val(5), 3.732050807568877d0)

    r_ans(1) = 0.28867513459481281d0
    r_ans(2) =-0.5d0
    r_ans(3) = 0.57735026918962640d0
    r_ans(4) =-0.5d0
    r_ans(5) = 0.28867513459481270d0

    call monolis_test_check_eq_R("monolis_get_eigen_pair_from_tridiag_test 1 b", dabs(r_ans), dabs(eig_mode(:,5)))

    r_ans(1) =-0.5d0
    r_ans(2) = 0.5d0
    r_ans(3) = 0.0d0
    r_ans(4) =-0.5d0
    r_ans(5) = 0.5d0

    call monolis_test_check_eq_R("monolis_get_eigen_pair_from_tridiag_test 2 b", dabs(r_ans), dabs(eig_mode(:,4)))

    r_ans(1) = 0.5773502691896257d0
    r_ans(2) = 0.0d0
    r_ans(3) =-0.5773502691896257d0
    r_ans(4) = 0.0d0
    r_ans(5) = 0.5773502691896257d0

    call monolis_test_check_eq_R("monolis_get_eigen_pair_from_tridiag_test 3 b", dabs(r_ans), dabs(eig_mode(:,3)))

    r_ans(1) = 0.5d0
    r_ans(2) = 0.5d0
    r_ans(3) = 0.0d0
    r_ans(4) =-0.5d0
    r_ans(5) =-0.5d0

    call monolis_test_check_eq_R("monolis_get_eigen_pair_from_tridiag_test 4 b", dabs(r_ans), dabs(eig_mode(:,2)))

    r_ans(1) =-0.28867513459481281d0
    r_ans(2) =-0.5d0
    r_ans(3) =-0.57735026918962640d0
    r_ans(4) =-0.5d0
    r_ans(5) =-0.28867513459481270d0

    call monolis_test_check_eq_R("monolis_get_eigen_pair_from_tridiag_test 5 b", dabs(r_ans), dabs(eig_mode(:,1)))
  end subroutine monolis_get_eigen_pair_from_tridiag_test
end module mod_monolis_eigen_lanczos_util_test
