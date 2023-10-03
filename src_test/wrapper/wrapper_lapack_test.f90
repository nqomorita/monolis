!> LAPACK ラッパーモジュール
module mod_monolis_lapack_test
  use mod_monolis_utils
  use mod_monolis_lapack

  implicit none

contains

  subroutine monolis_lapack_test()
    implicit none

    call monolis_lapack_dense_matvec_local_R_test()
    call monolis_lapack_dense_matvec_trans_local_R_test()
    call monolis_lapack_stev_R_test()
    call monolis_lapack_LU_fact_R_test()
  end subroutine monolis_lapack_test

  subroutine monolis_lapack_dense_matvec_local_R_test()
    implicit none
    real(kdouble) :: a(5), b(5), b_th(5), mat_dense(5,5), time

    call monolis_std_global_log_string("monolis_lapack_dense_matvec_local_R")

    mat_dense(1, 1) = 2.0d0
    mat_dense(1, 2) = 1.0d0
    mat_dense(2, 1) = 1.0d0
    mat_dense(2, 2) = 2.0d0
    mat_dense(2, 3) = 3.0d0
    mat_dense(3, 2) = 1.0d0
    mat_dense(3, 3) = 2.0d0
    mat_dense(3, 4) = 4.0d0
    mat_dense(4, 3) = 1.0d0
    mat_dense(4, 4) = 2.0d0
    mat_dense(4, 5) = 5.0d0
    mat_dense(5, 4) = 1.0d0
    mat_dense(5, 5) = 2.0d0

    a(1) = 1.0d0
    a(2) = 1.0d0
    a(3) = 1.0d0
    a(4) = 1.0d0
    a(5) = 1.0d0

    call monolis_lapack_dense_matvec_local_R(5, 5, mat_dense, a, b, time)

    b_th = matmul(mat_dense, a)

    call monolis_test_check_eq_R("monolis_lapack_dense_matvec_local_R", b, b_th)
  end subroutine monolis_lapack_dense_matvec_local_R_test

  subroutine monolis_lapack_dense_matvec_trans_local_R_test()
    implicit none
    real(kdouble) :: a(5), b(5), b_th(5), mat_dense(5,5), time

    call monolis_std_global_log_string("monolis_lapack_dense_matvec_trans_local_R")

    mat_dense(1, 1) = 2.0d0
    mat_dense(1, 2) = 1.0d0
    mat_dense(2, 1) = 1.0d0
    mat_dense(2, 2) = 2.0d0
    mat_dense(2, 3) = 3.0d0
    mat_dense(3, 2) = 1.0d0
    mat_dense(3, 3) = 2.0d0
    mat_dense(3, 4) = 4.0d0
    mat_dense(4, 3) = 1.0d0
    mat_dense(4, 4) = 2.0d0
    mat_dense(4, 5) = 5.0d0
    mat_dense(5, 4) = 1.0d0
    mat_dense(5, 5) = 2.0d0

    a(1) = 1.0d0
    a(2) = 1.0d0
    a(3) = 1.0d0
    a(4) = 1.0d0
    a(5) = 1.0d0

    call monolis_lapack_dense_matvec_trans_local_R(5, 5, mat_dense, a, b, time)

    b_th = matmul(transpose(mat_dense), a)

    call monolis_test_check_eq_R("monolis_lapack_dense_matvec_trans_local_R", b, b_th)
  end subroutine monolis_lapack_dense_matvec_trans_local_R_test

  subroutine monolis_lapack_LU_fact_R_test()
    implicit none
    !> 行列の大きさ
    integer(kint) :: N = 3
    integer(kint) :: IPV(3)
    real(kdouble) :: a(3,3), c(3,3), B(3), X(3), X_th(3)

    a(1,1) = 3.0d0; a(1,2) = 3.0d0; a(1,3) = 2.0d0
    a(2,1) = 2.0d0; a(2,2) = 3.0d0; a(2,3) = 3.0d0
    a(3,1) = 1.0d0; a(3,2) = 2.0d0; a(3,3) = 1.0d0

    call monolis_lapack_LU_fact_R(N, A, IPV)

    B(1) = 1.0d0
    B(2) = 2.0d0
    B(3) = 3.0d0

    call monolis_lapack_LU_solve_R(N, A, B, IPV, X)

    c(1,1) = 0.75d0; c(1,2) =-0.25d0; c(1,3) =-0.75d0
    c(2,1) =-0.25d0; c(2,2) =-0.25d0; c(2,3) = 1.25d0
    c(3,1) =-0.25d0; c(3,2) = 0.75d0; c(3,3) =-0.75d0

    X_th = matmul(c, B)

    call monolis_test_check_eq_R("monolis_lapack_LU_fact_R_test", X, X_th)
  end subroutine monolis_lapack_LU_fact_R_test

  subroutine monolis_lapack_stev_R_test()
    implicit none
    !> 行列の大きさ
    integer(kint) :: N
    !> 行列の対角成分
    real(kdouble) :: D(5)
    !> 行列の副対角成分
    real(kdouble) :: S(5)
    !> 固有値
    real(kdouble) :: eig_val(5)
    !> 固有ベクトル
    real(kdouble) :: eig_mode(5,5)
    real(kdouble) :: r_ans(5)

    call monolis_std_global_log_string("monolis_lapack_stev_R")

    N = 5

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

    call monolis_lapack_stev_R(N, D, S, eig_val, eig_mode)

    call monolis_test_check_eq_R1("monolis_lapack_stev_R_test 1", eig_val(1), 0.267949192431122d0)
    call monolis_test_check_eq_R1("monolis_lapack_stev_R_test 2", eig_val(2), 1.0d0)
    call monolis_test_check_eq_R1("monolis_lapack_stev_R_test 3", eig_val(3), 2.0d0)
    call monolis_test_check_eq_R1("monolis_lapack_stev_R_test 4", eig_val(4), 3.0d0)
    call monolis_test_check_eq_R1("monolis_lapack_stev_R_test 5", eig_val(5), 3.732050807568877d0)

    r_ans(1) = 0.28867513459481281d0
    r_ans(2) =-0.5d0
    r_ans(3) = 0.57735026918962640d0
    r_ans(4) =-0.5d0
    r_ans(5) = 0.28867513459481270d0

    call monolis_test_check_eq_R("monolis_lapack_stev_R_test 1 b", dabs(r_ans), dabs(eig_mode(:,1)))

    r_ans(1) =-0.5d0
    r_ans(2) = 0.5d0
    r_ans(3) = 0.0d0
    r_ans(4) =-0.5d0
    r_ans(5) = 0.5d0

    call monolis_test_check_eq_R("monolis_lapack_stev_R_test 2 b", dabs(r_ans), dabs(eig_mode(:,2)))

    r_ans(1) = 0.5773502691896257d0
    r_ans(2) = 0.0d0
    r_ans(3) =-0.5773502691896257d0
    r_ans(4) = 0.0d0
    r_ans(5) = 0.5773502691896257d0

    call monolis_test_check_eq_R("monolis_lapack_stev_R_test 3 b", dabs(r_ans), dabs(eig_mode(:,3)))

    r_ans(1) = 0.5d0
    r_ans(2) = 0.5d0
    r_ans(3) = 0.0d0
    r_ans(4) =-0.5d0
    r_ans(5) =-0.5d0

    call monolis_test_check_eq_R("monolis_lapack_stev_R_test 4 b", dabs(r_ans), dabs(eig_mode(:,4)))

    r_ans(1) =-0.28867513459481281d0
    r_ans(2) =-0.5d0
    r_ans(3) =-0.57735026918962640d0
    r_ans(4) =-0.5d0
    r_ans(5) =-0.28867513459481270d0

    call monolis_test_check_eq_R("monolis_lapack_stev_R_test 5 b", dabs(r_ans), dabs(eig_mode(:,5)))
  end subroutine monolis_lapack_stev_R_test

end module mod_monolis_lapack_test
