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
    call monolis_lapack_dgeqrf_test()
    call monolis_lapack_dgeev_test()
    call monolis_lapack_dsysv_test()
    call monolis_lapack_dsysv_with_select_test()
  end subroutine monolis_lapack_test

  subroutine monolis_lapack_dense_matvec_local_R_test()
    implicit none
    real(kdouble) :: a(3), b(5), b_th(5), mat_dense(5,3), time

    call monolis_std_global_log_string("monolis_lapack_dense_matvec_local_R")

    mat_dense = 0.0d0
    mat_dense(1, 1) = 2.0d0
    mat_dense(1, 2) = 1.0d0
    mat_dense(2, 1) = 1.0d0
    mat_dense(2, 2) = 2.0d0
    mat_dense(2, 3) = 3.0d0
    mat_dense(3, 2) = 1.0d0
    mat_dense(3, 3) = 2.0d0
    mat_dense(4, 3) = 4.0d0
    mat_dense(5, 3) = 5.0d0

    a(1) = 1.0d0
    a(2) = 1.0d0
    a(3) = 1.0d0

    call monolis_lapack_dense_matvec_local_R(5, 3, mat_dense, a, b, time)

    b_th = matmul(mat_dense, a)

    call monolis_test_check_eq_R("monolis_lapack_dense_matvec_local_R", b, b_th)
  end subroutine monolis_lapack_dense_matvec_local_R_test

  subroutine monolis_lapack_dense_matvec_trans_local_R_test()
    implicit none
    real(kdouble) :: a(3), b(5), b_th(5), mat_dense(3,5), time

    call monolis_std_global_log_string("monolis_lapack_dense_matvec_trans_local_R")

    mat_dense = 0.0d0
    mat_dense(1, 1) = 2.0d0
    mat_dense(2, 1) = 1.0d0
    mat_dense(1, 2) = 1.0d0
    mat_dense(2, 2) = 2.0d0
    mat_dense(3, 2) = 3.0d0
    mat_dense(2, 3) = 1.0d0
    mat_dense(3, 3) = 2.0d0
    mat_dense(3, 4) = 4.0d0
    mat_dense(3, 5) = 5.0d0

    a(1) = 1.0d0
    a(2) = 1.0d0
    a(3) = 1.0d0

    call monolis_lapack_dense_matvec_trans_local_R(5, 3, mat_dense, a, b, time)

    a(1) = 1.0d0
    a(2) = 1.0d0
    a(3) = 1.0d0

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

  subroutine monolis_lapack_dgeqrf_test()
    implicit none
    !> 行列の大きさ
    integer(kint) :: m, n
    !> 入力行列
    real(kdouble) :: A(3,2)
    !> 直交行列 Q
    real(kdouble) :: Q(3,2)
    !> 上三角行列 R
    real(kdouble) :: R(2,2)
    !> 再構成行列 Q * R
    real(kdouble) :: QR(3,2)
    !> 直交性確認用 Q^T * Q
    real(kdouble) :: QtQ(2,2)
    real(kdouble) :: I_ans(2,2)

    call monolis_std_global_log_string("monolis_lapack_dgeqrf")

    m = 3
    n = 2

    A(1,1) = 1.0d0; A(1,2) = 1.0d0
    A(2,1) = 1.0d0; A(2,2) = 0.0d0
    A(3,1) = 0.0d0; A(3,2) = 1.0d0

    call monolis_lapack_dgeqrf(m, n, A, Q, R)

    !> R は上三角（下三角部分は 0）
    call monolis_test_check_eq_R1("monolis_lapack_dgeqrf_test R lower", R(2,1), 0.0d0)

    !> A = Q * R を満たす
    QR = matmul(Q, R)
    call monolis_test_check_eq_R("monolis_lapack_dgeqrf_test QR = A", &
      & reshape(QR, [6]), reshape(A, [6]))

    !> Q は直交（Q^T * Q = I）
    QtQ = matmul(transpose(Q), Q)
    I_ans = 0.0d0
    I_ans(1,1) = 1.0d0
    I_ans(2,2) = 1.0d0
    call monolis_test_check_eq_R("monolis_lapack_dgeqrf_test Q orthogonal", &
      & reshape(QtQ, [4]), reshape(I_ans, [4]))
  end subroutine monolis_lapack_dgeqrf_test

  subroutine monolis_lapack_dgeev_test()
    implicit none
    !> 行列の大きさ
    integer(kint) :: n, j
    !> 一般実行列
    real(kdouble) :: A(2,2)
    !> 固有値の実部・虚部
    real(kdouble) :: eig_re(2), eig_im(2)
    !> 右固有ベクトルの実部・虚部
    real(kdouble) :: eig_vec_re(2,2), eig_vec_im(2,2)
    real(kdouble) :: res_re(2), res_im(2), zero(2)

    call monolis_std_global_log_string("monolis_lapack_dgeev")

    n = 2
    zero = 0.0d0

    !> 実固有値を持つ対称行列 [[2,1],[1,2]] -> 固有値 1, 3
    A(1,1) = 2.0d0; A(1,2) = 1.0d0
    A(2,1) = 1.0d0; A(2,2) = 2.0d0

    call monolis_lapack_dgeev(n, A, eig_re, eig_im, eig_vec_re, eig_vec_im)

    !> 虚部は 0、トレース（固有値の和）= 4、行列式（固有値の積）= 3
    call monolis_test_check_eq_R("monolis_lapack_dgeev_test real eig_im", eig_im, zero)
    call monolis_test_check_eq_R1("monolis_lapack_dgeev_test real trace", &
      & eig_re(1) + eig_re(2), 4.0d0)
    call monolis_test_check_eq_R1("monolis_lapack_dgeev_test real det", &
      & eig_re(1)*eig_re(2), 3.0d0)

    !> 各固有対が A v = lambda v を満たす（順序に依存しない確認）
    do j = 1, n
      res_re = matmul(A, eig_vec_re(:,j)) - eig_re(j)*eig_vec_re(:,j) + eig_im(j)*eig_vec_im(:,j)
      res_im = matmul(A, eig_vec_im(:,j)) - eig_re(j)*eig_vec_im(:,j) - eig_im(j)*eig_vec_re(:,j)
      call monolis_test_check_eq_R("monolis_lapack_dgeev_test real res_re", res_re, zero)
      call monolis_test_check_eq_R("monolis_lapack_dgeev_test real res_im", res_im, zero)
    enddo

    !> 複素共役固有値を持つ回転行列 [[0,-1],[1,0]] -> 固有値 +i, -i
    A(1,1) = 0.0d0; A(1,2) =-1.0d0
    A(2,1) = 1.0d0; A(2,2) = 0.0d0

    call monolis_lapack_dgeev(n, A, eig_re, eig_im, eig_vec_re, eig_vec_im)

    !> 実部は 0、虚部は ±1（和は 0）
    call monolis_test_check_eq_R("monolis_lapack_dgeev_test cmplx eig_re", eig_re, zero)
    call monolis_test_check_eq_R1("monolis_lapack_dgeev_test cmplx eig_im abs 1", dabs(eig_im(1)), 1.0d0)
    call monolis_test_check_eq_R1("monolis_lapack_dgeev_test cmplx eig_im abs 2", dabs(eig_im(2)), 1.0d0)
    call monolis_test_check_eq_R1("monolis_lapack_dgeev_test cmplx eig_im sum", &
      & eig_im(1) + eig_im(2), 0.0d0)

    !> 各固有対が A v = lambda v を満たす
    do j = 1, n
      res_re = matmul(A, eig_vec_re(:,j)) - eig_re(j)*eig_vec_re(:,j) + eig_im(j)*eig_vec_im(:,j)
      res_im = matmul(A, eig_vec_im(:,j)) - eig_re(j)*eig_vec_im(:,j) - eig_im(j)*eig_vec_re(:,j)
      call monolis_test_check_eq_R("monolis_lapack_dgeev_test cmplx res_re", res_re, zero)
      call monolis_test_check_eq_R("monolis_lapack_dgeev_test cmplx res_im", res_im, zero)
    enddo
  end subroutine monolis_lapack_dgeev_test

  subroutine monolis_lapack_dsysv_test()
    implicit none
    !> 行列の大きさ
    integer(kint) :: n
    !> 対称行列
    real(kdouble) :: A(3,3)
    !> 右辺ベクトル
    real(kdouble) :: b(3)
    !> 解ベクトル
    real(kdouble) :: x(3)
    real(kdouble) :: x_ans(3)

    call monolis_std_global_log_string("monolis_lapack_dsysv")

    n = 3

    !> 三重対角対称行列 [[2,-1,0],[-1,2,-1],[0,-1,2]]
    A = 0.0d0
    A(1,1) = 2.0d0; A(1,2) =-1.0d0
    A(2,1) =-1.0d0; A(2,2) = 2.0d0; A(2,3) =-1.0d0
    A(3,2) =-1.0d0; A(3,3) = 2.0d0

    b(1) = 1.0d0
    b(2) = 0.0d0
    b(3) = 1.0d0

    call monolis_lapack_dsysv(n, A, b, x)

    x_ans(1) = 1.0d0
    x_ans(2) = 1.0d0
    x_ans(3) = 1.0d0

    call monolis_test_check_eq_R("monolis_lapack_dsysv_test", x, x_ans)
  end subroutine monolis_lapack_dsysv_test

  subroutine monolis_lapack_dsysv_with_select_test()
    implicit none
    !> 行列の大きさ
    integer(kint) :: n
    !> 対称行列
    real(kdouble) :: A(3,3)
    !> 右辺ベクトル
    real(kdouble) :: b(3)
    !> 解ベクトル
    real(kdouble) :: x(3)
    real(kdouble) :: x_ans(3)
    !> 選択フラグ
    logical :: P(3)

    call monolis_std_global_log_string("monolis_lapack_dsysv_with_select")

    n = 3

    !> 三重対角対称行列 [[2,-1,0],[-1,2,-1],[0,-1,2]]
    A = 0.0d0
    A(1,1) = 2.0d0; A(1,2) =-1.0d0
    A(2,1) =-1.0d0; A(2,2) = 2.0d0; A(2,3) =-1.0d0
    A(3,2) =-1.0d0; A(3,3) = 2.0d0

    b(1) = 1.0d0
    b(2) = 0.0d0
    b(3) = 1.0d0

    !> 1, 3 番目の自由度のみを選択（部分行列 [[2,0],[0,2]], b = [1,1]）
    P(1) = .true.
    P(2) = .false.
    P(3) = .true.

    call monolis_lapack_dsysv_with_select(n, A, b, x, P)

    x_ans(1) = 0.5d0
    x_ans(2) = 0.0d0
    x_ans(3) = 0.5d0

    call monolis_test_check_eq_R("monolis_lapack_dsysv_with_select_test", x, x_ans)
  end subroutine monolis_lapack_dsysv_with_select_test

end module mod_monolis_lapack_test
