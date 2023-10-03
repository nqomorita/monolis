!> LAPACK ラッパーモジュール
module mod_monolis_lapack
  use mod_monolis_utils

  implicit none

contains

  !> @ingroup wrapper
  !> 密行列ベクトル積（DGEMV 関数、実数型）
  subroutine monolis_lapack_dense_matvec_local_R(N, M, MAT, X, Y, tdemv)
    implicit none
    !> [in] 行列サイズ N
    integer(kint) :: N
    !> [in] 行列サイズ M
    integer(kint) :: M
    !> [in] 入力行列（N x M）
    real(kdouble) :: MAT(:,:)
    !> [in] 入力ベクトル（M）
    real(kdouble) :: X(:)
    !> [out] 出力ベクトル（N）
    real(kdouble) :: Y(:)
    !> [in,out] 計算時間
    real(kdouble), optional :: tdemv
    real(kdouble) :: t1, t2

    t1 = monolis_get_time()

    Y = 0.0d0
    call dgemv("N", N, M, 1.0d0, MAT, N, X, 1, 1.0d0, Y, 1)

    t2 = monolis_get_time()

    if(present(tdemv)) tdemv = tdemv + t2 - t1
  end subroutine monolis_lapack_dense_matvec_local_R

  !> @ingroup wrapper
  !> 密行列ベクトル積（DGEMV 関数、実数型、転置）
  subroutine monolis_lapack_dense_matvec_trans_local_R(N, M, MAT, X, Y, tdemv)
    implicit none
    !> [in] 行列サイズ N
    integer(kint) :: N
    !> [in] 行列サイズ M
    integer(kint) :: M
    !> [in] 入力行列（N x M）
    real(kdouble) :: MAT(:,:)
    !> [in] 入力ベクトル（M）
    real(kdouble) :: X(:)
    !> [out] 出力ベクトル（N）
    real(kdouble) :: Y(:)
    !> [in,out] 計算時間
    real(kdouble), optional :: tdemv
    real(kdouble) :: t1, t2

    t1 = monolis_get_time()

    Y = 0.0d0
    call dgemv("T", N, M, 1.0d0, MAT, N, X, 1, 1.0d0, Y, 1)

    t2 = monolis_get_time()

    if(present(tdemv)) tdemv = tdemv + t2 - t1
  end subroutine monolis_lapack_dense_matvec_trans_local_R

  !> @ingroup wrapper
  !> 三重対角行列の固有値問題求解（DSTEV 関数、実数型）
  subroutine monolis_lapack_stev_R(N, D, S, eig_val, eig_mode)
    implicit none
    !> [in] 行列の大きさ
    integer(kint), intent(in) :: N
    !> [in] 行列の対角成分
    real(kdouble), intent(in) :: D(:)
    !> [in] 行列の副対角成分
    real(kdouble), intent(in) :: S(:)
    !> [out] 固有値
    real(kdouble), intent(out) :: eig_val(:)
    !> [out] 固有ベクトル
    real(kdouble), intent(out) :: eig_mode(:,:)
    real(kdouble), allocatable :: r1(:)
    real(kdouble), allocatable :: r2(:)
    integer(kint) :: ldz, info

    ldz = N

    call monolis_alloc_R_1d(r1, N - 1)
    call monolis_alloc_R_1d(r2, 2*N)

    eig_val(1:N) = D(1:N)
    r1 = S(1:N - 1)

    call dstev("V", N, eig_val, r1, eig_mode, ldz, r2, info)
  end subroutine monolis_lapack_stev_R

  !> @ingroup wrapper
  !> 密行列の LU 分解・分解フェーズ（DGETRF 関数、実数型）
  subroutine monolis_lapack_LU_fact_R(N, A, IPV)
    implicit none
    !> 行列の大きさ
    integer(4) :: N
    !> 行列
    real(8) :: A(N,N)
    !> ピボット情報配列
    integer(4) :: IPV(N)
    integer(4) :: LDA, info

    LDA = N
    call dgetrf(N, N, A, LDA, IPV, info)
  end subroutine monolis_lapack_LU_fact_R

  !> @ingroup wrapper
  !> 密行列の LU 分解・求解フェーズ（DGETRS 関数、実数型）
  subroutine monolis_lapack_LU_solve_R(N, A, B, IPV, X)
    implicit none
    !> 行列の大きさ
    integer(4) :: N
    !> 行列
    real(8) :: A(N,N)
    !> 右辺ベクトル
    real(8) :: B(N)
    !> ピボット情報配列
    integer(4) :: IPV(N)
    !> 解ベクトル
    real(8) :: X(N)
    integer(4) :: NRHS, LDA, info

    LDA = N
    NRHS = 1
    X = B
    call dgetrs("N", N, NRHS, A, LDA, IPV, X, LDA, info)
  end subroutine monolis_lapack_LU_solve_R
end module mod_monolis_lapack
