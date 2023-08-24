!> 密行列積関数群
module mod_monolis_matalg_dense
  use mod_monolis_utils
  implicit none

contains

  !> @ingroup linalg
  !> 密行列ベクトル積（実数型、通信なし）
  subroutine monolis_dense_matvec_local_R(N, M, MAT, X, Y, tdemv)
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

    Y(1:N) = matmul(MAT(1:N,1:M), X(1:M))

    t2 = monolis_get_time()

    if(present(tdemv)) tdemv = tdemv + t2 - t1
  end subroutine monolis_dense_matvec_local_R

  !> @ingroup linalg
  !> 密行列密行列積（実数型、通信なし）
  subroutine monolis_dense_matmul_local_R(N1, N2, N3, MAT1, MAT2, Y, tdemv)
    implicit none
    !> [in] 行列サイズ N1
    integer(kint) :: N1
    !> [in] 行列サイズ N2
    integer(kint) :: N2
    !> [in] 行列サイズ N3
    integer(kint) :: N3
    !> [in] 入力行列（N x M）
    real(kdouble) :: MAT1(:,:)
    !> [in] 入力行列（M x ）
    real(kdouble) :: MAT2(:,:)
    !> [out] 出力行列（N）
    real(kdouble) :: Y(:,:)
    !> [in,out] 計算時間
    real(kdouble), optional :: tdemv
    real(kdouble) :: t1, t2

    t1 = monolis_get_time()

    Y(1:N1,1:N3) = matmul(MAT1(1:N1,1:N2), MAT2(1:N2,1:N3))

    t2 = monolis_get_time()

    if(present(tdemv)) tdemv = tdemv + t2 - t1
  end subroutine monolis_dense_matmul_local_R

end module mod_monolis_matalg_dense