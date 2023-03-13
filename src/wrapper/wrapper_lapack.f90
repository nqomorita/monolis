!> LAPACK ラッパーモジュール
module mod_monolis_lapack
  use mod_monolis_utils

  implicit none

contains

  !> DSTEV 関数（実数型）
  subroutine monolis_lapack_dstev(n, D, S, eig_val, eig_mode)
    implicit none
    !> 行列の大きさ
    integer(kint), intent(in) :: n
    !> 行列の対角成分
    real(kdouble) :: D(:)
    !> 行列の副対角成分
    real(kdouble) :: S(:)
    !> 固有値
    real(kdouble) :: eig_val(:)
    !> 固有ベクトル
    real(kdouble) :: eig_mode(:,:)
    real(kdouble), allocatable :: r1(:)
    real(kdouble), allocatable :: r2(:)
    integer(kint) :: ldz, info

    ldz = n

    call monolis_alloc_R_1d(r1, n - 1)
    call monolis_alloc_R_1d(r2, 2*n)

    eig_val = D
    r1 = S(1:n - 1)

    call dstev("V", n, eig_val, r1, eig_mode, ldz, r2, info)
  end subroutine monolis_lapack_dstev

end module mod_monolis_lapack