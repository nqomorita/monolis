!> ベクトル演算関数群
module mod_monolis_vec_util
  use mod_monolis_utils
  implicit none

contains

  !> @ingroup dev_linalg
  !> ベクトル配列コピー（整数型）
  subroutine monolis_vec_copy_I(n, n_dof, X, Y)
    implicit none
    !> 自由度数
    integer(kint) :: n
    !> ブロックサイズ
    integer(kint) :: n_dof
    !> ベクトル 1
    integer(kint) :: X(:)
    !> ベクトル 2
    integer(kint) :: Y(:)
    integer(kint) :: i

    call monolis_std_debug_log_header("monolis_vec_copy_I")

!$omp parallel default(none) &
!$omp & shared(X, Y) &
!$omp & firstprivate(n, n_dof) &
!$omp & private(i)
!$omp do
    do i = 1, n * n_dof
      Y(i) = X(i)
    enddo
!$omp end do
!$omp end parallel
  end subroutine monolis_vec_copy_I

  !> @ingroup dev_linalg
  !> ベクトル配列コピー（実数型）
  subroutine monolis_vec_copy_R(n, n_dof, X, Y)
    implicit none
    !> 自由度数
    integer(kint) :: n
    !> ブロックサイズ
    integer(kint) :: n_dof
    !> ベクトル 1
    real(kdouble) :: X(:)
    !> ベクトル 2
    real(kdouble) :: Y(:)
    integer(kint) :: i

    call monolis_std_debug_log_header("monolis_vec_copy_R")

!$omp parallel default(none) &
!$omp & shared(X, Y) &
!$omp & firstprivate(n, n_dof) &
!$omp & private(i)
!$omp do
    do i = 1, n * n_dof
      Y(i) = X(i)
    enddo
!$omp end do
!$omp end parallel
  end subroutine monolis_vec_copy_R

  !> @ingroup dev_linalg
  !> ベクトル配列コピー（複素数型）
  subroutine monolis_vec_copy_C(n, n_dof, X, Y)
    implicit none
    !> 自由度数
    integer(kint) :: n
    !> ブロックサイズ
    integer(kint) :: n_dof
    !> ベクトル 1
    complex(kdouble) :: X(:)
    !> ベクトル 2
    complex(kdouble) :: Y(:)
    integer(kint) :: i

    call monolis_std_debug_log_header("monolis_vec_copy_C")

!$omp parallel default(none) &
!$omp & shared(X, Y) &
!$omp & firstprivate(n, n_dof) &
!$omp & private(i)
!$omp do
    do i = 1, n * n_dof
      Y(i) = X(i)
    enddo
!$omp end do
!$omp end parallel
  end subroutine monolis_vec_copy_C

  !> @ingroup dev_linalg
  !> ベクトル和 $z = \alpha * x + y$ （整数型）
  subroutine monolis_vec_AXPBY_I(n, n_dof, alpha, X, beta, Y, Z)
    implicit none
    !> 自由度数
    integer(kint) :: n
    !> ブロックサイズ
    integer(kint) :: n_dof
    !> 係数 1
    integer(kint) :: alpha
    !> ベクトル 1
    integer(kint) :: X(:)
    !> 係数 2
    integer(kint) :: beta
    !> ベクトル 2
    integer(kint) :: Y(:)
    !> 結果ベクトル
    integer(kint) :: Z(:)
    integer(kint) :: i

    call monolis_std_debug_log_header("monolis_vec_AXPBY_I")

!$omp parallel default(none) &
!$omp & shared(X, Y, Z) &
!$omp & firstprivate(n, n_dof, alpha) &
!$omp & private(i)
!$omp do
    do i = 1, n * n_dof
      Z(i) = alpha*X(i) + beta*Y(i)
    enddo
!$omp end do
!$omp end parallel
  end subroutine monolis_vec_AXPBY_I

  !> @ingroup dev_linalg
  !> ベクトル和 $z = \alpha * x + y$ （実数型）
  subroutine monolis_vec_AXPBY_R(n, n_dof, alpha, X, beta, Y, Z)
    implicit none
    !> 自由度数
    integer(kint) :: n
    !> ブロックサイズ
    integer(kint) :: n_dof
    !> 係数 1
    real(kdouble) :: alpha
    !> ベクトル 1
    real(kdouble) :: X(:)
    !> 係数 2
    real(kdouble) :: beta
    !> ベクトル 2
    real(kdouble) :: Y(:)
    !> 結果ベクトル
    real(kdouble) :: Z(:)
    integer(kint) :: i

    call monolis_std_debug_log_header("monolis_vec_AXPBY_R")

!$omp parallel default(none) &
!$omp & shared(X, Y, Z) &
!$omp & firstprivate(n, n_dof, alpha) &
!$omp & private(i)
!$omp do
    do i = 1, n * n_dof
      Z(i) = alpha*X(i) + beta*Y(i)
    enddo
!$omp end do
!$omp end parallel
  end subroutine monolis_vec_AXPBY_R

  !> @ingroup dev_linalg
  !> ベクトル和 $z = \alpha * x + y$ （複素数型）
  subroutine monolis_vec_AXPBY_C(n, n_dof, alpha, X, beta, Y, Z)
    implicit none
    !> 自由度数
    integer(kint) :: n
    !> ブロックサイズ
    integer(kint) :: n_dof
    !> 係数 1
    complex(kdouble) :: alpha
    !> ベクトル 1
    complex(kdouble) :: X(:)
    !> 係数 2
    complex(kdouble) :: beta
    !> ベクトル 2
    complex(kdouble) :: Y(:)
    !> 結果ベクトル
    complex(kdouble) :: Z(:)
    integer(kint) :: i

    call monolis_std_debug_log_header("monolis_vec_AXPBY_C")

!$omp parallel default(none) &
!$omp & shared(X, Y, Z) &
!$omp & firstprivate(n, n_dof, alpha) &
!$omp & private(i)
!$omp do
    do i = 1, n * n_dof
      Z(i) = alpha*X(i) + beta*Y(i)
    enddo
!$omp end do
!$omp end parallel
  end subroutine monolis_vec_AXPBY_C
end module mod_monolis_vec_util