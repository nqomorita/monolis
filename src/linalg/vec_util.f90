!> ベクトル演算関数群
module mod_monolis_vec_util
  use mod_monolis_utils
  use mod_monolis_inner_product
  implicit none

contains

  !> @ingroup linalg
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

  !> @ingroup linalg
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

  !> @ingroup linalg
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

  !> @ingroup linalg
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

  !> @ingroup linalg
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

  !> @ingroup linalg
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

  !> @ingroup linalg
  !> Gram Schmidt 法（実数型）
  subroutine monolis_gram_schmidt_R(monoCOM, n_vec, N, NDOF, p, q)
    implicit none
    !> 通信テーブル構造体
    type(monolis_com) :: monoCOM
    !> ベクトル本数
    integer(kint) :: n_vec
    !> 計算点数
    integer(kint) :: N
    !> 計算点上の自由度
    integer(kint) :: NDOF
    !> 入力ベクトル
    real(kdouble) :: p(:)
    !> 直交化されるベクトル
    real(kdouble) :: q(:,:)
    integer(kint) :: i, j
    real(kdouble) :: norm, tdotp, tcomm_dotp

    do i = 1, n_vec
      call monolis_inner_product_main_R(monoCOM, N, NDOF, p, q(:,i), norm, tdotp, tcomm_dotp)
      do j = 1, N*NDOF
        p(j) = p(j) - norm*q(j,i)
      enddo
    enddo
  end subroutine monolis_gram_schmidt_R
end module mod_monolis_vec_util