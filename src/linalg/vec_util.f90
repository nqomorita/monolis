!> ベクトル演算関数群
module mod_monolis_vec_util
  use mod_monolis_utils
  use mod_monolis_inner_product
  implicit none

contains

  !> @ingroup linalg
  !> ベクトル配列コピー（整数型）
  subroutine monolis_vec_copy_I(m, X, Y)
    implicit none
    !> [in] 計算点数✕計算点が持つ自由度
    integer(kint), intent(in) :: m
    !> [in] ベクトル 1 (コピー元)
    integer(kint), intent(in) :: X(:)
    !> [out] ベクトル 2 (コピー先)
    integer(kint), intent(out) :: Y(:)
    integer(kint) :: i

    call monolis_std_debug_log_header("monolis_vec_copy_I")

!$omp parallel default(none) &
!$omp & shared(X, Y) &
!$omp & firstprivate(m) &
!$omp & private(i)
!$omp do
    do i = 1, m
      Y(i) = X(i)
    enddo
!$omp end do
!$omp end parallel
  end subroutine monolis_vec_copy_I

  !> @ingroup linalg
  !> ベクトル配列コピー（実数型）
  subroutine monolis_vec_copy_R(m, X, Y)
    implicit none
    !> [in] 計算点数✕計算点が持つ自由度
    integer(kint), intent(in) :: m
    !> [in] ベクトル 1 (コピー元)
    real(kdouble), intent(in) :: X(:)
    !> [out] ベクトル 2 (コピー先)
    real(kdouble), intent(out) :: Y(:)
    integer(kint) :: i

    call monolis_std_debug_log_header("monolis_vec_copy_R")

!$omp parallel default(none) &
!$omp & shared(X, Y) &
!$omp & firstprivate(m) &
!$omp & private(i)
!$omp do
    do i = 1, m
      Y(i) = X(i)
    enddo
!$omp end do
!$omp end parallel
  end subroutine monolis_vec_copy_R

  !> @ingroup linalg
  !> ベクトル配列コピー（複素数型）
  subroutine monolis_vec_copy_C(m, X, Y)
    implicit none
    !> [in] 計算点数✕計算点が持つ自由度
    integer(kint), intent(in) :: m
    !> [in] ベクトル 1 (コピー元)
    complex(kdouble), intent(in) :: X(:)
    !> [out] ベクトル 2 (コピー先)
    complex(kdouble), intent(out) :: Y(:)
    integer(kint) :: i

    call monolis_std_debug_log_header("monolis_vec_copy_C")

!$omp parallel default(none) &
!$omp & shared(X, Y) &
!$omp & firstprivate(m) &
!$omp & private(i)
!$omp do
    do i = 1, m
      Y(i) = X(i)
    enddo
!$omp end do
!$omp end parallel
  end subroutine monolis_vec_copy_C

  !> @ingroup linalg
  !> ベクトル和 $z = alpha * x + y$ （整数型）
  subroutine monolis_vec_AXPBY_I(m, alpha, X, beta, Y, Z)
    implicit none
    !> [in] 計算点数✕計算点が持つ自由度
    integer(kint), intent(in) :: m
    !> [in] 係数 1
    integer(kint), intent(in) :: alpha
    !> [in] ベクトル 1
    integer(kint), intent(in) :: X(:)
    !> [in] 係数 2
    integer(kint), intent(in) :: beta
    !> [in] ベクトル 2
    integer(kint), intent(in) :: Y(:)
    !> [out] 結果ベクトル
    integer(kint), intent(out) :: Z(:)
    integer(kint) :: i

    call monolis_std_debug_log_header("monolis_vec_AXPBY_I")

!$omp parallel default(none) &
!$omp & shared(X, Y, Z) &
!$omp & firstprivate(m, alpha) &
!$omp & private(i)
!$omp do
    do i = 1, m
      Z(i) = alpha*X(i) + beta*Y(i)
    enddo
!$omp end do
!$omp end parallel
  end subroutine monolis_vec_AXPBY_I

  !> @ingroup linalg
  !> ベクトル和 $z = alpha * x + y$ （実数型）
  subroutine monolis_vec_AXPBY_R(m, alpha, X, beta, Y, Z)
    implicit none
    !> [in] 計算点数✕計算点が持つ自由度
    integer(kint), intent(in) :: m
    !> [in] 係数 1
    real(kdouble), intent(in) :: alpha
    !> [in] ベクトル 1
    real(kdouble), intent(in) :: X(:)
    !> [in] 係数 2
    real(kdouble), intent(in) :: beta
    !> [in] ベクトル 2
    real(kdouble), intent(in) :: Y(:)
    !> [out] 結果ベクトル
    real(kdouble), intent(out) :: Z(:)
    integer(kint) :: i

    call monolis_std_debug_log_header("monolis_vec_AXPBY_R")

!$omp parallel default(none) &
!$omp & shared(X, Y, Z) &
!$omp & firstprivate(m, alpha) &
!$omp & private(i)
!$omp do
    do i = 1, m
      Z(i) = alpha*X(i) + beta*Y(i)
    enddo
!$omp end do
!$omp end parallel
  end subroutine monolis_vec_AXPBY_R

  !> @ingroup linalg
  !> ベクトル和 $z = alpha * x + y$ （複素数型）
  subroutine monolis_vec_AXPBY_C(m, alpha, X, beta, Y, Z)
    implicit none
    !> [in] 計算点数✕計算点が持つ自由度
    integer(kint), intent(in) :: m
    !> [in] 係数 1
    complex(kdouble), intent(in) :: alpha
    !> [in] ベクトル 1
    complex(kdouble), intent(in) :: X(:)
    !> [in] 係数 2
    complex(kdouble), intent(in) :: beta
    !> [in] ベクトル 2
    complex(kdouble), intent(in) :: Y(:)
    !> [out] 結果ベクトル
    complex(kdouble), intent(out) :: Z(:)
    integer(kint) :: i

    call monolis_std_debug_log_header("monolis_vec_AXPBY_C")

!$omp parallel default(none) &
!$omp & shared(X, Y, Z) &
!$omp & firstprivate(m, alpha) &
!$omp & private(i)
!$omp do
    do i = 1, m
      Z(i) = alpha*X(i) + beta*Y(i)
    enddo
!$omp end do
!$omp end parallel
  end subroutine monolis_vec_AXPBY_C

  !> @ingroup linalg
  !> Gram Schmidt 法（実数型）
  subroutine monolis_gram_schmidt_R(monoCOM, n_vec, m, p, q)
    implicit none
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] ベクトル本数
    integer(kint), intent(in) :: n_vec
    !> [in] 計算点数✕計算点が持つ自由度
    integer(kint), intent(in) :: m
    !> [in,out] 入力ベクトル
    real(kdouble), intent(inout) :: p(:)
    !> [in] 直交化されるベクトル
    real(kdouble), intent(in) :: q(:,:)
    integer(kint) :: i, j
    real(kdouble) :: norm, tdotp, tcomm_dotp

    do i = 1, n_vec
      call monolis_inner_product_main_R(monoCOM, m, p, q(:,i), norm, tdotp, tcomm_dotp)
      do j = 1, m
        p(j) = p(j) - norm*q(j,i)
      enddo
    enddo
  end subroutine monolis_gram_schmidt_R
end module mod_monolis_vec_util