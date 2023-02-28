module mod_monolis_vec_util
  use mod_monolis_utils
  implicit none

contains

  !> @ingroup dev_linalg
  subroutine monolis_vec_copy_R(n, ndof, X, Y)
    implicit none
    integer(kint) :: i, n, ndof
    real(kdouble) :: X(:), Y(:)

#ifdef DEBUG
    call monolis_std_debug_log_header("monolis_vec_copy_R")
#endif

!$omp parallel default(none) &
!$omp & shared(X, Y) &
!$omp & firstprivate(n, ndof) &
!$omp & private(i)
!$omp do
    do i = 1, n * ndof
      Y(i) = X(i)
    enddo
!$omp end do
!$omp end parallel
  end subroutine monolis_vec_copy_R

  !> @ingroup dev_linalg
  subroutine monolis_vec_AXPY(n, ndof, alpha, X, Y, Z)
    implicit none
    integer(kint) :: i, n, ndof
    real(kdouble) :: alpha
    real(kdouble) :: X(:), Y(:), Z(:)

#ifdef DEBUG
    call monolis_std_debug_log_header("monolis_vec_AXPY")
#endif

!$omp parallel default(none) &
!$omp & shared(X, Y, Z) &
!$omp & firstprivate(n, ndof, alpha) &
!$omp & private(i)
!$omp do
    do i = 1, n * ndof
      Z(i) = alpha*X(i) + Y(i)
    enddo
!$omp end do
!$omp end parallel
  end subroutine monolis_vec_AXPY

end module mod_monolis_vec_util