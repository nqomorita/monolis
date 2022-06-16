module mod_monolis_linalg
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_linalg_com
  use mod_monolis_util
  implicit none

  private
  public :: monolis_vec_copy_R
  public :: monolis_vec_AXPY
  public :: monolis_inner_product_I
  public :: monolis_inner_product_R
  public :: monolis_inner_product_R_local

contains

  subroutine monolis_vec_copy_R(n, ndof, X, Y)
    implicit none
    integer(kint) :: i, n, ndof
    real(kdouble) :: X(:), Y(:)

#ifdef DEBUG
    call monolis_debug_header("monolis_vec_copy_R")
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

  subroutine monolis_vec_AXPY(n, ndof, alpha, X, Y, Z)
    implicit none
    integer(kint) :: i, n, ndof
    real(kdouble) :: alpha
    real(kdouble) :: X(:), Y(:), Z(:)

#ifdef DEBUG
    call monolis_debug_header("monolis_vec_AXPY")
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

  subroutine monolis_inner_product_I(monoCOM, n, ndof, X, Y, sum, tdotp, tcomm)
    implicit none
    type(monolis_com) :: monoCOM
    integer(kint) :: i, n, ndof, sum
    integer(kint) :: X(:), Y(:)
    real(kdouble) :: t1, t2, t3
    real(kdouble), optional :: tdotp, tcomm

#ifdef DEBUG
    call monolis_debug_header("monolis_inner_product_I")
#endif

    t1 = monolis_get_time()
    sum = 0
!$omp parallel default(none) &
!$omp & shared(X, Y, sum) &
!$omp & firstprivate(n, ndof) &
!$omp & private(i)
!$omp do reduction(+:sum)
    do i = 1, n * ndof
      sum = sum + X(i)*Y(i)
    enddo
!$omp end do
!$omp end parallel

    t2 = monolis_get_time()
    call monolis_allreduce_I1(sum, monolis_sum, monoCOM%comm)
    t3 = monolis_get_time()

    if(present(tdotp)) tdotp = tdotp + t3 - t1
    if(present(tcomm)) tcomm = tcomm + t3 - t2
  end subroutine monolis_inner_product_I

  subroutine monolis_inner_product_R(monoCOM, n, ndof, X, Y, sum, tdotp, tcomm)
    implicit none
    type(monolis_com) :: monoCOM
    integer(kint) :: i, n, ndof
    real(kdouble) :: X(:), Y(:)
    real(kdouble) :: t1, t2, t3, sum
    real(kdouble), optional :: tdotp, tcomm

#ifdef DEBUG
    call monolis_debug_header("monolis_inner_product_R")
#endif

    t1 = monolis_get_time()
    sum = 0.0d0
!$omp parallel default(none) &
!$omp & shared(X, Y, sum) &
!$omp & firstprivate(n, ndof) &
!$omp & private(i)
!$omp do reduction(+:sum)
    do i = 1, n * ndof
      sum = sum + X(i)*Y(i)
    enddo
!$omp end do
!$omp end parallel

    t2 = monolis_get_time()
    call monolis_allreduce_R1(sum, monolis_sum, monoCOM%comm)
    t3 = monolis_get_time()

    if(present(tdotp)) tdotp = tdotp + t3 - t1
    if(present(tcomm)) tcomm = tcomm + t3 - t2
  end subroutine monolis_inner_product_R

  subroutine monolis_inner_product_R_local(monoCOM, n, ndof, X, Y, sum)
    implicit none
    type(monolis_com) :: monoCOM
    integer(kint) :: i, n, ndof
    real(kdouble) :: X(:), Y(:)
    real(kdouble) :: sum

    sum = 0.0d0
!$omp parallel default(none) &
!$omp & shared(X, Y, sum) &
!$omp & firstprivate(n, ndof) &
!$omp & private(i)
!$omp do reduction(+:sum)
    do i = 1, n * ndof
      sum = sum + X(i)*Y(i)
    enddo
!$omp end do
!$omp end parallel
  end subroutine monolis_inner_product_R_local
end module mod_monolis_linalg