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
    integer(kind=kint) :: i, n, ndof
    real(kind=kdouble) :: X(:), Y(:)

#ifdef DEBUG
    call monolis_debug_header("monolis_vec_copy_R")
#endif

    do i = 1, n * ndof
      Y(i) = X(i)
    enddo
  end subroutine monolis_vec_copy_R

  subroutine monolis_vec_AXPY(n, ndof, alpha, X, Y, Z)
    implicit none
    integer(kind=kint) :: i, n, ndof
    real(kind=kdouble) :: alpha
    real(kind=kdouble) :: X(:), Y(:), Z(:)

#ifdef DEBUG
    call monolis_debug_header("monolis_vec_AXPY")
#endif

    do i = 1, n * ndof
      Z(i) = alpha*X(i) + Y(i)
    enddo
  end subroutine monolis_vec_AXPY

  subroutine monolis_inner_product_I(monoCOM, n, ndof, X, Y, sum, tdotp, tcomm)
    implicit none
    type(monolis_com) :: monoCOM
    integer(kind=kint) :: i, n, ndof, sum
    integer(kind=kint) :: X(:), Y(:)
    real(kind=kdouble) :: t1, t2, t3
    real(kind=kdouble) :: tdotp, tcomm

#ifdef DEBUG
    call monolis_debug_header("monolis_inner_product_I")
#endif

    t1 = monolis_get_time()
    sum = 0
    do i = 1, n * ndof
      sum = sum + X(i)*Y(i)
    enddo

    t2 = monolis_get_time()
    call monolis_allreduce_I1(sum, monolis_sum, monoCOM%comm)
    t3 = monolis_get_time()

    tdotp = tdotp + t3 - t1
    tcomm = tcomm + t3 - t2
  end subroutine monolis_inner_product_I

  subroutine monolis_inner_product_R(monoCOM, n, ndof, X, Y, sum, tdotp, tcomm)
    implicit none
    type(monolis_com) :: monoCOM
    integer(kind=kint) :: i, n, ndof
    real(kind=kdouble) :: X(:), Y(:)
    real(kind=kdouble) :: t1, t2, t3, sum
    real(kind=kdouble) :: tdotp, tcomm

#ifdef DEBUG
    call monolis_debug_header("monolis_inner_product_R")
#endif

    t1 = monolis_get_time()
    sum = 0.0d0
    do i = 1, n * ndof
      sum = sum + X(i)*Y(i)
    enddo

    t2 = monolis_get_time()
    call monolis_allreduce_R1(sum, monolis_sum, monoCOM%comm)
    t3 = monolis_get_time()

    tdotp = tdotp + t3 - t1
    tcomm = tcomm + t3 - t2
  end subroutine monolis_inner_product_R

  subroutine monolis_inner_product_R_local(monoCOM, n, ndof, X, Y, sum)
    implicit none
    type(monolis_com) :: monoCOM
    integer(kind=kint) :: i, n, ndof
    real(kind=kdouble) :: X(:), Y(:)
    real(kind=kdouble) :: sum

    sum = 0.0d0
    do i = 1, n * ndof
      sum = sum + X(i)*Y(i)
    enddo
  end subroutine monolis_inner_product_R_local
end module mod_monolis_linalg