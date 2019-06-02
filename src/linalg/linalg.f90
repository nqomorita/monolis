module mod_monolis_linalg
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_linalg_com
  implicit none

contains

  subroutine monolis_vec_copy_R(n, ndof, X, Y)
    implicit none
    integer(kind=kint) :: i, n, ndof
    real(kind=kdouble) :: X(:), Y(:)

    do i = 1, n * ndof
      Y(i) = X(i)
    enddo
  end subroutine monolis_vec_copy_R

  subroutine monolis_vec_AXPY(n, ndof, alpha, X, Y, Z)
    implicit none
    integer(kind=kint) :: i, n, ndof
    real(kind=kdouble) :: alpha
    real(kind=kdouble) :: X(:), Y(:), Z(:)

    do i = 1, n * ndof
      Z(i) = alpha*X(i) + Y(i)
    enddo
  end subroutine monolis_vec_AXPY

  subroutine monolis_inner_product_I(monoCOM, n, ndof, X, Y, sum, tcomm)
    implicit none
    type(monolis_com) :: monoCOM
    integer(kind=kint) :: i, n, ndof, sum
    integer(kind=kint) :: X(:), Y(:)
    real(kind=kdouble) :: t1, t2
    real(kind=kdouble), optional :: tcomm

    sum = 0
    do i = 1, n * ndof
      sum = sum + X(i)*Y(i)
    enddo

    t1 = monolis_wtime()
    call monolis_allreduce_I1(sum, monolis_sum, monoCOM%comm)
    t2 = monolis_wtime()
    if(present(tcomm)) tcomm = tcomm + t2 - t1
  end subroutine monolis_inner_product_I

  subroutine monolis_inner_product_R(monoCOM, n, ndof, X, Y, sum, tcomm)
    implicit none
    type(monolis_com) :: monoCOM
    integer(kind=kint) :: i, n, ndof
    real(kind=kdouble) :: X(:), Y(:)
    real(kind=kdouble) :: t1, t2, sum
    real(kind=kdouble), optional :: tcomm

    sum = 0.0d0
    do i = 1, n * ndof
      sum = sum + X(i)*Y(i)
    enddo

    t1 = monolis_wtime()
    call monolis_allreduce_R1(sum, monolis_sum, monoCOM%comm)
    t2 = monolis_wtime()
    if (present(tcomm)) tcomm = tcomm + t2 - t1
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