module mod_monolis_linalg
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_linalg_com
  implicit none

contains

  subroutine monolis_inner_product_I(monoCOM, monoMAT, ndof, X, Y, sum, tcomm)
    use mod_monolis_com
    use mod_monolis_mat
    use mod_monolis_linalg_com
    implicit none
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: i, ndof, sum
    integer(kind=kint) :: X(:), Y(:)
    real(kind=kdouble) :: t1, t2
    real(kind=kdouble), optional :: tcomm

    sum = 0
    do i = 1, monoMAT%N * ndof
      sum = sum + X(i)*Y(i)
    end do

    t1 = monolis_wtime()
    call monolis_allreduce_I1(sum, monolis_sum, monoCOM%comm)
    t2 = monolis_wtime()
    if(present(tcomm)) tcomm = tcomm + t2 - t1
  end subroutine monolis_inner_product_I

  subroutine monolis_inner_product_R(monoCOM, monoMAT, ndof, X, Y, sum, tcomm)
    use mod_monolis_com
    use mod_monolis_mat
    use mod_monolis_linalg_com
    implicit none
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: i, ndof
    real(kind=kdouble) :: X(:), Y(:)
    real(kind=kdouble) :: t1, t2, sum
    real(kind=kdouble), optional :: tcomm

    sum = 0.0d0
    do i = 1, monoMAT%N * ndof
      sum = sum + X(i)*Y(i)
    end do

    t1= monolis_wtime()
    call monolis_allreduce_R1(sum, monolis_sum, monoCOM%comm)
    t2= monolis_wtime()
    if (present(tcomm)) tcomm = tcomm + t2 - t1
  end subroutine monolis_inner_product_R

  subroutine monolis_inner_product_R_local(monoCOM, monoMAT, ndof, X, Y, sum)
    use mod_monolis_com
    use mod_monolis_mat
    implicit none
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: i, ndof
    real(kind=kdouble) :: X(:), Y(:)
    real(kind=kdouble) :: sum

    sum = 0.0d0
    do i = 1, monoMAT%N * ndof
      sum = sum + X(i)*Y(i)
    end do
  end subroutine monolis_inner_product_R_local

end module mod_monolis_linalg