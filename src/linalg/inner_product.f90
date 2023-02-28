module mod_monolis_linalg
  use mod_monolis_utils
  implicit none

contains

  subroutine monolis_inner_product_main_I(monoCOM, n, ndof, X, Y, sum, tdotp, tcomm)
    implicit none
    type(monolis_com) :: monoCOM
    integer(kint) :: i, n, ndof, sum
    integer(kint) :: X(:), Y(:)
    real(kdouble) :: t1, t2, t3
    real(kdouble), optional :: tdotp, tcomm

#ifdef DEBUG
    call monolis_std_debug_log_header("monolis_inner_product_main_I")
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
    call monolis_allreduce_I1(sum, monolis_mpi_sum, monoCOM%comm)
    t3 = monolis_get_time()

    if(present(tdotp)) tdotp = tdotp + t3 - t1
    if(present(tcomm)) tcomm = tcomm + t3 - t2
  end subroutine monolis_inner_product_main_I

  subroutine monolis_inner_product_main_R(monoCOM, n, ndof, X, Y, sum, tdotp, tcomm)
    implicit none
    type(monolis_com) :: monoCOM
    integer(kint) :: i, n, ndof
    real(kdouble) :: X(:), Y(:)
    real(kdouble) :: t1, t2, t3, sum
    real(kdouble), optional :: tdotp, tcomm

#ifdef DEBUG
    call monolis_std_debug_log_header("monolis_inner_product_main_R")
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
    call monolis_allreduce_R1(sum, monolis_mpi_sum, monoCOM%comm)
    t3 = monolis_get_time()

    if(present(tdotp)) tdotp = tdotp + t3 - t1
    if(present(tcomm)) tcomm = tcomm + t3 - t2
  end subroutine monolis_inner_product_main_R

  subroutine monolis_inner_product_main_R_no_comm(monoCOM, n, ndof, X, Y, sum)
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
  end subroutine monolis_inner_product_main_R_no_comm
end module mod_monolis_linalg