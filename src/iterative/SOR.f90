!> SOR 法モジュール
module mod_monolis_solver_SOR
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_precond
  use mod_monolis_matvec
  use mod_monolis_inner_product
  use mod_monolis_vec_util
  use mod_monolis_converge

  implicit none
  private
  public monolis_solver_SOR

contains

  subroutine monolis_solver_SOR(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in,out] 行列構造体
    type(monolis_mat), target, intent(inout) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), target, intent(inout) :: monoPREC
    integer(kint) :: N, NP, NDOF
    integer(kint) :: iter
    real(kdouble) :: R2, B2
    real(kdouble) :: tspmv, tdotp, tcomm_spmv, tcomm_dotp
    real(kdouble), pointer :: B(:), X(:)
    real(kdouble), allocatable :: R(:)
    logical :: is_converge

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF

    X => monoMAT%R%X
    B => monoMAT%R%B

    tspmv = monoPRM%Rarray(monolis_R_time_spmv)
    tcomm_spmv = monoPRM%Rarray(monolis_R_time_comm_spmv)
    tdotp = monoPRM%Rarray(monolis_R_time_dotp)
    tcomm_dotp = monoPRM%Rarray(monolis_R_time_comm_dotp)

    if(monoPRM%Iarray(monolis_prm_I_is_init_x) == monolis_I_true)then
      X = 0.0d0
    endif

    call monolis_alloc_R_1d(R, NDOF*NP)
    call monolis_palloc_R_1d(monoMAT%R%D, NDOF*N)

    call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
    call monolis_set_converge_R(monoCOM, monoMAT, R, B2, is_converge, tdotp, tcomm_dotp)
    if(is_converge) return

    call monolis_solver_SOR_setup(monoMAT)
    call monolis_inner_product_main_R(monoCOM, N*NDOF, B, B, B2, tdotp, tcomm_dotp)

    do iter = 1, monoPRM%Iarray(monolis_prm_I_max_iter)
      call monolis_solver_SOR_matvec(monoCOM, monoMAT, NDOF, X, B, tspmv, tcomm_spmv)
      call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
      call monolis_inner_product_main_R(monoCOM, N*NDOF, R, R, R2, tdotp, tcomm_dotp)
      call monolis_check_converge_R(monoPRM, monoCOM, monoMAT, R, B2, iter, is_converge, tdotp, tcomm_dotp)
      if(is_converge) exit
    enddo

    call monolis_mpi_update_R(monoCOM, NDOF, X, tcomm_spmv)

    call monolis_dealloc_R_1d(R)
    call monolis_pdealloc_R_1d(monoMAT%R%D)
  end subroutine monolis_solver_SOR

  subroutine monolis_solver_SOR_setup(monoMAT)
    implicit none
    type(monolis_mat) :: monoMAT
    integer(kint) :: i, ii, j, jS, jE, in, N, NDOF, NDOF2
    integer(kint), pointer :: index(:), item(:)
    real(kdouble), pointer :: A(:), D(:)

    N =  monoMAT%N
    NDOF =  monoMAT%NDOF
    NDOF2 = NDOF*NDOF
    A => monoMAT%R%A
    D => monoMAT%R%D
    index => monoMAT%CSR%index
    item => monoMAT%CSR%item

!$omp parallel default(none) &
!$omp & shared(A, D, index, item) &
!$omp & firstprivate(N, NDOF, NDOF2) &
!$omp & private(i, j, jS, jE, in)
!$omp do
    do i = 1, N
      jS = index(i) + 1
      jE = index(i + 1)
      do ii = jS, jE
        in = item(ii)
        if(i == in)then
          do j = 1, NDOF
            D(NDOF*(i-1)+j) = A(NDOF2*(ii-1) + (NDOF+1)*(j-1) + 1)
          enddo
        endif
      enddo
    enddo
!$omp end do
!$omp end parallel
  end subroutine monolis_solver_SOR_setup

  subroutine monolis_solver_SOR_matvec(monoCOM, monoMAT, NDOF, X, B, tspmv, tcomm)
    implicit none
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: N, NP, NDOF
    integer(kint) :: i
    real(kdouble) :: X(:), B(:)
    real(kdouble), allocatable :: Y(:)
    real(kdouble), pointer :: D(:)
    real(kdouble) :: omega
    real(kdouble) :: tspmv, tcomm

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    D => monoMAT%R%D

    omega = 1.0d0

    call monolis_alloc_R_1d(Y, NDOF*NP)

    call monolis_matvec_product_main_R(monoCOM, monoMAT, X, Y, tspmv, tcomm)

    do i = 1, N*NDOF
      X(i) = (1.0d0 - omega)*X(i) + omega*(B(i) - Y(i) + D(i)*X(i)) / D(i)
    enddo

    call monolis_dealloc_R_1d(Y)
  end subroutine monolis_solver_SOR_matvec
end module mod_monolis_solver_SOR
