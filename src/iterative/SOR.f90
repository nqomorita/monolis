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
    integer(kint) :: NNDOF, NPNDOF
    integer(kint) :: iter
    real(kdouble) :: R2, B2
    real(kdouble) :: tspmv, tdotp, tcomm_spmv, tcomm_dotp
    real(kdouble), pointer, contiguous :: B(:), X(:)
    real(kdouble), allocatable :: R(:)
    logical :: is_converge

    X => monoMAT%R%X
    B => monoMAT%R%B

    tspmv = monoPRM%Rarray(monolis_R_time_spmv)
    tcomm_spmv = monoPRM%Rarray(monolis_R_time_comm_spmv)
    tdotp = monoPRM%Rarray(monolis_R_time_dotp)
    tcomm_dotp = monoPRM%Rarray(monolis_R_time_comm_dotp)

    if(monoPRM%Iarray(monolis_prm_I_is_init_x) == monolis_I_true)then
      X = 0.0d0
    endif

    call monolis_get_vec_size(monoMAT%N, monoMAT%NP, monoMAT%NDOF, &
      monoMAT%n_dof_index, NNDOF, NPNDOF)

    call monolis_alloc_R_1d(R, NPNDOF)
    call monolis_palloc_R_1d(monoMAT%R%D, NNDOF)

    call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
    call monolis_set_converge_R(monoCOM, monoMAT, R, B2, is_converge, tdotp, tcomm_dotp)
    if(is_converge) return

    call monolis_solver_SOR_setup(monoMAT)
    call monolis_inner_product_main_R(monoCOM, NNDOF, B, B, B2, tdotp, tcomm_dotp)

    do iter = 1, monoPRM%Iarray(monolis_prm_I_max_iter)
      call monolis_solver_SOR_matvec(monoCOM, monoMAT, NNDOF, NPNDOF, X, B, tspmv, tcomm_spmv)
      call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
      call monolis_inner_product_main_R(monoCOM, NNDOF, R, R, R2, tdotp, tcomm_dotp)
      call monolis_check_converge_R(monoPRM, monoCOM, monoMAT, R, B2, iter, is_converge, tdotp, tcomm_dotp)
      if(is_converge) exit
    enddo

    call monolis_mpi_update_R_wrapper(monoCOM, monoMAT%NDOF, monoMAT%n_dof_index, X, tcomm_spmv)

    call monolis_dealloc_R_1d(R)
    call monolis_pdealloc_R_1d(monoMAT%R%D)
  end subroutine monolis_solver_SOR

  subroutine monolis_solver_SOR_setup(monoMAT)
    implicit none
    type(monolis_mat) :: monoMAT
    integer(kint) :: i, ii, j, jS, jE, in, N, n1, n2, nz
    integer(kint), pointer :: index(:), item(:)
    real(kdouble), pointer :: A(:), D(:)

    A => monoMAT%R%A
    D => monoMAT%R%D
    index => monoMAT%CSR%index
    item => monoMAT%CSR%item

    do i = 1, monoMAT%N
      jS = index(i) + 1
      jE = index(i + 1)
      n1 = monoMAT%n_dof_list(i)
      n2 = monoMAT%n_dof_index(i)
      do ii = jS, jE
        in = item(ii)
        nz = monoMAT%n_dof_index2(ii)
        if(i == in)then
          do j = 1, n1
            D(n2 + j) = A(nz + (n1+1)*(j-1) + 1)
          enddo
        endif
      enddo
    enddo
  end subroutine monolis_solver_SOR_setup

  subroutine monolis_solver_SOR_matvec(monoCOM, monoMAT, NNDOF, NPNDOF, X, B, tspmv, tcomm)
    implicit none
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: NNDOF, NPNDOF
    integer(kint) :: i
    real(kdouble) :: X(:), B(:)
    real(kdouble), allocatable :: Y(:)
    real(kdouble), pointer :: D(:)
    real(kdouble) :: omega
    real(kdouble) :: tspmv, tcomm

    D => monoMAT%R%D

    omega = 1.0d0

    call monolis_alloc_R_1d(Y, NPNDOF)

    call monolis_matvec_product_main_R(monoCOM, monoMAT, X, Y, tspmv, tcomm)

    do i = 1, NNDOF
      X(i) = (1.0d0 - omega)*X(i) + omega*(B(i) - Y(i) + D(i)*X(i)) / D(i)
    enddo

    call monolis_dealloc_R_1d(Y)
  end subroutine monolis_solver_SOR_matvec
end module mod_monolis_solver_SOR
