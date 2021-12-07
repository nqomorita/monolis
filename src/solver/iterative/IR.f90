module mod_monolis_solver_IR
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_precond
  use mod_monolis_matvec
  use mod_monolis_linalg
  use mod_monolis_linalg_util
  use mod_monolis_converge

  implicit none

contains

  subroutine monolis_solver_IR(monoPRM, monoCOM, monoMAT)
    use mod_monolis_prm
    use mod_monolis_com
    use mod_monolis_mat
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: N, NP, NDOF, NNDOF
    integer(kind=kint) :: i, iter
    real(kind=kdouble) :: tol, resid, R2, B2
    real(kind=kdouble), pointer :: B(:), X(:)
    real(kind=kdouble), allocatable :: R(:), D(:)

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NNDOF = N*NDOF
    X => monoMAT%X
    B => monoMAT%B
    tol = monoPRM%tol

    if(monoPRM%is_init_x) X = 0.0d0

    allocate(R(NDOF*NP)); R = 0.0d0
    allocate(D(NDOF*NP)); D = 0.0d0

    call monolis_IR_setup(monoPRM, monoCOM, monoMAT)
    call monolis_residual(monoCOM, monoMAT, X, B, R, monoPRM%tspmv, monoPRM%tcomm_spmv)
    call monolis_inner_product_R(monoCOM, N, NDOF, R, R, B2, monoPRM%tdotp, monoPRM%tcomm_dotp)

    call monolis_vec_copy_R(N, NDOF, B, R)

    do iter = 1, monoPRM%maxiter
      call monolis_IR_apply(monoPRM, monoCOM, monoMAT, R, D)

      call monolis_vec_AXPY(N, NDOF, 1.0d0, X, D, X)

      call monolis_residual(monoCOM, monoMAT, X, B, R, monoPRM%tspmv, monoPRM%tcomm_spmv)
      call monolis_inner_product_R(monoCOM, N, NDOF, R, R, R2, monoPRM%tdotp, monoPRM%tcomm_dotp)

      resid = dsqrt(R2/B2)

      if(monoCOM%myrank == 0 .and. monoPRM%show_iterlog) write (*,"(i7, 1pe16.6)") iter, resid
      if(resid <= tol) exit
    enddo

    call monolis_IR_clear(monoPRM, monoCOM, monoMAT)
    call monolis_update_R(monoCOM, NDOF, X, monoPRM%tcomm_spmv)

    deallocate(R)
    deallocate(D)
  end subroutine monolis_solver_IR

  subroutine monolis_IR_setup(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    if(monoPRM%precond == monolis_prec_DIAG)then
      call monolis_precond_diag_setup(monoPRM, monoCOM, monoMAT)
    endif
  end subroutine monolis_IR_setup

  subroutine monolis_IR_apply(monoPRM, monoCOM, monoMAT, X, Y)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: i
    real(kind=kdouble) :: X(:), Y(:)

    if(monoPRM%precond == monolis_prec_DIAG)then
      call monolis_precond_diag_apply(monoPRM, monoCOM, monoMAT, X, Y)
    else
      call monolis_vec_copy_R(monoMAT%N, monoMAT%NDOF, X, Y)
    endif
  end subroutine monolis_IR_apply

  subroutine monolis_IR_clear(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    if(monoPRM%precond == monolis_prec_DIAG)then
      call monolis_precond_diag_clear(monoPRM, monoCOM, monoMAT)
    endif
  end subroutine monolis_IR_clear

end module mod_monolis_solver_IR
