module mod_monolis_solver_CG
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_precond
  use mod_monolis_matvec
  use mod_monolis_linalg
  use mod_monolis_linalg_util
  use mod_monolis_converge
  use mod_monolis_util

  implicit none

contains

  subroutine monolis_solver_CG(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: N, NP, NDOF, NNDOF
    integer(kind=kint) :: i, iter, iter_RR
    real(kind=kdouble) :: alpha, beta, rho, rho1, omega, B2
    real(kind=kdouble), allocatable :: R(:), Z(:), Q(:), P(:)
    real(kind=kdouble), pointer :: B(:), X(:)
    logical :: is_converge

    if(monoPRM%is_debug) call monolis_debug_header("monolis_solver_CG")

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NNDOF = N*NDOF
    X => monoMAT%X
    B => monoMAT%B
    iter_RR = 50

    if(monoPRM%is_init_x) X = 0.0d0

    allocate(R(NDOF*NP), source = 0.0d0)
    allocate(Z(NDOF*NP), source = 0.0d0)
    allocate(Q(NDOF*NP), source = 0.0d0)
    allocate(P(NDOF*NP), source = 0.0d0)

    call monolis_set_converge(monoPRM, monoCOM, monoMAT, B, B2, is_converge)
    if(is_converge) return
    call monolis_residual(monoCOM, monoMAT, X, B, R)

    do iter = 1, monoPRM%maxiter
      call monolis_precond_apply(monoPRM, monoCOM, monoMAT, R, Z)
      call monolis_inner_product_R(monoCOM, N, NDOF, R, Z, rho, monoPRM%tdotp)

      if(1 < iter)then
        beta = rho/rho1
        call monolis_vec_AXPY(N, NDOF, beta, P, Z, P)
      else
        call monolis_vec_copy_R(N, NDOF, Z, P)
      endif

      call monolis_matvec(monoCOM, monoMAT, P, Q, monoPRM%tspmv)
      call monolis_inner_product_R(monoCOM, N, NDOF, P, Q, omega, monoPRM%tdotp)
      alpha = rho/omega

      call monolis_vec_AXPY(N, NDOF, alpha, P, X, X)

      if(mod(iter, iter_RR) == 0)then
        call monolis_residual(monoCOM, monoMAT, X, B, R)
      else
        call monolis_vec_AXPY(N, NDOF, -alpha, Q, R, R)
      endif

      call monolis_check_converge(monoPRM, monoCOM, monoMAT, R, B2, iter, is_converge)
      if(is_converge) exit

      rho1 = rho
    enddo

    call monolis_update_R(monoCOM, NDOF, X)

    deallocate(R)
    deallocate(Z)
    deallocate(Q)
    deallocate(P)
  end subroutine monolis_solver_CG

end module mod_monolis_solver_CG
