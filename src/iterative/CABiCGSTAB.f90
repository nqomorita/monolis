module mod_monolis_solver_CABiCGSTAB_noprec
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_precond
  use mod_monolis_matvec
  use mod_monolis_linalg
  use mod_monolis_linalg_util

  implicit none

contains

  subroutine monolis_solver_CABiCGSTAB_noprec(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: N, NP, NDOF, NNDOF
    integer(kind=kint) :: i, j, k, l, iter, iter_RR
    real(kind=kdouble) :: tol, resid, R2, B2, D2, tcomm
    real(kind=kdouble) :: t1, t2, tset, tsol
    real(kind=kdouble) :: alpha, beta, rho, rho1, C1, CG(5), omega
    real(kind=kdouble), allocatable :: W(:,:)
    real(kind=kdouble), pointer :: B(:), X(:)
    integer(kind=kint), parameter :: R0= 1
    integer(kind=kint), parameter :: R = 2
    integer(kind=kint), parameter :: V = 3
    integer(kind=kint), parameter :: S = 4
    integer(kind=kint), parameter :: P = 5
    integer(kind=kint), parameter :: Q = 6
    integer(kind=kint), parameter :: Y = 7
    integer(kind=kint), parameter :: Z = 8

    t1 = monolis_wtime()

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NNDOF = N*NDOF
    X => monoMAT%X; X = 1.0d0
    B => monoMAT%B

    allocate(W(NDOF*NP,9))
    W = 0.0d0

    iter_RR = 50
    tol = monoPRM%tol

    call monolis_residual(monoCOM, monoMAT, X, B, W(:,R0), tcomm)

    do i=1, NNDOF
      W(i,R) = W(i,R0)
    enddo

    call monolis_matvec(monoCOM, monoMAT, W(:,R0), W(:,V), tcomm)
    call monolis_inner_product_R(monoCOM, monoMAT, NDOF, W(:,R0), W(:,R0), CG(1), tcomm)
    call monolis_inner_product_R(monoCOM, monoMAT, NDOF, W(:,R0), W(:,V) , CG(2), tcomm)

    alpha = CG(1)/CG(2)
    beta  = 0.0d0
    omega = 0.0d0
    rho   = CG(1)

    call monolis_inner_product_R(monoCOM, monoMAT, NDOF, B, B, B2, tcomm)

    do iter=1, monoPRM%maxiter
      do i = 1, NNDOF
        W(i,P) = W(i,R) + beta * (W(i,P) - omega * W(i,S))
        W(i,S) = W(i,V) + beta * (W(i,S) - omega * W(i,Z))
      enddo

      call monolis_matvec(monoCOM, monoMAT, W(:,S), W(:,Z), tcomm)

      do i = 1, NNDOF
        W(i,Q) = W(i,R) - alpha * W(i,S)
        W(i,Y) = W(i,V) - alpha * W(i,Z)
      enddo

      call monolis_inner_product_R(monoCOM, monoMAT, NDOF, W(:,Q), W(:,Y), CG(1), tcomm)
      call monolis_inner_product_R(monoCOM, monoMAT, NDOF, W(:,Y), W(:,Y), CG(2), tcomm)

      omega = CG(1)/CG(2)

      do i = 1, NNDOF
        X(i)   = X(i)   + alpha * W(i,P) + omega * W(i,Q)
        W(i,R) = W(i,Q) - omega * W(i,Y)
      enddo

      call monolis_matvec(monoCOM, monoMAT, W(:,R), W(:,V), tcomm)

      call monolis_inner_product_R(monoCOM, monoMAT, NDOF, W(:,R0), W(:,R), CG(1), tcomm)
      call monolis_inner_product_R(monoCOM, monoMAT, NDOF, W(:,R0), W(:,V), CG(2), tcomm)
      call monolis_inner_product_R(monoCOM, monoMAT, NDOF, W(:,R0), W(:,S), CG(3), tcomm)
      call monolis_inner_product_R(monoCOM, monoMAT, NDOF, W(:,R0), W(:,Z), CG(4), tcomm)
      call monolis_inner_product_R(monoCOM, monoMAT, NDOF, W(:,R ), W(:,R), CG(5), tcomm)

      beta  = (alpha/omega) * CG(1) / rho
      alpha = CG(1) / (CG(2) + beta * CG(3) - beta * omega * CG(4))
      R2    = CG(5)
      resid = dsqrt(R2/B2)

      if(monoCOM%myrank == 0) write (*,"(i7, 1pe16.6)") iter, resid
      if(resid <= tol) exit

      rho = CG(1)
    enddo

    call monolis_update_R(monoCOM, NDOF, X, tcomm)

    deallocate(W)

    t2 = monolis_wtime()
    tsol = t2 - t1
  end subroutine monolis_solver_CABiCGSTAB_noprec

end module mod_monolis_solver_CABiCGSTAB_noprec
