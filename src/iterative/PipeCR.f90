module mod_monolis_solver_PipeCR
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

  subroutine monolis_solver_PipeCR(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: N, NP, NDOF, NNDOF
    integer(kind=kint) :: i, iter, iter_RR
    !integer(kind=kint) :: requests(1)
    !integer(kind=kint) :: statuses(monolis_status_size,1)
    real(kind=kdouble) :: tol, resid, R2, B2, U2
    real(kind=kdouble) :: t1, t2, tsol, tcomm
    real(kind=kdouble) :: alpha, alpha1, beta, gamma, gamma1, delta, phi, utol
    real(kind=kdouble) :: buf(3), CG(3)
    real(kind=kdouble), allocatable :: W(:,:)
    real(kind=kdouble), pointer :: B(:), X(:)
    integer(kind=kint), parameter :: R = 1
    integer(kind=kint), parameter :: U = 2
    integer(kind=kint), parameter :: V = 3
    integer(kind=kint), parameter :: Q = 4
    integer(kind=kint), parameter :: P = 5
    integer(kind=kint), parameter :: Z = 6
    integer(kind=kint), parameter :: L = 7
    integer(kind=kint), parameter :: M = 8
    integer(kind=kint), parameter :: S = 9
    logical :: is_converge

    t1 = monolis_wtime()

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NNDOF = N*NDOF
    X => monoMAT%X; X = 1.0d0
    B => monoMAT%B

    allocate(W(NDOF*NP, 9))
    W = 0.0d0

    iter_RR = 50
    tol = monoPRM%tol

    call monolis_residual(monoCOM, monoMAT, X, B, W(:,R), tcomm)
    call monolis_inner_product_R(monoCOM, monoMAT, NDOF, B, B, B2, tcomm)
    call monolis_inner_product_R(monoCOM, monoMAT, NDOF, W(:,R), W(:,R), R2, tcomm)
    call monolis_precond_apply(monoPRM, monoCOM, monoMAT, W(:,R), W(:,U))
    call monolis_inner_product_R(monoCOM, monoMAT, NDOF, W(:,U), W(:,U), U2, tcomm)
    call monolis_matvec(monoCOM, monoMAT, W(:,U), W(:,V), tcomm)

    phi  = dsqrt(R2/U2)
    utol = tol/phi

    do iter=1, monoPRM%maxiter
      call monolis_precond_apply(monoPRM, monoCOM, monoMAT, W(:,V), W(:,M))

      call monolis_inner_product_R_local(monoCOM, monoMAT, NDOF, W(:,V), W(:,U), CG(1))
      call monolis_inner_product_R_local(monoCOM, monoMAT, NDOF, W(:,V), W(:,M), CG(2))
      call monolis_inner_product_R_local(monoCOM, monoMAT, NDOF, W(:,U), W(:,U), CG(3))
      call monolis_allreduce_R(3, CG, monolis_sum, monoCOM%comm)

      call monolis_matvec(monoCOM, monoMAT, W(:,M), W(:,L), tcomm)

      !call monolis_wait(requests, statuses)
      !gamma = buf(1)
      !delta = buf(2)
      !U2    = buf(3)
      gamma = CG(1)
      delta = CG(2)
      U2    = CG(3)

      if(1 < iter)then
        beta  = gamma*gamma1
        alpha = gamma/(delta-beta*gamma*alpha1)
      else
        beta  = 0.0d0
        alpha = gamma/delta
      endif

      do i = 1, NNDOF
         W(i,Z) = W(i,L) + beta * W(i,Z)
         W(i,Q) = W(i,M) + beta * W(i,Q)
         W(i,P) = W(i,U) + beta * W(i,P)
         X(i)   = X(i)   + alpha * W(i,P)
      enddo

      if(mod(iter, iter_RR) == 0)then
        call monolis_residual(monoCOM, monoMAT, X, B, W(:,R), tcomm)
        call monolis_precond_apply(monoPRM, monoCOM, monoMAT, W(:,R), W(:,U))
        call monolis_matvec(monoCOM, monoMAT, W(:,U), W(:,V), tcomm)
      else
        do i = 1, NNDOF
          W(i,U) = W(i,U) - alpha * W(i,Q)
          W(i,V) = W(i,V) - alpha * W(i,Z)
        enddo
      endif

      resid = dsqrt(U2/B2)
      if(monoCOM%myrank == 0) write (*,"(i7, 1pe16.6)") iter, resid*phi
      if(resid <= utol) exit

      gamma1 = 1.0d0/gamma
      alpha1 = 1.0d0/alpha
    enddo

    call monolis_update_R(monoCOM, NDOF, X, tcomm)

    deallocate(W)

    t2 = monolis_wtime()
    tsol = t2 - t1
  end subroutine monolis_solver_PipeCR

end module mod_monolis_solver_PipeCR
