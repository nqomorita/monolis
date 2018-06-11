module mod_monolis_solver_BiCGSTAB
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

  subroutine monolis_solver_BiCGSTAB(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: N, NP, NDOF, NNDOF
    integer(kind=kint) :: i, iter, iter_RR
    real(kind=kdouble) :: t1, t2, tsol, tcomm
    real(kind=kdouble) :: alpha, beta, rho, rho1, c2, omega
    real(kind=kdouble) :: CG(2)
    real(kind=kdouble), allocatable :: W(:,:)
    real(kind=kdouble), pointer :: B(:), X(:)
    integer(kind=kint), parameter :: R  = 1
    integer(kind=kint), parameter :: RT = 2
    integer(kind=kint), parameter :: P  = 3
    integer(kind=kint), parameter :: PT = 4
    integer(kind=kint), parameter :: S  = 5
    integer(kind=kint), parameter :: ST = 1
    integer(kind=kint), parameter :: T  = 6
    integer(kind=kint), parameter :: V  = 7
    logical :: is_converge

    t1 = monolis_wtime()

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NNDOF = N*NDOF
    X => monoMAT%X; X = 1.0d0
    B => monoMAT%B

    allocate(W(NDOF*NP,7))
    W = 0.0d0

    iter_RR = 50

    call monolis_set_converge(monoPRM, monoCOM, monoMAT, B, tcomm)
    call monolis_residual(monoCOM, monoMAT, X, B, W(:,R), tcomm)

    do i = 1, NNDOF
      W(i,RT) = W(i,R)
    enddo

    do iter = 1, monoPRM%maxiter
      call monolis_inner_product_R(monoCOM, monoMAT, NDOF, W(:,R), W(:,RT), rho, tcomm)

      if(1 < iter)then
        beta = (rho/rho1) * (alpha/omega)
        do i = 1, NNDOF
          W(i,P) = W(i,R) + beta * (W(i,P) - omega * W(i,V))
        enddo
      else
        do i = 1, NNDOF
          W(i,P) = W(i,R)
        enddo
      endif

      call monolis_precond_apply(monoPRM, monoCOM, monoMAT, W(:,P), W(:,PT))
      call monolis_matvec(monoCOM, monoMAT, W(:,PT), W(:,V), tcomm)
      call monolis_inner_product_R(monoCOM, monoMAT, NDOF, W(:,RT), W(:,V), c2, tcomm)

      alpha = rho / c2

      do i = 1, NNDOF
        W(i,S) = W(i,R) - alpha * W(i,V)
      enddo

      call monolis_precond_apply(monoPRM, monoCOM, monoMAT, W(:,S), W(:,ST))
      call monolis_matvec(monoCOM, monoMAT, W(:,ST), W(:,T), tcomm)

      call monolis_inner_product_R_local(monoCOM, monoMAT, NDOF, W(:,T), W(:,S), CG(1))
      call monolis_inner_product_R_local(monoCOM, monoMAT, NDOF, W(:,T), W(:,T), CG(2))
      call monolis_allreduce_R(2, CG, monolis_sum, monoCOM%comm)

      omega = CG(1) / CG(2)

      do i = 1, NNDOF
        X(i) = X(i) + alpha*W(i,PT) + omega*W(i,ST)
      enddo

      if(mod(iter, iter_RR) == 0)then
        call monolis_residual(monoCOM, monoMAT, X, B, W(:,R), tcomm)
      else
        do i = 1, NNDOF
          W(i,R) = W(i,S) - omega * W(i,T)
        enddo
      endif

      call monolis_check_converge(monoPRM, monoCOM, monoMAT, W(:,R), iter, is_converge, tcomm)
      if(is_converge) exit

      rho1 = rho
    enddo

    call monolis_update_R(monoCOM, NDOF, X, tcomm)

    deallocate(W)

    t2 = monolis_wtime()
    tsol = t2 - t1
  end subroutine monolis_solver_BiCGSTAB

end module mod_monolis_solver_BiCGSTAB
