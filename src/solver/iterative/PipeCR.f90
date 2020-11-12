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
    real(kind=kdouble) :: alpha, alpha1, beta, gamma, gamma1, delta, phi, utol
    real(kind=kdouble) :: buf(3), CG(3)
    real(kind=kdouble), allocatable :: R(:), U(:), V(:), Q(:), P(:), Z(:), L(:), M(:), S(:)
    real(kind=kdouble), pointer :: B(:), X(:)
    logical :: is_converge

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NNDOF = N*NDOF
    X => monoMAT%X
    B => monoMAT%B
    iter_RR = 50
    tol = monoPRM%tol

    if(monoPRM%is_init_x) X = 0.0d0

    allocate(R(NDOF*NP), source = 0.0d0)
    allocate(U(NDOF*NP), source = 0.0d0)
    allocate(V(NDOF*NP), source = 0.0d0)
    allocate(Q(NDOF*NP), source = 0.0d0)
    allocate(P(NDOF*NP), source = 0.0d0)
    allocate(Z(NDOF*NP), source = 0.0d0)
    allocate(L(NDOF*NP), source = 0.0d0)
    allocate(M(NDOF*NP), source = 0.0d0)
    allocate(S(NDOF*NP), source = 0.0d0)

    call monolis_residual(monoCOM, monoMAT, X, B, R, monoPRM%tspmv, monoPRM%tcomm_spmv)
    call monolis_inner_product_R(monoCOM, N, NDOF, B, B, B2, monoPRM%tdotp, monoPRM%tcomm_dotp)
    call monolis_inner_product_R(monoCOM, N, NDOF, R, R, R2, monoPRM%tdotp, monoPRM%tcomm_dotp)
    call monolis_precond_apply(monoPRM, monoCOM, monoMAT, R, U)
    call monolis_inner_product_R(monoCOM, N, NDOF, U, U, U2, monoPRM%tdotp, monoPRM%tcomm_dotp)
    call monolis_matvec(monoCOM, monoMAT, U, V, monoPRM%tspmv, monoPRM%tcomm_spmv)

    phi  = dsqrt(R2/U2)
    utol = tol/phi

    do iter = 1, monoPRM%maxiter
      call monolis_precond_apply(monoPRM, monoCOM, monoMAT, V, M)

      call monolis_inner_product_R_local(monoCOM, N, NDOF, V, U, CG(1))
      call monolis_inner_product_R_local(monoCOM, N, NDOF, V, M, CG(2))
      call monolis_inner_product_R_local(monoCOM, N, NDOF, U, U, CG(3))
      call monolis_allreduce_R(3, CG, monolis_sum, monoCOM%comm)

      call monolis_matvec(monoCOM, monoMAT, M, L, monoPRM%tspmv, monoPRM%tcomm_spmv)

      !call monolis_wait(requests, statuses)
      !gamma = buf(1)
      !delta = buf(2)
      !U2    = buf(3)
      gamma = CG(1)
      delta = CG(2)
      U2    = CG(3)

      resid = dsqrt(U2/B2)
      if(monoCOM%myrank == 0 .and. monoPRM%show_iterlog) write (*,"(i7, 1pe16.6)") iter, resid*phi
      if(resid <= utol) exit

      if(1 < iter)then
        beta  = gamma*gamma1
        alpha = gamma/(delta-beta*gamma*alpha1)
      else
        beta  = 0.0d0
        alpha = gamma/delta
      endif

      do i = 1, NNDOF
         Z(i) = L(i) + beta*Z(i)
         Q(i) = M(i) + beta*Q(i)
         P(i) = U(i) + beta*P(i)
         X(i) = X(i) + alpha*P(i)
      enddo

      if(mod(iter, iter_RR) == 0)then
        call monolis_residual(monoCOM, monoMAT, X, B, R, monoPRM%tspmv, monoPRM%tcomm_spmv)
        call monolis_precond_apply(monoPRM, monoCOM, monoMAT, R, U)
        call monolis_matvec(monoCOM, monoMAT, U, V, monoPRM%tspmv, monoPRM%tcomm_spmv)
      else
        do i = 1, NNDOF
          U(i) = U(i) - alpha*Q(i)
          V(i) = V(i) - alpha*Z(i)
        enddo
      endif

      gamma1 = 1.0d0/gamma
      alpha1 = 1.0d0/alpha
    enddo

    call monolis_update_R(monoCOM, NDOF, X, monoPRM%tcomm_spmv)

    deallocate(R)
    deallocate(U)
    deallocate(V)
    deallocate(Q)
    deallocate(P)
    deallocate(Z)
    deallocate(L)
    deallocate(M)
    deallocate(S)
  end subroutine monolis_solver_PipeCR

end module mod_monolis_solver_PipeCR
