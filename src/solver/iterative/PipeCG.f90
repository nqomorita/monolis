module mod_monolis_solver_PipeCG
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

  subroutine monolis_solver_PipeCG(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: N, NP, NDOF, NNDOF
    integer(kind=kint) :: i, iter, iter_RR
    !integer(kind=kint) :: requests(1)
    !integer(kind=kint) :: statuses(monolis_status_size,1)
    real(kind=kdouble) :: alpha, alpha1, beta, gamma, gamma1, delta
    real(kind=kdouble) :: buf(3), CG(3), B2, R2
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
    call monolis_set_converge(monoPRM, monoCOM, monoMAT, R, B2, is_converge, monoPRM%tdotp, monoPRM%tcomm_dotp)
    if(is_converge) return
    call monolis_inner_product_R(monoCOM, N, NDOF, B, B, B2, monoPRM%tdotp, monoPRM%tcomm_dotp)
    call monolis_precond_apply(monoPRM, monoCOM, monoMAT, R, U)
    call monolis_matvec(monoCOM, monoMAT, U, V, monoPRM%tspmv, monoPRM%tcomm_spmv)

    gamma = 1.0d0
    alpha = 1.0d0

    do iter = 1, monoPRM%maxiter
      call monolis_inner_product_R_local(monoCOM, N, NDOF, R, U, CG(1))
      call monolis_inner_product_R_local(monoCOM, N, NDOF, V, U, CG(2))
      call monolis_inner_product_R_local(monoCOM, N, NDOF, R, R, CG(3))
      call monolis_allreduce_R(3, CG, monolis_sum, monoCOM%comm)

      call monolis_precond_apply(monoPRM, monoCOM, monoMAT, V, M)
      call monolis_matvec(monoCOM, monoMAT, M, L, monoPRM%tspmv, monoPRM%tcomm_spmv)

      gamma1 = 1.0d0/gamma
      alpha1 = 1.0d0/alpha

      !call monolis_wait(requests, statuses)
      !gamma = buf(1)
      !delta = buf(2)
      !R2    = buf(3)
      gamma = CG(1)
      delta = CG(2)
      R2    = CG(3)

      call monolis_check_converge_2(monoPRM, monoCOM, monoMAT, R2, B2, iter, is_converge)
      if(is_converge) exit

      if(1 < iter)then
        beta  = gamma*gamma1
        alpha = gamma/(delta-beta*gamma*alpha1)
      else
        alpha = gamma/delta
        beta  = 0.0d0
      endif

      do i = 1, NNDOF
         Z(i) = L(i) + beta*Z(i)
         Q(i) = M(i) + beta*Q(i)
         S(i) = V(i) + beta*S(i)
         P(i) = U(i) + beta*P(i)
      enddo

      if(mod(iter, iter_RR) == 0)then
        call monolis_vec_AXPY(N, NDOF, alpha, P, X, X)
        call monolis_residual(monoCOM, monoMAT, X, B, R, monoPRM%tspmv, monoPRM%tcomm_spmv)
        call monolis_precond_apply(monoPRM, monoCOM, monoMAT, R, U)
        call monolis_matvec(monoCOM, monoMAT, U, V, monoPRM%tspmv, monoPRM%tcomm_spmv)
      else
        do i = 1, NNDOF
          R(i) = R(i) - alpha*S(i)
          U(i) = U(i) - alpha*Q(i)
          V(i) = V(i) - alpha*Z(i)
          X(i) = X(i) + alpha*P(i)
        enddo
      endif
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
  end subroutine monolis_solver_PipeCG

end module mod_monolis_solver_PipeCG
