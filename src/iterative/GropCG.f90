module mod_monolis_solver_GropCG
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

  subroutine monolis_solver_GropCG(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: N, NP, NDOF, NNDOF
    integer(kind=kint) :: i, iter, iter_RR
    !integer(kind=kint) :: reqs1(1), reqs2(1)
    !integer(kind=kint) :: statuses(monolis_status_size,1)
    real(kind=kdouble) :: R2
    real(kind=kdouble) :: t1, t2, tsol, tcomm
    real(kind=kdouble) :: alpha, beta, delta, gamma, gamma1
    !real(kind=kdouble) :: buf1(1), buf2(2)
    real(kind=kdouble), allocatable :: R(:), U(:), V(:), Q(:), P(:), S(:)
    real(kind=kdouble), pointer :: B(:), X(:)
    logical :: is_converge

    t1 = monolis_wtime()

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NNDOF = N*NDOF
    X => monoMAT%X
    B => monoMAT%B
    iter_RR = 50

    if(monoPRM%is_init_x) X = 0.0d0

    allocate(R(NDOF*NP)); R = 0.0d0
    allocate(U(NDOF*NP)); U = 0.0d0
    allocate(V(NDOF*NP)); V = 0.0d0
    allocate(Q(NDOF*NP)); Q = 0.0d0
    allocate(P(NDOF*NP)); P = 0.0d0
    allocate(S(NDOF*NP)); S = 0.0d0

    call monolis_set_converge(monoPRM, monoCOM, monoMAT, B, tcomm)
    call monolis_residual(monoCOM, monoMAT, X, B, R, tcomm)
    call monolis_precond_apply(monoPRM, monoCOM, monoMAT, R, U)

    do i = 1,NDOF*NP
      P(i) = U(i)
    enddo

    call monolis_matvec(monoCOM, monoMAT, P, S, tcomm)
    call monolis_inner_product_R(monoCOM, N, NDOF, R, U, gamma, tcomm)

    do iter = 1, monoPRM%maxiter
      call monolis_inner_product_R(monoCOM, N, NDOF, P, S, delta, tcomm)
      call monolis_precond_apply(monoPRM, monoCOM, monoMAT, S, Q)

      alpha = gamma/delta

      do i = 1, NNDOF
        X(i) = X(i) + alpha*P(i)
        R(i) = R(i) - alpha*S(i)
        U(i) = U(i) - alpha*Q(i)
      enddo

      call monolis_inner_product_R(monoCOM, N, NDOF, R, U, gamma1, tcomm)
      call monolis_inner_product_R(monoCOM, N, NDOF, R, R, R2, tcomm)

      call monolis_matvec(monoCOM, monoMAT, U, V, tcomm)

      beta  = gamma1/gamma
      gamma = gamma1

      do i = 1, NNDOF
        P(i) = U(i) + beta * P(i)
        S(i) = V(i) + beta * S(i)
      enddo

      call monolis_check_converge_2(monoPRM, monoCOM, monoMAT, R2, iter, is_converge, tcomm)
      if(is_converge) exit
    enddo

    call monolis_update_R(monoCOM, NDOF, X, tcomm)

    deallocate(R)
    deallocate(U)
    deallocate(V)
    deallocate(Q)
    deallocate(P)
    deallocate(S)

    t2 = monolis_wtime()
    tsol = t2 - t1
  end subroutine monolis_solver_GropCG

end module mod_monolis_solver_GropCG
