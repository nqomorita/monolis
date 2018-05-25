module mod_monolis_solver_GropCG

  implicit none

contains

  subroutine monolis_solver_GropCG(monoPRM, monoCOM, monoMAT)
    use mod_monolis_prm
    use mod_monolis_com
    use mod_monolis_mat
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: N, NP, NDOF, NNDOF
    integer(kind=kint) :: i, j, k, l, iter, iter_RR
    integer(kind=kint) :: reqs1(1)
    integer(kind=kint) :: reqs2(1)
    integer(kind=kint) :: statuses(monolis_status_size,1)
    real(kind=kdouble) :: tol, resid, R2, B2
    real(kind=kdouble) :: t1, t2, tset, tsol
    real(kind=kdouble) :: alpha, beta, delta, gamma, gamma1, C1
    real(kind=kdouble) :: buf1(1), buf2(2)
    real(kind=kdouble), pointer :: B(:), X(:)
    real(kind=kdouble), allocatable :: R(:)
    real(kind=kdouble), allocatable :: U(:)
    real(kind=kdouble), allocatable :: W(:)
    real(kind=kdouble), allocatable :: Q(:)
    real(kind=kdouble), allocatable :: P(:)
    real(kind=kdouble), allocatable :: S(:)
    real(kind=kdouble), allocatable :: WK(:)

    t1 = monolis_wtime()

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NNDOF = N*NDOF
    X => monoMAT%X; X = 0.0d0
    B => monoMAT%B

    allocate(R (NDOF*NP)); R = 0.0d0
    allocate(U (NDOF*NP)); U = 0.0d0
    allocate(W (NDOF*NP)); W = 0.0d0
    allocate(Q (NDOF*NP)); Q = 0.0d0
    allocate(P (NDOF*NP)); P = 0.0d0
    allocate(S (NDOF*NP)); S = 0.0d0
    allocate(WK(NDOF*NP)); WK= 0.0d0

    !call monolis_precond_setup()
    !call monolis_residual()
    !call monolis_inner_product_R()
    !call monolis_precond_apply()

    do i=1,NDOF*NP
      P(i) = U(i)
    enddo

    !call monolis_matvec()
    !call monolis_inner_product_R()

    do iter = 1, monoPRM%maxiter
      !call monolis_inner_product_R()
      !call monolis_precond_apply()

      alpha = gamma/delta

      do i = 1, NNDOF
        X(i) = X(i) + alpha * P(i)
        R(i) = R(i) - alpha * S(i)
        U(i) = U(i) - alpha * Q(i)
      enddo

      !call monolis_inner_product_R()
      !call monolis_inner_product_R()

      !call monolis_matvec()

      beta  = gamma1/gamma
      gamma = gamma1

      do i = 1, NNDOF
        P(i) = U(i) + beta * P(i)
        S(i) = W(i) + beta * S(i)
      enddo

      resid = dsqrt(R2/B2)
      if(monoCOM%myrank == 0) write (*,'(i7, 1pe16.6)') iter, resid
      if(resid <= tol) exit
    enddo

    !call monolis_update_R()

    deallocate(R,U,W,Q,P,S,WK)

    t2 = monolis_wtime()
    tsol = t2 - t1
  end subroutine monolis_solver_GropCG

end module mod_monolis_solver_GropCG
