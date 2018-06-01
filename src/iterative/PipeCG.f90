module mod_monolis_solver_PipeCG

  implicit none

contains

  subroutine monolis_solver_PipeCG(monoPRM, monoCOM, monoMAT)
    use mod_monolis_prm
    use mod_monolis_com
    use mod_monolis_mat
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: N, NP, NDOF, NNDOF
    integer(kind=kint) :: i, j, k, iter, iter_RR
    integer(kind=kint ) :: requests(1)
    integer(kind=kint ) :: statuses(monolis_status_size,1)
    real(kind=kdouble) :: tol, resid, R2, B2
    real(kind=kdouble) :: t1, t2, tset, tsol, D2
    real(kind=kdouble) :: alpha, alpha1, beta, rho, rho1, gamma, gamma1, delta, C1
    real(kind=kdouble) :: buf(3), CG(3)
    real(kind=kdouble), allocatable :: WW(:,:)
    real(kind=kdouble), pointer :: B(:), X(:)
    !integer(kind=kint), parameter ::  R = 1
    !integer(kind=kint), parameter ::  Z = 2
    !integer(kind=kint), parameter ::  Q = 2
    !integer(kind=kint), parameter ::  P = 3
    !integer(kind=kint), parameter :: WK = 4

    real(kind=kdouble), allocatable :: R(:)
    real(kind=kdouble), allocatable :: U(:)
    real(kind=kdouble), allocatable :: W(:)
    real(kind=kdouble), allocatable :: Q(:)
    real(kind=kdouble), allocatable :: P(:)
    real(kind=kdouble), allocatable :: Z(:)
    real(kind=kdouble), allocatable :: L(:)
    real(kind=kdouble), allocatable :: M(:)
    real(kind=kdouble), allocatable :: S(:)
    real(kind=kdouble), allocatable :: WK(:)
    real(kind=kdouble), allocatable  :: D(:), E(:)


    !integer(kind=kint), parameter ::  R= 1
    !integer(kind=kint), parameter ::  U= 2
    !integer(kind=kint), parameter ::  W= 3
    !integer(kind=kint), parameter ::  Q= 4
    !integer(kind=kint), parameter ::  P= 5
    !integer(kind=kint), parameter ::  Z= 6
    !integer(kind=kint), parameter ::  L= 7
    !integer(kind=kint), parameter ::  M= 8
    !integer(kind=kint), parameter ::  S= 9
    !integer(kind=kint), parameter :: WK=10

    t1 = monolis_wtime()

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NNDOF = N*NDOF
    X => monoMAT%X; X = 0.0d0
    B => monoMAT%B

    allocate(R(NDOF*NP))
    allocate(U(NDOF*NP))
    allocate(W(NDOF*NP))
    allocate(Q(NDOF*NP))
    allocate(P(NDOF*NP))
    allocate(Z(NDOF*NP))
    allocate(L(NDOF*NP))
    allocate(M(NDOF*NP))
    allocate(S(NDOF*NP))
    allocate(WK(NDOF*NP))
    R = 0.d0; U = 0.d0; W = 0.d0; Q = 0.d0; P = 0.d0
    Z = 0.d0; L = 0.d0; M = 0.d0; S = 0.d0; WK= 0.d0

    !call monolis_precond_setup()
    !call monolis_residual()
    !call monolis_inner_product_R()
    !call monolis_precond_apply()
    !call monolis_matvec()

    gamma = 1.0d0
    alpha = 1.0d0

    do iter=1, monoPRM%maxiter
      !call monolis_inner_product_R_nocomm()
      !call monolis_inner_product_R_nocomm()
      !call monolis_inner_product_R_nocomm()
      !call monolis_Iallreduce_R()

      !call monolis_precond_apply()
      !call monolis_matvec()

      gamma1 = 1.0d0/gamma
      alpha1 = 1.0d0/alpha

      !call monolis_wait(requests, statuses)
      gamma = buf(1)
      delta = buf(2)
      D2    = buf(3)

      if(iter > 1)then
        beta  = gamma*gamma1
        alpha = gamma/(delta-beta*gamma*alpha1)
      else
        alpha = gamma/delta
        beta  = 0.0d0
      endif

      do i = 1, NNDOF
         Z(i) = L(i) + beta * Z(i)
         Q(i) = M(i) + beta * Q(i)
         S(i) = W(i) + beta * S(i)
         P(i) = U(i) + beta * P(i)
      enddo

      if(mod(iter, iter_RR) == 0)then
        do i = 1, NNDOF
         X(i) = X(i) + alpha * P(i)
        enddo
        !call monolis_residual()
        !call monolis_precond_apply()
      else
        do i = 1, NNDOF
          R(i) = R(i) - alpha * S(i)
          U(i) = U(i) - alpha * Q(i)
          W(i) = W(i) - alpha * Z(i)
          X(i) = X(i) + alpha * P(i)
        enddo
      endif

      resid = dsqrt(R2/B2)

      if(monoCOM%myrank == 0) write (*,'(i7, 1pe16.6)') iter, resid
    enddo

    !call monolis_update_R()

    t2 = monolis_wtime()
    tsol = t2 - t1
  end subroutine monolis_solver_PipeCG

end module mod_monolis_solver_PipeCG
