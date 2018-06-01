module mod_monolis_solver_BiCGSTAB

  implicit none

contains

  subroutine monolis_solver_BiCGSTAB(monoPRM, monoCOM, monoMAT)
    use mod_monolis_prm
    use mod_monolis_com
    use mod_monolis_mat
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: N, NP, NDOF, NNDOF
    integer(kind=kint) :: i, j, k, l, iter, iter_RR
    real(kind=kdouble) :: tol, resid, R2, B2
    real(kind=kdouble) :: t1, t2, tset, tsol
    real(kind=kdouble) :: alpha, beta, rho, rho1, C1, C2, omega
    real(kind=kdouble) :: CG(2)
    real(kind=kdouble), allocatable :: WW(:,:)
    real(kind=kdouble), pointer :: B(:), X(:)
    integer(kind=kint), parameter :: R = 1
    integer(kind=kint), parameter :: RT= 2
    integer(kind=kint), parameter :: P = 3
    integer(kind=kint), parameter :: PT= 4
    integer(kind=kint), parameter :: S = 5
    integer(kind=kint), parameter :: ST= 1
    integer(kind=kint), parameter :: T = 6
    integer(kind=kint), parameter :: V = 7
    integer(kind=kint), parameter :: WK= 8

    t1 = monolis_wtime()

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NNDOF = N*NDOF
    X => monoMAT%X; X = 0.0d0
    B => monoMAT%B

    allocate(WW(NDOF*NP,4))
    WW = 0.0d0

    !call monolis_precond_setup()
    !call monolis_residual()

    do i=1, NNDOF
      WW(i,RT) = WW(i,R)
    enddo

    !call monolis_inner_product_R()

    do iter=1, monoPRM%maxiter
      !call monolis_inner_product_R()

      if ( iter.gt.1 ) then
        beta = (RHO/RHO1) * (alpha/omega)
        do i = 1, NNDOF
          WW(i,P) = WW(i,R) + beta * (WW(i,P) - omega * WW(i,V))
        enddo
       else
        do i = 1, NNDOF
          WW(i,P) = WW(i,R)
        enddo
      endif

      !call monolis_precond_apply()
      !call monolis_matvec()
      !call hecmw_inner_product_R()

      alpha = RHO / C2

      do i = 1, NNDOF
        WW(i,S) = WW(i,R) - alpha * WW(i,V)
      enddo

      !call monolis_precond_apply()
      !call monolis_matvec()

      !call hecmw_inner_product_R_nocomm()
      !call hecmw_inner_product_R_nocomm()

      omega = CG(1) / CG(2)

      do i = 1, NNDOF
        X (i) = X(i) + alpha * WW(i,PT) + omega * WW(i,ST)
      enddo

      if(mod(iter, iter_RR) == 0)then
        !call monolis_residual()
      else
        do i=1, NNDOF
          !WW(i,R) = WW(i,R) - alpha * WW(i,Q)
        enddo
      endif

      !call monolis_inner_product()
      resid = dsqrt(R2/B2)

      if(monoCOM%myrank == 0) write (*,'(i7, 1pe16.6)') iter, resid
      if(resid <= tol) exit

      rho1 = rho
    enddo

    !call monolis_update_R()

    t2 = monolis_wtime()
    tsol = t2 - t1

  end subroutine monolis_solver_BiCGSTAB

end module mod_monolis_solver_BiCGSTAB
