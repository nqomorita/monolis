module mod_monolis_solver_CABiCGSTAB

  implicit none

contains

  subroutine monolis_solver_CABiCGSTAB(monoPRM, monoCOM, monoMAT)
    use mod_monolis_prm
    use mod_monolis_com
    use mod_monolis_mat
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: N, NP, NDOF, NNDOF
    integer(kind=kint) :: i, j, k, l, iter, iter_RR
    real(kind=kdouble) :: tol, resid, R2, B2, D2
    real(kind=kdouble) :: t1, t2, tset, tsol
    real(kind=kdouble) :: alpha, beta, rho, rho1, C1, CG(5), omega
    real(kind=kdouble), allocatable :: WW(:,:)
    real(kind=kdouble), pointer :: B(:), X(:)
    integer(kind=kint), parameter :: R0= 1
    integer(kind=kint), parameter :: R = 2
    integer(kind=kint), parameter :: W = 3
    integer(kind=kint), parameter :: S = 4
    integer(kind=kint), parameter :: P = 5
    integer(kind=kint), parameter :: Q = 6
    integer(kind=kint), parameter :: Y = 7
    integer(kind=kint), parameter :: Z = 8
    integer(kind=kint), parameter :: WK= 9

    t1 = monolis_wtime()

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NNDOF = N*NDOF
    X => monoMAT%X; X = 0.0d0
    B => monoMAT%B

    allocate(WW(NDOF*NP,9))
    WW = 0.0d0

    !call monolis_residual()

    do i=1, NNDOF
      WW(i,R) = WW(i,R0)
    enddo

    !call monolis_matvec()
    !call monolis_inner_product_R()
    !call monolis_inner_product_R()

    ALPHA = CG(1)/CG(2)
    BETA  = 0.0d0
    OMEGA = 0.0d0
    RHO   = CG(1)

    !call monolis_inner_product_R()

    do iter=1, monoPRM%maxiter
      do i = 1, NNDOF
        WW(i,P) = WW(i,R) + BETA * (WW(i,P) - OMEGA * WW(i,S))
        WW(i,S) = WW(i,W) + BETA * (WW(i,S) - OMEGA * WW(i,Z))
      enddo

      !call monolis_matvec()

      do i = 1, NNDOF
        WW(i,Q) = WW(i,R) - ALPHA * WW(i,S)
        WW(i,Y) = WW(i,W) - ALPHA * WW(i,Z)
      enddo

      !call monolis_inner_product_R_nocomm()
      !call monolis_inner_product_R_nocomm()

      OMEGA = CG(1)/CG(2)

      do i = 1, NNDOF
        X (i)   = X (i)   + ALPHA * WW(i,P) + OMEGA * WW(i,Q)
        WW(i,R) = WW(i,Q) - OMEGA * WW(i,Y)
      enddo

      !call monolis_matvec()

      !call monolis_inner_product_R_nocomm()
      !call monolis_inner_product_R_nocomm()
      !call monolis_inner_product_R_nocomm()
      !call monolis_inner_product_R_nocomm()
      !call monolis_inner_product_R_nocomm()

      BETA  = (ALPHA/OMEGA) * CG(1) / RHO
      ALPHA = CG(1) / (CG(2) + BETA * CG(3) - BETA * OMEGA * CG(4))
      D2 = CG(5)
      resid = dsqrt(R2/B2)

      if(monoCOM%myrank == 0) write (*,'(i7, 1pe16.6)') iter, resid
      if(resid <= tol) exit

      rho1 = CG(1)
    enddo

    !call monolis_update_R()

    t2 = monolis_wtime()
    tsol = t2 - t1
  end subroutine monolis_solver_CABiCGSTAB

end module mod_monolis_solver_CABiCGSTAB
