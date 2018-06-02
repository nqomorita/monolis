module mod_monolis_solver_CG

  implicit none

contains

  subroutine monolis_solver_CG(monoPRM, monoCOM, monoMAT)
    use mod_monolis_prm
    use mod_monolis_com
    use mod_monolis_mat
    use mod_monolis_precond
    use mod_monolis_matvec
    use mod_monolis_linalg
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: N, NP, NDOF, NNDOF
    integer(kind=kint) :: i, j, k, l, iter, iter_RR
    real(kind=kdouble) :: tol, resid, R2, B2
    real(kind=kdouble) :: t1, t2, tset, tsol, tcomm
    real(kind=kdouble) :: alpha, beta, rho, rho1, C1
    real(kind=kdouble), allocatable :: WW(:,:)
    real(kind=kdouble), pointer :: B(:), X(:)
    integer(kind=kint), parameter ::  R = 1
    integer(kind=kint), parameter ::  Z = 2
    integer(kind=kint), parameter ::  Q = 2
    integer(kind=kint), parameter ::  P = 3
    integer(kind=kint), parameter :: WK = 4

    t1 = monolis_wtime()

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NNDOF = N*NDOF
    X => monoMAT%X; X = 0.0d0
    B => monoMAT%B

    allocate(WW(NDOF*NP,4))
    WW = 0.0d0

    call monolis_precond_setup(monoPRM, monoCOM, monoMAT)
    call monolis_residual(monoCOM, monoMAT, X, B, tcomm)
    call monolis_inner_product_R(monoCOM, monoMAT, ndof, X, B, beta, tcomm)

    do iter=1, monoPRM%maxiter
      call monolis_precond_apply(monoPRM, monoCOM, monoMAT, X, B)
      call monolis_inner_product_R(monoCOM, monoMAT, ndof, X, B, beta, tcomm)

      if ( ITER.eq.1 ) then
        do i=1, NNDOF
          WW(i,P) = WW(i,Z)
        enddo
      else
        beta = rho/rho1
        do i = 1, NNDOF
          WW(i,P) = WW(i,Z) + beta*WW(i,P)
        enddo
      endif

      call monolis_matvec(monoCOM, monoMAT, X, B, tcomm)
      call monolis_inner_product_R(monoCOM, monoMAT, ndof, X, B, beta, tcomm)
      alpha = rho/C1

      do i=1, NNDOF
         X(i) = X(i) + alpha * WW(i,P)
      enddo

      if(mod(iter, iter_RR) == 0)then
        call monolis_residual(monoCOM, monoMAT, X, B, tcomm)
      else
        do i=1, NNDOF
          WW(i,R) = WW(i,R) - alpha * WW(i,Q)
        enddo
      endif

      call monolis_inner_product_R(monoCOM, monoMAT, ndof, X, B, beta, tcomm)
      resid = dsqrt(R2/B2)

      if(monoCOM%myrank == 0) write (*,'(i7, 1pe16.6)') iter, resid
      if(resid <= tol) exit

      rho1 = rho
    enddo

    !call monolis_update_R()

    t2 = monolis_wtime()
    tsol = t2 - t1
  end subroutine monolis_solver_CG

end module mod_monolis_solver_CG
