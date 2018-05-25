module mod_monolis_solver_GS

  implicit none

contains

  subroutine monolis_solver_GS(monoPRM, monoCOM, monoMAT)
    use mod_monolis_prm
    use mod_monolis_com
    use mod_monolis_mat
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: N, NP, NDOF, NDOF2, NNDOF
    integer(kind=kint) :: i, j, k, l, iter
    real(kind=kdouble) :: tol, resid, R2, B2
    real(kind=kdouble) :: t1, t2, tset, tsol
    real(kind=kdouble), pointer :: B(:), X(:)
    real(kind=kdouble), allocatable :: R(:), Dinv(:), LU(:,:), T(:)

    t1 = monolis_wtime()

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NDOF2 = NDOF*NDOF
    NNDOF = N*NDOF
    X => monoMAT%X; X = 0.0d0
    B => monoMAT%B

    allocate(R(NDOF*NP));     R    = 0.0d0
    allocate(T(NDOF));        T    = 0.0d0
    allocate(Dinv(NDOF2*NP)); Dinv = 0.0d0
    allocate(LU(NDOF,NDOF));  LU   = 0.0d0

    do i=1,NP
      do j=1,NDOF
        do k=1,NDOF
          LU(k,j) = monoMAT%D(NDOF2*(i-1)+NDOF*(j-1)+k)
        enddo
      enddo
      do k= 1, NDOF
        LU(k,k)= 1.0d0/LU(k,k)
        do l= k+1, NDOF
          LU(l,k)= LU(l,k) * LU(k,k)
          do j= k+1, NDOF
            T(j)= LU(l,j) - LU(l,k)*LU(k,j)
          enddo
          do j= k+1, NDOF
            LU(l,j)= T(j)
          enddo
        enddo
      enddo
      do j=1,NDOF
        do k=1,NDOF
          Dinv(NDOF2*(i-1)+NDOF*(j-1)+k) = LU(k,j)
        enddo
      enddo
    enddo

    !call monolis_inner_product_R()

    do iter=1, monoPRM%maxiter
      !call monolis_matvec_gs()
      !call monolis_residual()
      !call monolis_inner_product_R()
      resid = dsqrt(R2/B2)

      if(monoCOM%myrank == 0) write (*,"(i7, 1pe16.6)") iter, resid
      if(resid <= tol) exit
    enddo

    !call monolis_update_R()

    deallocate(R)
    deallocate(X)
    deallocate(Dinv)
    deallocate(LU)
    deallocate(T)

    t2 = monolis_wtime()
    tsol = t2 - t1
  end subroutine monolis_solver_GS

end module mod_monolis_solver_GS
