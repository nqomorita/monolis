module mod_monolis_solver_GS

  implicit none

contains

  subroutine monolis_solve_GS(monoPRM, monoCOM, monoMAT)
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
    real(kind=kdouble), allocatable :: R(:), X(:), Dinv(:), ALUtmp(:,:), T(:)

    t1 = monolis_wtime()

    N       = monoMAT%N
    NP      = monoMAT%NP
    NDOF    = monoMAT%NDOF
    NDOF2   = NDOF*NDOF
    NNDOF   = N*NDOF

    allocate(R(NDOF*NP))
    allocate(X(NDOF*NP))
    allocate(Dinv(NDOF2*NP))
    allocate(ALUtmp(NDOF,NDOF))
    allocate(T(NDOF))
    R = 0.0d0
    X = 0.0d0
    Dinv = 0.0d0
    ALUtmp = 0.0d0
    T = 0.0d0

    do i=1,NP
      do j=1,NDOF
        do k=1,NDOF
          ALUtmp(k,j) = monoMAT%D(NDOF2*(i-1)+NDOF*(j-1)+k)
        enddo
      enddo
      do k= 1, NDOF
        ALUtmp(k,k)= 1.0d0/ALUtmp(k,k)
        do l= k+1, NDOF
          ALUtmp(l,k)= ALUtmp(l,k) * ALUtmp(k,k)
          do j= k+1, NDOF
            T(j)= ALUtmp(l,j) - ALUtmp(l,k)*ALUtmp(k,j)
          enddo
          do j= k+1, NDOF
            ALUtmp(l,j)= T(j)
          enddo
        enddo
      enddo
      do j=1,NDOF
        do k=1,NDOF
          Dinv(NDOF2*(i-1)+NDOF*(j-1)+k) = ALUtmp(k,j)
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

    do i=1,NNDOF
      monoMAT%X(i) = X(i)
    enddo

    !call monolis_update_R()

    deallocate(R)
    deallocate(X)
    deallocate(Dinv)
    deallocate(ALUtmp)
    deallocate(T)

    t2 = monolis_wtime()
    tsol = t2 - t1
  end subroutine monolis_solve_GS

end module mod_monolis_solver_GS
