module mod_monolis_solver_IR

  implicit none

contains

  subroutine monolis_solver_IR(monoPRM, monoCOM, monoMAT)
    use mod_monolis_prm
    use mod_monolis_com
    use mod_monolis_mat
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: N, NP, NDOF, NNDOF
    integer(kind=kint) :: i, iter
    real(kind=kdouble) :: tol, resid, R2, B2, D2
    real(kind=kdouble) :: t1, t2, tset, tsol
    real(kind=kdouble), allocatable :: R(:), D(:), X(:), T(:)

    t1 = monolis_wtime()

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NNDOF = N*NDOF

    allocate(R(NDOF*NP))
    allocate(D(NDOF*NP))
    allocate(X(NDOF*NP))
    allocate(T(NDOF*NP))
    R = 0.0d0
    D = 0.0d0
    X = 0.0d0
    T = 0.0d0

    !call monolis_inner_product_R()
    !call monolis_precond_setup()

    do i=1,NNDOF
      R(i) = monoMAT%B(i)
    enddo

    do iter=1, monoPRM%maxiter
      !call monolis_precond_apply()

      do i=1,NNDOF
        X(i) = X(i) + D(i)
      enddo

      !call monolis_residual()
      !call monolis_inner_product_R()
      resid = dsqrt(R2/B2)

      if(monoCOM%myrank == 0) write (*,"(i7, 1pe16.6)") iter, resid
      if(resid <= tol) exit
    enddo

    do i=1,NNDOF
      monoMAT%X(i) = X(i)
    enddo

    !call hecmw_update_R()

    deallocate(R)
    deallocate(D)
    deallocate(X)
    deallocate(T)

    t2 = monolis_wtime()
    tsol = t2 - t1
  end subroutine monolis_solver_IR

end module mod_monolis_solver_IR
