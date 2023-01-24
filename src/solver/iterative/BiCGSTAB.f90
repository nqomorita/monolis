module mod_monolis_solver_BiCGSTAB
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

  subroutine monolis_solver_BiCGSTAB(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: N, NP, NDOF, NNDOF
    integer(kint) :: i, iter, iter_RR
    real(kdouble) :: alpha, beta, rho, rho1, c2, omega
    real(kdouble) :: B2, CG(2)
    real(kdouble), allocatable :: R(:), RT(:), P(:), PT(:), S(:), ST(:), T(:), V(:)
    real(kdouble), pointer :: B(:), X(:)
    logical :: is_converge

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NNDOF = N*NDOF
    X => monoMAT%X
    B => monoMAT%B
    iter_RR = 200
    omega = 0.0d0

    if(monoPRM%is_init_x) X = 0.0d0

    allocate(R (NDOF*NP), source = 0.0d0)
    allocate(RT(NDOF*NP), source = 0.0d0)
    allocate(P (NDOF*NP), source = 0.0d0)
    allocate(PT(NDOF*NP), source = 0.0d0)
    allocate(S (NDOF*NP), source = 0.0d0)
    allocate(ST(NDOF*NP), source = 0.0d0)
    allocate(T (NDOF*NP), source = 0.0d0)
    allocate(V (NDOF*NP), source = 0.0d0)

    call monolis_residual(monoCOM, monoMAT, X, B, R, monoPRM%tspmv, monoPRM%tcomm_spmv)
    call monolis_set_converge(monoPRM, monoCOM, monoMAT, R, B2, is_converge, monoPRM%tdotp, monoPRM%tcomm_dotp)
    if(is_converge) return

    call monolis_vec_copy_R(N, NDOF, R, RT)

    do iter = 1, monoPRM%maxiter
      call monolis_inner_product_R(monoCOM, N, NDOF, R, RT, rho, monoPRM%tdotp, monoPRM%tcomm_dotp)

      if(1 < iter)then
        beta = (rho/rho1) * (alpha/omega)
        do i = 1, NNDOF
          P(i) = R(i) + beta * (P(i) - omega * V(i))
        enddo
      else
        call monolis_vec_copy_R(N, NDOF, R, P)
      endif

      call monolis_precond_apply(monoPRM, monoCOM, monoMAT, P, PT)
      call monolis_matvec(monoCOM, monoMAT, PT, V, monoPRM%tspmv, monoPRM%tcomm_spmv)
      call monolis_inner_product_R(monoCOM, N, NDOF, RT, V, c2, monoPRM%tdotp, monoPRM%tcomm_dotp)

      alpha = rho / c2
      call monolis_vec_AXPY(N, NDOF, -alpha, V, R, S)

      call monolis_precond_apply(monoPRM, monoCOM, monoMAT, S, ST)
      call monolis_matvec(monoCOM, monoMAT, ST, T, monoPRM%tspmv, monoPRM%tcomm_spmv)

      call monolis_inner_product_R_local(monoCOM, N, NDOF, T, S, CG(1))
      call monolis_inner_product_R_local(monoCOM, N, NDOF, T, T, CG(2))
      call monolis_allreduce_R(2, CG, monolis_sum, monoCOM%comm)

      if(CG(2) /= 0.0d0)then
        omega = CG(1) / CG(2)
      else
        omega = 0.0d0
      endif

      do i = 1, NNDOF
        X(i) = X(i) + alpha*PT(i) + omega*ST(i)
      enddo

      if(mod(iter, iter_RR) == 0)then
        call monolis_residual(monoCOM, monoMAT, X, B, R, monoPRM%tspmv, monoPRM%tcomm_spmv)
      else
        call monolis_vec_AXPY(N, NDOF, -omega, T, S, R)
      endif

      call monolis_check_converge(monoPRM, monoCOM, monoMAT, R, B2, iter, is_converge, monoPRM%tdotp, monoPRM%tcomm_dotp)
      if(is_converge) exit

      rho1 = rho
    enddo

    call monolis_update_R(monoCOM, NDOF, X, monoPRM%tcomm_spmv)

    deallocate(R )
    deallocate(RT)
    deallocate(P )
    deallocate(PT)
    deallocate(S )
    deallocate(ST)
    deallocate(T )
    deallocate(V )
  end subroutine monolis_solver_BiCGSTAB

end module mod_monolis_solver_BiCGSTAB
