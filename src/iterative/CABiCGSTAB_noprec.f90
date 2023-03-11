module mod_monolis_solver_CABiCGSTAB_noprec
  use mod_monolis_prm
  use mod_monolis_mat
  use mod_monolis_precond
  use mod_monolis_matvec
  use mod_monolis_inner_product
  use mod_monolis_converge

  implicit none

contains

  subroutine monolis_solver_CABiCGSTAB_noprec(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: N, NP, NDOF, NNDOF
    integer(kind=kint) :: i, iter, iter_RR
    real(kind=kdouble) :: alpha, beta, rho, CG(5), omega, B2, R2
    real(kind=kdouble), allocatable :: R0(:), R(:), V(:), S(:), P(:), Q(:), Y(:), Z(:)
    real(kind=kdouble), pointer :: B(:), X(:)
    logical :: is_converge

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NNDOF = N*NDOF
    X => monoMAT%X
    B => monoMAT%B
    iter_RR = 50

    if(monoPRM%is_init_x) X = 0.0d0

    allocate(R0(NDOF*NP), source = 1.0d0)
    allocate(R (NDOF*NP), source = 1.0d0)
    allocate(V (NDOF*NP), source = 1.0d0)
    allocate(S (NDOF*NP), source = 1.0d0)
    allocate(P (NDOF*NP), source = 1.0d0)
    allocate(Q (NDOF*NP), source = 1.0d0)
    allocate(Y (NDOF*NP), source = 1.0d0)
    allocate(Z (NDOF*NP), source = 1.0d0)

    call monolis_residual(monoCOM, monoMAT, X, B, R0, monoPRM%tspmv, monoPRM%tcomm_spmv)
    call monolis_set_converge(monoPRM, monoCOM, monoMAT, R0, B2, is_converge, monoPRM%tdotp, monoPRM%tcomm_dotp)
    if(is_converge) return

    call monolis_vec_copy_R(N, NDOF, R0, R)

    call monolis_matvec(monoCOM, monoMAT, R0, V, monoPRM%tspmv, monoPRM%tcomm_spmv)
    call monolis_inner_product_R(monoCOM, N, NDOF, R0, R0, CG(1), monoPRM%tdotp, monoPRM%tcomm_dotp)
    call monolis_inner_product_R(monoCOM, N, NDOF, R0, V , CG(2), monoPRM%tdotp, monoPRM%tcomm_dotp)

    alpha = CG(1)/CG(2)
    beta  = 0.0d0
    omega = 0.0d0
    rho   = CG(1)

    call monolis_inner_product_R(monoCOM, N, NDOF, B, B, B2, monoPRM%tdotp, monoPRM%tcomm_dotp)

    do iter = 1, monoPRM%maxiter
      do i = 1, NNDOF
        P(i) = R(i) + beta * (P(i) - omega * S(i))
        S(i) = V(i) + beta * (S(i) - omega * Z(i))
      enddo

      call monolis_matvec(monoCOM, monoMAT, S, Z, monoPRM%tspmv, monoPRM%tcomm_spmv)

      do i = 1, NNDOF
        Q(i) = R(i) - alpha * S(i)
        Y(i) = V(i) - alpha * Z(i)
      enddo

      call monolis_inner_product_R(monoCOM, N, NDOF, Q, Y, CG(1), monoPRM%tdotp, monoPRM%tcomm_dotp)
      call monolis_inner_product_R(monoCOM, N, NDOF, Y, Y, CG(2), monoPRM%tdotp, monoPRM%tcomm_dotp)

      omega = CG(1)/CG(2)

      do i = 1, NNDOF
        X(i) = X(i) + alpha*P(i) + omega*Q(i)
        R(i) = Q(i) - omega*Y(i)
      enddo

      call monolis_matvec(monoCOM, monoMAT, R, V, monoPRM%tspmv, monoPRM%tcomm_spmv)

      call monolis_inner_product_R(monoCOM, N, NDOF, R0, R, CG(1), monoPRM%tdotp, monoPRM%tcomm_dotp)
      call monolis_inner_product_R(monoCOM, N, NDOF, R0, V, CG(2), monoPRM%tdotp, monoPRM%tcomm_dotp)
      call monolis_inner_product_R(monoCOM, N, NDOF, R0, S, CG(3), monoPRM%tdotp, monoPRM%tcomm_dotp)
      call monolis_inner_product_R(monoCOM, N, NDOF, R0, Z, CG(4), monoPRM%tdotp, monoPRM%tcomm_dotp)
      call monolis_inner_product_R(monoCOM, N, NDOF, R , R, CG(5), monoPRM%tdotp, monoPRM%tcomm_dotp)

      beta  = (alpha/omega) * CG(1) / rho
      alpha = CG(1) / (CG(2) + beta * CG(3) - beta * omega * CG(4))
      R2    = CG(5)

      call monolis_check_converge_2(monoPRM, monoCOM, monoMAT, R2, B2, iter, is_converge)
      if(is_converge) exit

      rho = CG(1)
    enddo

    call monolis_update_R(monoCOM, NDOF, X, monoPRM%tcomm_spmv)

    deallocate(R0)
    deallocate(R )
    deallocate(V )
    deallocate(S )
    deallocate(P )
    deallocate(Q )
    deallocate(Y )
    deallocate(Z )
  end subroutine monolis_solver_CABiCGSTAB_noprec

end module mod_monolis_solver_CABiCGSTAB_noprec
