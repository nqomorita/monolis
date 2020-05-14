module mod_monolis_solver_PipeBiCGSTAB_noprec
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

  subroutine monolis_solver_PipeBiCGSTAB_noprec(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: N, NP, NDOF, NNDOF
    integer(kind=kint) :: i, iter, iter_RR
    real(kind=kdouble) :: t1, t2, tsol, tcomm, CG(5), RR, RW, RR1, RS, RZ, R2, QY, YY
    real(kind=kdouble) :: alpha, beta, omega, omega1, B2
    real(kind=kdouble), allocatable :: R(:), R0(:), W0(:), T(:), S(:), P(:), Z(:), Q(:), Y(:), V(:)
    real(kind=kdouble), pointer :: B(:), X(:)
    logical :: is_converge

    t1 = monolis_get_time()

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NNDOF = N*NDOF
    X => monoMAT%X
    B => monoMAT%B
    iter_RR = 50

    if(monoPRM%is_init_x) X = 0.0d0

    allocate(R (NDOF*NP)); R  = 0.0d0
    allocate(R0(NDOF*NP)); R0 = 0.0d0
    allocate(W0(NDOF*NP)); W0 = 0.0d0
    allocate(T (NDOF*NP)); T  = 0.0d0
    allocate(S (NDOF*NP)); S  = 0.0d0
    allocate(P (NDOF*NP)); P  = 0.0d0
    allocate(Z (NDOF*NP)); Z  = 0.0d0
    allocate(Q (NDOF*NP)); Q  = 0.0d0
    allocate(Y (NDOF*NP)); Y  = 0.0d0
    allocate(V (NDOF*NP)); V  = 0.0d0

    call monolis_set_converge(monoPRM, monoCOM, monoMAT, B, B2, is_converge, tcomm)
    if(is_converge) return
    call monolis_residual(monoCOM, monoMAT, X, B, R, tcomm)

    call monolis_vec_copy_R(N, NDOF, R, R0)

    call monolis_matvec(monoCOM, monoMAT, R , W0, tcomm)
    call monolis_matvec(monoCOM, monoMAT, W0, T , tcomm)
    call monolis_inner_product_R(monoCOM, N, NDOF, R, R , RR, tcomm)
    call monolis_inner_product_R(monoCOM, N, NDOF, R, W0, RW, tcomm)

    alpha = RR / RW
    beta  = 0.0d0
    omega = 0.0d0

    do iter = 1, monoPRM%maxiter
      do i = 1, NNDOF
        P(i) = R (i) + beta *(P(i) - omega*S(i))
        S(i) = W0(i) + beta *(S(i) - omega*Z(i))
        Z(i) = T (i) + beta *(Z(i) - omega*V(i))
        Q(i) = R (i) - alpha* S(i)
        Y(i) = W0(i) - alpha* Z(i)
      enddo

      call monolis_inner_product_R(monoCOM, N, NDOF, Q, Y, CG(1), tcomm)
      call monolis_inner_product_R(monoCOM, N, NDOF, Y, Y, CG(2), tcomm)

      call monolis_matvec(monoCOM, monoMAT, Z, V, tcomm)

      QY = CG(1)
      YY = CG(2)
      omega1 = QY / YY

      if(mod(iter, iter_RR) == 0)then
        do i = 1, NNDOF
          X(i) = X(i) + alpha*P(i) + omega1*T(i)
        enddo
        call monolis_residual(monoCOM, monoMAT, X, B, R, tcomm)
        call monolis_matvec(monoCOM, monoMAT, R, W0, tcomm)
      else
        do i = 1, NNDOF
          X (i) = X(i) + alpha * P(i) + omega1*Q(i)
          R (i) = Q(i) - omega1* Y(i)
          W0(i) = Y(i) - omega1*(T(i) - alpha*V(i))
        enddo
      endif

      call monolis_inner_product_R(monoCOM, N, NDOF, R0, R, CG(1), tcomm)
      call monolis_inner_product_R(monoCOM, N, NDOF, R0, W0,CG(2), tcomm)
      call monolis_inner_product_R(monoCOM, N, NDOF, R0, S, CG(3), tcomm)
      call monolis_inner_product_R(monoCOM, N, NDOF, R0, Z, CG(4), tcomm)
      call monolis_inner_product_R(monoCOM, N, NDOF, R,  R, CG(5), tcomm)

      call monolis_matvec(monoCOM, monoMAT, W0, T, tcomm)

      RR1= CG(1)
      RW = CG(2)
      RS = CG(3)
      RZ = CG(4)
      R2 = CG(5)

      call monolis_check_converge_2(monoPRM, monoCOM, monoMAT, R2, B2, iter, is_converge, tcomm)
      if(is_converge) exit

      beta  = (alpha*RR1) / (RR*omega1)
      alpha = RR1 / (RW + beta*RS - beta*omega1*RZ)
      omega = omega1
      RR    = RR1
    enddo

    call monolis_update_R(monoCOM, NDOF, X, tcomm)

    deallocate(R )
    deallocate(R0)
    deallocate(W0)
    deallocate(T )
    deallocate(S )
    deallocate(P )
    deallocate(Z )
    deallocate(Q )
    deallocate(Y )
    deallocate(V )

    t2 = monolis_get_time()
    tsol = t2 - t1
  end subroutine monolis_solver_PipeBiCGSTAB_noprec
end module mod_monolis_solver_PipeBiCGSTAB_noprec
