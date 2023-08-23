!> PipeBiCGSTAB 法（前処理なし）モジュール
module mod_monolis_solver_PipeBiCGSTAB_noprec
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_precond
  use mod_monolis_matvec
  use mod_monolis_inner_product
  use mod_monolis_vec_util
  use mod_monolis_converge

  implicit none

contains

  !> @ingroup solver
  !> PipeBiCGSTAB 法（前処理なし）
  subroutine monolis_solver_PipeBiCGSTAB_noprec(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in,out] 行列構造体
    type(monolis_mat), target, intent(inout) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC
    integer(kint) :: N, NP, NDOF, NNDOF
    integer(kint) :: i, iter, iter_RR
    real(kdouble) :: CG(5), RR, RW, RR1, RS, RZ, R2, QY, YY
    real(kdouble) :: alpha, beta, omega, omega1, B2
    real(kdouble) :: tspmv, tdotp, tcomm_spmv, tcomm_dotp
    real(kdouble), allocatable :: R(:), R0(:), W0(:), T(:), S(:), P(:), Z(:), Q(:), Y(:), V(:)
    real(kdouble), pointer :: B(:), X(:)
    logical :: is_converge

    call monolis_std_debug_log_header("monolis_solver_PipeBiCGSTAB_noprec")

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NNDOF = N*NDOF
    X => monoMAT%R%X
    B => monoMAT%R%B
    iter_RR = 200

    tspmv = monoPRM%Rarray(monolis_R_time_spmv)
    tcomm_spmv = monoPRM%Rarray(monolis_R_time_comm_spmv)
    tdotp = monoPRM%Rarray(monolis_R_time_dotp)
    tcomm_dotp = monoPRM%Rarray(monolis_R_time_comm_dotp)

    if(monoPRM%Iarray(monolis_prm_I_is_init_x) == monolis_I_true)then
      X = 0.0d0
    endif

    call monolis_alloc_R_1d(R , NDOF*NP)
    call monolis_alloc_R_1d(R0, NDOF*NP)
    call monolis_alloc_R_1d(W0, NDOF*NP)
    call monolis_alloc_R_1d(T , NDOF*NP)
    call monolis_alloc_R_1d(S , NDOF*NP)
    call monolis_alloc_R_1d(P , NDOF*NP)
    call monolis_alloc_R_1d(Z , NDOF*NP)
    call monolis_alloc_R_1d(Q , NDOF*NP)
    call monolis_alloc_R_1d(Y , NDOF*NP)
    call monolis_alloc_R_1d(V , NDOF*NP)

    call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
    call monolis_set_converge_R(monoCOM, monoMAT, R, B2, is_converge, tdotp, tcomm_dotp)
    if(is_converge) return

    call monolis_vec_copy_R(N, NDOF, R, R0)

    call monolis_matvec_product_main_R(monoCOM, monoMAT, R , W0, tspmv, tcomm_spmv)
    call monolis_matvec_product_main_R(monoCOM, monoMAT, W0, T , tspmv, tcomm_spmv)
    call monolis_inner_product_main_R(monoCOM, N, NDOF, R, R , RR, tdotp, tcomm_dotp)
    call monolis_inner_product_main_R(monoCOM, N, NDOF, R, W0, RW, tdotp, tcomm_dotp)

    alpha = RR / RW
    beta  = 0.0d0
    omega = 0.0d0

    do iter = 1, monoPRM%Iarray(monolis_prm_I_max_iter)
      do i = 1, NNDOF
        P(i) = R (i) + beta *(P(i) - omega*S(i))
        S(i) = W0(i) + beta *(S(i) - omega*Z(i))
        Z(i) = T (i) + beta *(Z(i) - omega*V(i))
        Q(i) = R (i) - alpha* S(i)
        Y(i) = W0(i) - alpha* Z(i)
      enddo

      call monolis_inner_product_main_R(monoCOM, N, NDOF, Q, Y, CG(1), tdotp, tcomm_dotp)
      call monolis_inner_product_main_R(monoCOM, N, NDOF, Y, Y, CG(2), tdotp, tcomm_dotp)

      call monolis_matvec_product_main_R(monoCOM, monoMAT, Z, V, tspmv, tcomm_spmv)

      QY = CG(1)
      YY = CG(2)
      omega1 = QY / YY

      if(mod(iter, iter_RR) == 0)then
        do i = 1, NNDOF
          X(i) = X(i) + alpha*P(i) + omega1*T(i)
        enddo
        call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
        call monolis_matvec_product_main_R(monoCOM, monoMAT, R, W0, tspmv, tcomm_spmv)
      else
        do i = 1, NNDOF
          X (i) = X(i) + alpha * P(i) + omega1*Q(i)
          R (i) = Q(i) - omega1* Y(i)
          W0(i) = Y(i) - omega1*(T(i) - alpha*V(i))
        enddo
      endif

      call monolis_inner_product_main_R(monoCOM, N, NDOF, R0, R, CG(1), tdotp, tcomm_dotp)
      call monolis_inner_product_main_R(monoCOM, N, NDOF, R0, W0,CG(2), tdotp, tcomm_dotp)
      call monolis_inner_product_main_R(monoCOM, N, NDOF, R0, S, CG(3), tdotp, tcomm_dotp)
      call monolis_inner_product_main_R(monoCOM, N, NDOF, R0, Z, CG(4), tdotp, tcomm_dotp)
      call monolis_inner_product_main_R(monoCOM, N, NDOF, R,  R, CG(5), tdotp, tcomm_dotp)

      call monolis_matvec_product_main_R(monoCOM, monoMAT, W0, T, tspmv, tcomm_spmv)

      RR1= CG(1)
      RW = CG(2)
      RS = CG(3)
      RZ = CG(4)
      R2 = CG(5)

      call monolis_check_converge_R(monoPRM, monoCOM, monoMAT, R, B2, iter, is_converge, tspmv, tcomm_spmv)
      if(is_converge) exit

      beta  = (alpha*RR1) / (RR*omega1)
      alpha = RR1 / (RW + beta*RS - beta*omega1*RZ)
      omega = omega1
      RR    = RR1
    enddo

    call monolis_mpi_update_R(monoCOM, NDOF, X, tcomm_spmv)

    call monolis_dealloc_R_1d(R )
    call monolis_dealloc_R_1d(R0)
    call monolis_dealloc_R_1d(W0)
    call monolis_dealloc_R_1d(T )
    call monolis_dealloc_R_1d(S )
    call monolis_dealloc_R_1d(P )
    call monolis_dealloc_R_1d(Z )
    call monolis_dealloc_R_1d(Q )
    call monolis_dealloc_R_1d(Y )
    call monolis_dealloc_R_1d(V )
  end subroutine monolis_solver_PipeBiCGSTAB_noprec
end module mod_monolis_solver_PipeBiCGSTAB_noprec
