!> Deflated CG 法モジュール
module mod_monolis_solver_DeflatedCG
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
  !> Deflated CG 法
  subroutine monolis_solver_DeflatedCG(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in,out] 行列構造体
    type(monolis_mat), target, intent(inout) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC
    type(monolis_prm) :: monoPRM_deflated_eq
    type(monolis_com) :: monoCOM_deflated_eq
    type(monolis_mat) :: monoMAT_deflated_eq
    integer(kint) :: N, NP, NDOF, NNDOF, NPNDOF, M, M_neib
    integer(kint) :: i, iter, iter_RR
    real(kdouble) :: alpha, beta, rho, rho1, omega, B2
    real(kdouble) :: tspmv, tdotp, tcomm_spmv, tcomm_dotp
    logical :: is_converge
    real(kdouble), allocatable :: R(:), Z(:), Q(:), P(:)
    real(kdouble), allocatable :: W(:,:), WtA(:,:), WtW(:,:), g(:), X0(:)
    real(kdouble), pointer :: B(:), X(:)

    call monolis_std_debug_log_header("monolis_solver_CG")

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NNDOF = N*NDOF
    NPNDOF= NP*NDOF
    X => monoMAT%R%X
    B => monoMAT%R%B
    iter_RR = 100
    M = monoPRM%Iarray(monolis_prm_I_n_local_deflation_mode)
    M_neib = M*(monoCOM%recv_n_neib + 1)

    tspmv = monoPRM%Rarray(monolis_R_time_spmv)
    tcomm_spmv = monoPRM%Rarray(monolis_R_time_comm_spmv)
    tdotp = monoPRM%Rarray(monolis_R_time_dotp)
    tcomm_dotp = monoPRM%Rarray(monolis_R_time_comm_dotp)

    if(monoPRM%Iarray(monolis_prm_I_is_init_x) == monolis_I_true)then
      X = 0.0d0
    endif

    call monolis_alloc_R_1d(R, NPNDOF)
    call monolis_alloc_R_1d(Z, NPNDOF)
    call monolis_alloc_R_1d(Q, NPNDOF)
    call monolis_alloc_R_1d(P, NPNDOF)
    call monolis_alloc_R_1d(X0, NPNDOF)

    call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
    call monolis_set_converge_R(monoCOM, monoMAT, R, B2, is_converge, tdotp, tcomm_dotp)
    if(is_converge) return

    if(M > 0)then
      call deflatedCG_initialize()
      call deflatedCG_pre()
    endif

    !call monolis_set_converge_R(monoCOM, monoMAT, R, B2, is_converge, tdotp, tcomm_dotp)

    do iter = 1, monoPRM%Iarray(monolis_prm_I_max_iter)
      call monolis_precond_apply_R(monoPRM, monoCOM, monoMAT, monoPREC, R, Z)
      call monolis_inner_product_main_R(monoCOM, N, NDOF, R, Z, rho, tdotp, tcomm_dotp)

      if(1 < iter)then
        beta = rho/rho1
        call monolis_vec_AXPBY_R(N, NDOF, beta, P, 1.0d0, Z, P)
      else
        call monolis_vec_copy_R(N, NDOF, Z, P)
      endif

      if(M > 0)then
        call deflatedCG_omit()
      endif

      call monolis_matvec_product_main_R(monoCOM, monoMAT, P, Q, tspmv, tcomm_spmv)
      call monolis_inner_product_main_R(monoCOM, N, NDOF, P, Q, omega, tdotp, tcomm_dotp)
      alpha = rho/omega

      call monolis_vec_AXPBY_R(N, NDOF, alpha, P, 1.0d0, X, X)
      call monolis_vec_AXPBY_R(N, NDOF,-alpha, Q, 1.0d0, R, R)

      if(mod(iter, iter_RR) == 0)then
        if(M > 0)then
          call deflatedCG_residual_replacement()
        endif
      endif

      call monolis_check_converge_R(monoPRM, monoCOM, monoMAT, R, B2, iter, is_converge, tdotp, tcomm_dotp)
      if(is_converge) exit

      rho1 = rho
    enddo

    if(M > 0)then
      call deflatedCG_post()
    endif

    call monolis_mpi_update_R(monoCOM, NDOF, X, tcomm_spmv)

    monoPRM%Rarray(monolis_R_time_spmv) = tspmv
    monoPRM%Rarray(monolis_R_time_comm_spmv) = tcomm_spmv
    monoPRM%Rarray(monolis_R_time_dotp) = tdotp
    monoPRM%Rarray(monolis_R_time_comm_dotp) = tcomm_dotp

    call monolis_dealloc_R_1d(R)
    call monolis_dealloc_R_1d(Z)
    call monolis_dealloc_R_1d(Q)
    call monolis_dealloc_R_1d(P)
    call monolis_dealloc_R_1d(X0)

    if(M > 0)then
      call deflatedCG_finalize()
    endif
  end subroutine monolis_solver_DeflatedCG

  !> @ingroup dev_solver
  !> Deflated CG 法
  subroutine deflatedCG_initialize()
    implicit none

!    call monolis_prm_initialize(monoPRM_deflated_eq)
!    call monolis_mat_initialize(monoMAT_deflated_eq)
!    call monolis_com_initialize(monoCOM_deflated_eq)
!
!    call monolis_alloc_R_2d(W, NPNDOF, M_neib)
!    call monolis_alloc_R_2d(WtA, M_neib, NNDOF)
!    call monolis_alloc_R_2d(WtA, M, M)
!    call monolis_alloc_R_1d(g, M_neib)
!    call monolis_alloc_R_1d(IPV_R, M)

    !call deflatedCG_set_deflation_mode()
    !call deflatedCG_set_defeq_initialize()
    !call deflatedCG_set_defeq_mat_gen()
  end subroutine deflatedCG_initialize

  !> @ingroup dev_solver
  !> Deflated CG 法
  subroutine deflatedCG_pre()
    implicit none

!    N     = monoMAT%N
!    NP    = monoMAT%NP
!    NDOF  = monoMAT%NDOF
!    NNDOF = N*NDOF
!
!    allocate(WTR(M_neib), source = 0.0d0)
!    allocate(Wg(NNDOF), source = 0.0d0)
!
!    monoMAT_deflated_eq%B = 0.0d0
!    call monolis_dense_matvec_local(M, NNDOF, transpose(W(1:NNDOF,1:M)), R(1:NNDOF), monoMAT_deflated_eq%B(1:M), monoPRM%tdemv)
!
!    call interface_monolis_solve_(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq)
!    call monolis_DeflatedCG_timer_totalling(monoPRM, monoPRM_deflated_eq)
!
!    g = monoMAT_deflated_eq%X
!
!    call monolis_dense_matvec_local(NNDOF, M, W, g, Wg, monoPRM%tdemv)
!    call monolis_vec_AXPY(N, NDOF, 1.0d0, X, Wg, X)
!    call monolis_vec_copy_R(N, NDOF, X, X0)
!
!    call monolis_residual(monoCOM, monoMAT, X, B, R, time, time)
!
!    call monolis_dense_matvec_local(M, NNDOF, transpose(W(1:NNDOF,1:M)), R(1:NNDOF), monoMAT_deflated_eq%B(1:M), monoPRM%tdemv)
!
!    call interface_monolis_solve_(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq)
!    call monolis_DeflatedCG_timer_totalling(monoPRM, monoPRM_deflated_eq)
!
!    g = monoMAT_deflated_eq%X
!
!    call monolis_dense_matmul_local(M, NNDOF, M, transpose(W(1:NNDOF,1:M)), W(1:NNDOF,1:M), WTW, monoPRM%tdemv)
!
!    call monolis_lapack_LU_fact(M, WTW, IPV_R)
  end subroutine deflatedCG_pre

  !> @ingroup dev_solver
  !> Deflated CG 法
  subroutine deflatedCG_omit()
    implicit none

!    N     = monoMAT%N
!    NDOF  = monoMAT%NDOF
!    NNDOF = N*NDOF
!
!    allocate(WTAZ(M_neib), source = 0.0d0)
!    allocate(WLinvWTAZ(NNDOF), source = 0.0d0)
!
!    call monolis_dense_matvec_local(M_neib, NNDOF, WTA(1:M_neib,1:NNDOF), Z(1:NNDOF), WTAZ, monoPRM%tdemv)
!
!    monoMAT_deflated_eq%B(1:M) = WTAZ(1:M)
!
!    call monolis_update_R_reverse(monoCOM_deflated_eq, M, WTAZ, monoPRM%tcomm_urev)
!
!    monoMAT_deflated_eq%B(1:M) = monoMAT_deflated_eq%B(1:M) + WTAZ(1:M)
!
!    call interface_monolis_solve_(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq)
!    call monolis_DeflatedCG_timer_totalling(monoPRM, monoPRM_deflated_eq)
!
!    call monolis_dense_matvec_local(NNDOF, M, W(1:NNDOF,1:M), monoMAT_deflated_eq%X(1:M), WLinvWTAZ, monoPRM%tdemv)
!    call monolis_vec_AXPY(N, NDOF, -1.0d0, WLinvWTAZ, P, P)
  end subroutine deflatedCG_omit

  !> @ingroup dev_solver
  !> Deflated CG 法
  subroutine deflatedCG_residual_replacement()
    implicit none

!    N     = monoMAT%N
!    NDOF  = monoMAT%NDOF
!    NNDOF = N*NDOF
!
!    allocate(WWTWinvWTR(NNDOF), source = 0.0d0)
!
!    call monolis_dense_matvec_local(M, NNDOF, transpose(W(1:NNDOF,1:M)), R(1:NNDOF), WTR, monoPRM%tdemv)
!
!    call monolis_lapack_LU_solve(M, WTW, WTR, IPV_R, WTR)
!
!    call monolis_dense_matvec_local(NNDOF, M, W, WTR, WWTWinvWTR, monoPRM%tdemv)
!
!    call monolis_vec_AXPY(N, NDOF, -1.0d0, WWTWinvWTR, R, R)
  end subroutine deflatedCG_residual_replacement

  !> @ingroup dev_solver
  !> Deflated CG 法
  subroutine deflatedCG_post()
    implicit none

!    N     = monoMAT%N
!    NP    = monoMAT%NP
!    NDOF  = monoMAT%NDOF
!    NNDOF = N*NDOF
!    NPNDOF= NP*NDOF
!
!    allocate(AX(NPNDOF), source = 0.0d0)
!    allocate(h(M_neib), source = 0.0d0)
!    allocate(WgMh(NNDOF), source = 0.0d0)
!
!    call monolis_matvec(monoCOM, monoMAT, X - X0, AX, time, time)
!
!    call monolis_dense_matvec_local(M, NNDOF, transpose(W(1:NNDOF,1:M)), AX(1:NNDOF), monoMAT_deflated_eq%B(1:M), monoPRM%tdemv)
!
!    call interface_monolis_solve_(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq)
!    call monolis_DeflatedCG_timer_totalling(monoPRM, monoPRM_deflated_eq)
!
!    h = monoMAT_deflated_eq%X
!
!    call monolis_dense_matvec_local(NNDOF, M_neib, W, g - h, WgMh, monoPRM%tdemv)
!
!    call monolis_vec_AXPY(N, NDOF, 1.0d0, WgMh, X, X)
!
!    call monolis_update_R(monoCOM, NDOF, X, time)
  end subroutine deflatedCG_post

  !> @ingroup dev_solver
  !> Deflated CG 法
  subroutine deflatedCG_finalize()
    implicit none

!    deallocate(W)
!    deallocate(WTA)
!    deallocate(g)
!    deallocate(WTW)
!    deallocate(IPV_R)
!
!    call monolis_prm_finalize(monoPRM_deflated_eq)
!    call monolis_mat_finalize(monoMAT_deflated_eq)
!    call monolis_com_finalize(monoCOM_deflated_eq)
  end subroutine deflatedCG_finalize
end module mod_monolis_solver_DeflatedCG
