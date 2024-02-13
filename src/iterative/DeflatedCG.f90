!> Deflated CG 法モジュール
module mod_monolis_solver_DeflatedCG
  use mod_monolis_def_solver
  use mod_monolis_def_mat
  use mod_monolis_utils_define_com_init
  use mod_monolis_def_solver_util
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_precond
  use mod_monolis_matvec
  use mod_monolis_matalg_dense
  use mod_monolis_inner_product
  use mod_monolis_vec_util
  use mod_monolis_converge
  use mod_monolis_spmat_handler
  use mod_monolis_spmat_nonzero_pattern_util
  use mod_monolis_solver_DeflatedCG_util

  implicit none

contains

  !> @ingroup solver
  !> Deflated CG 法
  subroutine monolis_solver_DeflatedCG1(monoPRM, monoCOM, monoMAT, monoPREC)
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
    type(monolis_mat) :: monoPRE_deflated_eq
    integer(kint) :: N, NP, NDOF, NNDOF, NPNDOF, M, M_neib
    integer(kint) :: i, iter, iter_RR
    real(kdouble) :: alpha, beta, rho, rho1, omega, B2
    real(kdouble) :: tspmv, tdotp, tcomm_spmv, tcomm_dotp, tdemv
    logical :: is_converge
    integer(kint), allocatable :: IPV_R(:)
    real(kdouble), allocatable :: R(:), Z(:), Q(:), P(:), X0(:), Qb(:), PtX(:)
    real(kdouble), allocatable :: W(:,:), AW(:,:), WtA(:,:), WtW(:,:)
    real(kdouble), pointer :: B(:), X(:)

    call monolis_std_debug_log_header("monolis_solver_DeflatedCG1")

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
    call monolis_alloc_R_1d(Qb, NPNDOF)
    call monolis_alloc_R_1d(PtX, NPNDOF)

    call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
    call monolis_set_converge_R(monoCOM, monoMAT, R, B2, is_converge, tdotp, tcomm_dotp)
    if(is_converge) return

    if(M > 0)then
      call deflatedCG_set_deflation_mode(monoPRM, monoCOM, monoMAT, NPNDOF, M, M_neib, W)
      call deflatedCG_E_initialize(monoCOM, monoPRM, monoMAT, monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, &
        & M, M_neib, NNDOF, W, AW, WtA)
      call deflatedCG_residual_replacement_initialize(M, NNDOF, W, WtW, IPV_R)

      call deflatedCG_P(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, &
        & M, M_neib, NNDOF, W, AW, R, R, tdemv)
    endif

    call monolis_precond_apply_R(monoPRM, monoCOM, monoMAT, monoPREC, R, Z)
    call monolis_vec_copy_R(N, NDOF, Z, P)

    call monolis_inner_product_main_R(monoCOM, N, NDOF, R, Z, rho, tdotp, tcomm_dotp)

    do iter = 1, monoPRM%Iarray(monolis_prm_I_max_iter)
      call monolis_matvec_product_main_R(monoCOM, monoMAT, P, Q, tspmv, tcomm_spmv)

      call deflatedCG_P(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, &
        & M, M_neib, NNDOF, W, AW, Q, Q, tdemv)

      call monolis_inner_product_main_R(monoCOM, N, NDOF, P, Q, omega, tdotp, tcomm_dotp)

      alpha = rho/omega

      call monolis_vec_AXPBY_R(N, NDOF, alpha, P, 1.0d0, X, X)
      call monolis_vec_AXPBY_R(N, NDOF,-alpha, Q, 1.0d0, R, R)

      if(mod(iter, iter_RR) == 0)then
        if(M > 0)then
          call deflatedCG_residual_replacement(M, N, NDOF, W, R, WtW, IPV_R)
        endif
      endif

      call monolis_check_converge_R(monoPRM, monoCOM, monoMAT, R, B2, iter, is_converge, tdotp, tcomm_dotp)
      if(is_converge) exit

      rho1 = rho

      call monolis_precond_apply_R(monoPRM, monoCOM, monoMAT, monoPREC, R, Z)
      call monolis_inner_product_main_R(monoCOM, N, NDOF, R, Z, rho, tdotp, tcomm_dotp)

      beta = rho/rho1

      call monolis_vec_AXPBY_R(N, NDOF, beta, P, 1.0d0, Z, P)
    enddo

    call deflatedCG_Q(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, &
      & M, NNDOF, W, B, Qb, tdemv)
    call deflatedCG_Pt(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, &
      & M, M_neib, NNDOF, W, WtA, X, PtX, tdemv)
    call monolis_vec_AXPBY_R(N, NDOF, 1.0d0, Qb, 1.0d0, PtX, X)

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
    call monolis_dealloc_R_1d(Qb)
    call monolis_dealloc_R_1d(PtX)

    if(M > 0)then
      call deflatedCG_finalize(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, &
        & IPV_R, W, AW, WtA, WtW)
    endif
  end subroutine monolis_solver_DeflatedCG1

  subroutine monolis_solver_DeflatedCG2(monoPRM, monoCOM, monoMAT, monoPREC)
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
    type(monolis_mat) :: monoPRE_deflated_eq
    integer(kint) :: N, NP, NDOF, NNDOF, NPNDOF, M, M_neib
    integer(kint) :: i, iter, iter_RR
    real(kdouble) :: alpha, beta, rho, rho1, omega, B2
    real(kdouble) :: tspmv, tdotp, tcomm_spmv, tcomm_dotp, tdemv
    logical :: is_converge
    integer(kint), allocatable :: IPV_R(:)
    real(kdouble), allocatable :: R(:), Z(:), Q(:), P(:), X0(:), Qb(:), PtX(:)
    real(kdouble), allocatable :: W(:,:), AW(:,:), WtA(:,:), WtW(:,:)
    real(kdouble), pointer :: B(:), X(:)

    call monolis_std_debug_log_header("monolis_solver_DeflatedCG2")

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
    call monolis_alloc_R_1d(Qb, NPNDOF)
    call monolis_alloc_R_1d(PtX, NPNDOF)

    call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
    call monolis_set_converge_R(monoCOM, monoMAT, R, B2, is_converge, tdotp, tcomm_dotp)
    if(is_converge) return

    if(M > 0)then
      call deflatedCG_set_deflation_mode(monoPRM, monoCOM, monoMAT, NPNDOF, M, M_neib, W)
      call deflatedCG_E_initialize(monoCOM, monoPRM, monoMAT, monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, &
        & M, M_neib, NNDOF, W, AW, WtA)

      call deflatedCG_Q(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, &
        & M, NNDOF, W, B, Qb, tdemv)

      if(monoPRM%Iarray(monolis_prm_I_is_init_x) /= monolis_I_true)then
        call deflatedCG_Pt(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, &
          & M, M_neib, NNDOF, W, WtA, X, PtX, tdemv)
      endif

      call monolis_vec_AXPBY_R(N, NDOF, 1.0d0, Qb, 1.0d0, PtX, X)
      call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)

      call deflatedCG_residual_replacement_initialize(M, NNDOF, W, WtW, IPV_R)
    endif

    call monolis_inner_product_main_R(monoCOM, N, NDOF, R, R, rho, tdotp, tcomm_dotp)

    if(rho/B2 < monoPRM%Rarray(monolis_prm_R_tol))then
      monoPRM%Iarray(monolis_prm_I_cur_iter) = 0
      monoPRM%Rarray(monolis_prm_R_cur_resid) = rho/B2
      return
    endif

    call monolis_precond_apply_R(monoPRM, monoCOM, monoMAT, monoPREC, R, Z)

    call deflatedCG_Pt(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, &
      & M, M_neib, NNDOF, W, WtA, Z, P, tdemv)

    call monolis_inner_product_main_R(monoCOM, N, NDOF, R, Z, rho, tdotp, tcomm_dotp)

    do iter = 1, monoPRM%Iarray(monolis_prm_I_max_iter)
      call monolis_matvec_product_main_R(monoCOM, monoMAT, P, Q, tspmv, tcomm_spmv)
      call monolis_inner_product_main_R(monoCOM, N, NDOF, P, Q, omega, tdotp, tcomm_dotp)

      alpha = rho/omega

      call monolis_vec_AXPBY_R(N, NDOF, alpha, P, 1.0d0, X, X)
      call monolis_vec_AXPBY_R(N, NDOF,-alpha, Q, 1.0d0, R, R)

      if(mod(iter, iter_RR) == 0)then
        if(M > 0)then
          call deflatedCG_residual_replacement(M, N, NDOF, W, R, WtW, IPV_R)
        endif
      endif

      call monolis_check_converge_R(monoPRM, monoCOM, monoMAT, R, B2, iter, is_converge, tdotp, tcomm_dotp)
      if(is_converge) exit

      rho1 = rho

      call monolis_precond_apply_R(monoPRM, monoCOM, monoMAT, monoPREC, R, Z)
      call monolis_inner_product_main_R(monoCOM, N, NDOF, R, Z, rho, tdotp, tcomm_dotp)

      beta = rho/rho1

      call deflatedCG_Pt(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, &
        & M, M_neib, NNDOF, W, WtA, Z, PtX, tdemv)

      call monolis_vec_AXPBY_R(N, NDOF, beta, P, 1.0d0, PtX, P)
    enddo

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
    call monolis_dealloc_R_1d(Qb)
    call monolis_dealloc_R_1d(PtX)

    if(M > 0)then
      call deflatedCG_finalize(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, &
        & IPV_R, W, AW, WtA, WtW)
    endif
  end subroutine monolis_solver_DeflatedCG2

  subroutine monolis_solver_ADeflatedCG2(monoPRM, monoCOM, monoMAT, monoPREC)
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
    type(monolis_mat) :: monoPRE_deflated_eq
    integer(kint) :: N, NP, NDOF, NNDOF, NPNDOF, M, M_neib
    integer(kint) :: i, iter, iter_RR
    real(kdouble) :: alpha, beta, rho, rho1, omega, B2
    real(kdouble) :: tspmv, tdotp, tcomm_spmv, tcomm_dotp, tdemv
    logical :: is_converge
    integer(kint), allocatable :: IPV_R(:)
    real(kdouble), allocatable :: R(:), Z(:), Q(:), P(:), X0(:), Qb(:), PtX(:)
    real(kdouble), allocatable :: W(:,:), AW(:,:), WtA(:,:), WtW(:,:)
    real(kdouble), pointer :: B(:), X(:)

    call monolis_std_debug_log_header("monolis_solver_ADeflatedCG2")

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
    call monolis_alloc_R_1d(Qb, NPNDOF)
    call monolis_alloc_R_1d(PtX, NPNDOF)

    call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
    call monolis_set_converge_R(monoCOM, monoMAT, R, B2, is_converge, tdotp, tcomm_dotp)
    if(is_converge) return

    if(M > 0)then
      call deflatedCG_set_deflation_mode(monoPRM, monoCOM, monoMAT, NPNDOF, M, M_neib, W)
      call deflatedCG_E_initialize(monoCOM, monoPRM, monoMAT, monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, &
        & M, M_neib, NNDOF, W, AW, WtA)

      call deflatedCG_Q(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, &
        & M, NNDOF, W, B, Qb, tdemv)

      if(monoPRM%Iarray(monolis_prm_I_is_init_x) /= monolis_I_true)then
        call deflatedCG_Pt(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, &
          & M, M_neib, NNDOF, W, WtA, X, PtX, tdemv)
      endif

      call monolis_vec_AXPBY_R(N, NDOF, 1.0d0, Qb, 1.0d0, PtX, X)
      call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)

      call deflatedCG_residual_replacement_initialize(M, NNDOF, W, WtW, IPV_R)
    endif

    call monolis_inner_product_main_R(monoCOM, N, NDOF, R, R, rho, tdotp, tcomm_dotp)

    if(rho/B2 < monoPRM%Rarray(monolis_prm_R_tol))then
      monoPRM%Iarray(monolis_prm_I_cur_iter) = 0
      monoPRM%Rarray(monolis_prm_R_cur_resid) = rho/B2
      return
    endif

    call monolis_precond_apply_R(monoPRM, monoCOM, monoMAT, monoPREC, R, Z)

    call deflatedCG_Pt(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, &
      & M, M_neib, NNDOF, W, WtA, Z, PtX, tdemv)

    call deflatedCG_Q(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, &
      & M, NNDOF, W, R, Qb, tdemv)

    call monolis_vec_AXPBY_R(N, NDOF, 1.0d0, PtX, 1.0d0, Qb, Z)

    call monolis_vec_copy_R(1, NNDOF, Z, P)

    call monolis_inner_product_main_R(monoCOM, N, NDOF, R, Z, rho, tdotp, tcomm_dotp)

    do iter = 1, monoPRM%Iarray(monolis_prm_I_max_iter)
      call monolis_matvec_product_main_R(monoCOM, monoMAT, P, Q, tspmv, tcomm_spmv)
      call monolis_inner_product_main_R(monoCOM, N, NDOF, P, Q, omega, tdotp, tcomm_dotp)

      alpha = rho/omega

      call monolis_vec_AXPBY_R(N, NDOF, alpha, P, 1.0d0, X, X)
      call monolis_vec_AXPBY_R(N, NDOF,-alpha, Q, 1.0d0, R, R)

      if(mod(iter, iter_RR) == 0)then
        if(M > 0)then
          call deflatedCG_residual_replacement(M, N, NDOF, W, R, WtW, IPV_R)
        endif
      endif

      call monolis_check_converge_R(monoPRM, monoCOM, monoMAT, R, B2, iter, is_converge, tdotp, tcomm_dotp)
      if(is_converge) exit

      call monolis_precond_apply_R(monoPRM, monoCOM, monoMAT, monoPREC, R, Z)

      call deflatedCG_Pt(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, &
        & M, M_neib, NNDOF, W, WtA, Z, PtX, tdemv)

      call deflatedCG_Q(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, &
        & M, NNDOF, W, R, Qb, tdemv)

      call monolis_vec_AXPBY_R(N, NDOF, 1.0d0, PtX, 1.0d0, Qb, Z)

      rho1 = rho

      call monolis_inner_product_main_R(monoCOM, N, NDOF, R, Z, rho, tdotp, tcomm_dotp)

      beta = rho/rho1

      call monolis_vec_AXPBY_R(N, NDOF, beta, P, 1.0d0, Z, P)
    enddo

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
    call monolis_dealloc_R_1d(Qb)
    call monolis_dealloc_R_1d(PtX)

    if(M > 0)then
      call deflatedCG_finalize(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, &
        & IPV_R, W, AW, WtA, WtW)
    endif
  end subroutine monolis_solver_ADeflatedCG2

end module mod_monolis_solver_DeflatedCG
