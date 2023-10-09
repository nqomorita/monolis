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
    type(monolis_mat) :: monoPRE_deflated_eq
    integer(kint) :: N, NP, NDOF, NNDOF, NPNDOF, M, M_neib
    integer(kint) :: i, iter, iter_RR
    real(kdouble) :: alpha, beta, rho, rho1, omega, B2
    real(kdouble) :: tspmv, tdotp, tcomm_spmv, tcomm_dotp, tdemv
    logical :: is_converge
    integer(kint), allocatable :: IPV_R(:)
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
      call deflatedCG_set_deflation_mode(monoPRM, monoCOM, monoMAT, NPNDOF, M, M_neib, W)
      call deflatedCG_initialize(monoCOM, monoPRM, monoMAT, monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, &
        & M, M_neib, NNDOF, NPNDOF, IPV_R, W, WtA, WtW, g)
      call deflatedCG_pre(monoPRM, monoCOM, monoMAT, &
        & monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, &
        & M, M_neib, NNDOF, W, X, X0, g, B, R, WTW, IPV_R)
    endif

    call monolis_set_converge_R(monoCOM, monoMAT, R, B2, is_converge, tdotp, tcomm_dotp)

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
        call deflatedCG_omit(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, &
          & M, M_neib, NNDOF, W, WtA, P, Z, tdemv)
      endif

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
    enddo

    if(M > 0)then
      call deflatedCG_post(monoPRM, monoCOM, monoMAT, &
        & monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, M, M_neib, W, X, X0, g)
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
      call deflatedCG_finalize(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, &
        & IPV_R, W, WtA, WtW, g)
    endif
  end subroutine monolis_solver_DeflatedCG

  !> @ingroup dev_solver
  !> Deflated CG 法
  subroutine deflatedCG_initialize(monoCOM, monoPRM, monoMAT, monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, &
    & M, M_neib, NNDOF, NPNDOF, IPV_R, W, WtA, WtW, g)
    implicit none
    !> [in] 通信テーブル構造体
    type(monolis_com) :: monoCOM
    !> [in] パラメータ構造体
    type(monolis_prm) :: monoPRM
    !> [in] 行列構造体
    type(monolis_mat) :: monoMAT
    type(monolis_prm) :: monoPRM_deflated_eq
    type(monolis_com) :: monoCOM_deflated_eq
    type(monolis_com) :: monoCOM_deflated_eq_self
    type(monolis_mat) :: monoMAT_deflated_eq
    integer(kint) :: NNDOF, NPNDOF, M, M_neib
    integer(kint) :: i, NP
    real(kdouble) :: tdemv, time
    integer(kint), allocatable :: global_id(:)
    integer(kint), allocatable :: IPV_R(:)
    real(kdouble), allocatable :: W(:,:), WtA(:,:), WtW(:,:), g(:)
    real(kdouble), allocatable :: AW(:,:), L(:,:)

    !# allocation
    call monolis_alloc_R_2d(WtA, M_neib, NNDOF)
    call monolis_alloc_R_2d(WtW, M, M)
    call monolis_alloc_R_1d(g, M_neib)
    call monolis_alloc_I_1d(IPV_R, M)

    call monolis_alloc_R_2d(AW, NNDOF, M_neib)
    call monolis_alloc_R_2d(L, M, M_neib)

    !# initizalization
    call monolis_prm_initialize(monoPRM_deflated_eq)
    call monolis_mat_initialize(monoMAT_deflated_eq)
    call monolis_com_initialize_by_self(monoCOM_deflated_eq_self)

    monoPRM_deflated_eq%Iarray(monolis_prm_I_show_iterlog) = .false.
    monoPRM_deflated_eq%Iarray(monolis_prm_I_show_time) = .false.
    monoPRM_deflated_eq%Iarray(monolis_prm_I_show_summary) = .false.
    monoPRM_deflated_eq%Iarray(monolis_prm_I_show_time_statistics) = .false.

    NP = monoCOM%recv_n_neib + 1
    call monolis_alloc_I_1d(global_id, NP)

    global_id(1) = monolis_mpi_get_local_my_rank(monoCOM%comm)
    call monolis_mpi_update_I(monoCOM, 1, global_id, time)
    call monolis_com_initialize_by_global_id(monoCOM_deflated_eq, monoCOM%comm, 1, NP, global_id)

    monoMAT_deflated_eq%N = NP
    monoMAT_deflated_eq%NP = NP
    monoMAT_deflated_eq%NDOF = M

    call monolis_palloc_I_1d(monoMAT_deflated_eq%CSR%index, 2)
    call monolis_palloc_I_1d(monoMAT_deflated_eq%CSR%item, NP)

    monoMAT_deflated_eq%CSR%index(1) = 0
    monoMAT_deflated_eq%CSR%index(2) = NP

    do i = 1, NP
      monoMAT_deflated_eq%CSR%item(i) = i
    enddo

    call monolis_alloc_nonzero_pattern_mat_val_R(monoMAT_deflated_eq)

    !# matrix value assign
    call monolis_matmat_product_main_local_R(monoCOM_deflated_eq_self, monoMAT, M_neib, W, AW, time, time)
    call monolis_dense_matmul_local_R(M, NNDOF, M_neib, transpose(W), AW, L, tdemv)
    WtA = transpose(AW)

    call monolis_set_block_to_sparse_matrix_main_R(monoMAT_deflated_eq%CSR%index, monoMAT_deflated_eq%CSR%item, &
      & monoMAT_deflated_eq%R%A, monoMAT_deflated_eq%ndof, 1, 1, L(1:M,1:M))

    do i = 1, monoCOM%recv_n_neib  - 1
      call monolis_set_block_to_sparse_matrix_main_R(monoMAT_deflated_eq%CSR%index, monoMAT_deflated_eq%CSR%item, &
        & monoMAT_deflated_eq%R%A, monoMAT_deflated_eq%ndof, 1, i + 1, L(1:M, i*M + 1:i*M + M))
    enddo
  end subroutine deflatedCG_initialize

  subroutine deflatedCG_set_deflation_mode(monoPRM, monoCOM, monoMAT, NPNDOF, M, M_neib, W)
    implicit none
    !> monolis 構造体
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    !> 縮約次数 (local)
    integer(kint) :: NPNDOF, M, M_neib
    !> 縮約基底ベクトル
    real(kdouble), allocatable :: W(:,:)
    real(kdouble) :: tcomm
    integer(kint) :: iS, in, i, j, k, ierr, id, NDOF

    call monolis_alloc_R_2d(W, NPNDOF, M_neib)
    W = monoPRM%deflation_mode

    NDOF  = monoMAT%NDOF
    W(:,1:M) = monoPRM%deflation_mode(:,1:M)

    do i = 1, M
      call monolis_mpi_update_R(monoCOM, NDOF, W(:,i), tcomm)
    enddo

    id = M
    do i = 1, monoCOM%recv_n_neib
      iS = monoCOM%recv_index(i-1)
      in = monoCOM%recv_index(i  ) - iS
      do j = 1, M
        id = id + 1
        do k = iS + 1, iS + in
          W(NDOF*(monoCOM%recv_item(k)-1)+1:NDOF*(monoCOM%recv_item(k)),id) = &
        & W(NDOF*(monoCOM%recv_item(k)-1)+1:NDOF*(monoCOM%recv_item(k)),j)

          W(NDOF*(monoCOM%recv_item(k)-1)+1:NDOF*(monoCOM%recv_item(k)),j) = 0.0d0
        enddo
      enddo
    enddo
  end subroutine deflatedCG_set_deflation_mode

  !> @ingroup dev_solver
  !> Deflated CG 法
  subroutine deflatedCG_pre(monoPRM, monoCOM, monoMAT, &
    & monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, &
    & M, M_neib, NNDOF, W, X, X0, g, B, R, WTW, IPV_R)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    type(monolis_prm) :: monoPRM_deflated_eq
    type(monolis_com) :: monoCOM_deflated_eq
    type(monolis_mat) :: monoMAT_deflated_eq
    type(monolis_mat) :: monoPRE_deflated_eq
    !> 縮約次数 (local)
    integer(kint), intent(in) :: M
    !> 縮約次数 (local+隣接領域)
    integer(kint), intent(in) :: M_neib
    !> 縮約基底ベクトル
    real(kdouble), intent(in) :: W(:,:)
    !> X
    real(kdouble), intent(inout) :: X(:)
    !> X0
    real(kdouble), intent(inout) :: X0(:)
    !> g
    real(kdouble), intent(inout) :: g(:)
    !> 右辺ベクトル
    real(kdouble), intent(in) :: B(:)
    !> 残差ベクトル
    real(kdouble), intent(inout) :: R(:)
    !> WTW
    real(kdouble), intent(inout) :: WTW(:,:)
    !> LUピボット情報(残差再計算用)
    integer(kint), intent(inout) :: IPV_R(:)
    !> 計算時間
    real(kdouble) :: tdemv
    integer(kint) :: NNDOF
    real(kdouble), allocatable :: WTR(:), Wg(:)
    real(kdouble) :: time

    call monolis_alloc_R_1d(WTR, M_neib)
    call monolis_alloc_R_1d(Wg, NNDOF)

    monoMAT_deflated_eq%R%B = 0.0d0
    call monolis_dense_matvec_local_R(M, NNDOF, transpose(W), R, monoMAT_deflated_eq%R%B, tdemv)
    !> call monolis_lapack_dense_matvec_local_R(M_neib, NNDOF, WtA, Z, WtAZ, tdemv)

    call interface_monolis_solve_main_R(monoPRM_deflated_eq, monoCOM_deflated_eq, &
      & monoMAT_deflated_eq, monoPRE_deflated_eq)

    g = monoMAT_deflated_eq%R%X

    call monolis_dense_matvec_local_R(NNDOF, M, W, g, Wg, tdemv)
    call monolis_vec_AXPBY_R(NNDOF, 1, 1.0d0, X, 1.0d0, Wg, X)
    call monolis_vec_copy_R(NNDOF, 1, X, X0)

    call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, time, time)

    call monolis_dense_matvec_local_R(M, NNDOF, transpose(W), R, monoMAT_deflated_eq%R%B, tdemv)
    !> call monolis_lapack_dense_matvec_trans_local_R(M, NNDOF, W, R, monoMAT_deflated_eq%R%B, tdemv)

    call interface_monolis_solve_main_R(monoPRM_deflated_eq, monoCOM_deflated_eq, &
      & monoMAT_deflated_eq, monoPRE_deflated_eq)

    g = monoMAT_deflated_eq%R%X

    call monolis_dense_matmul_local_R(M, NNDOF, M, transpose(W), W, WTW, tdemv)
    call monolis_lapack_LU_fact_R(M, WTW, IPV_R)
  end subroutine deflatedCG_pre

  !> @ingroup dev_solver
  !> Deflated CG 法
  subroutine deflatedCG_omit(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, &
      & M, M_neib, NNDOF, W, WtA, P, Z, tdemv)
    implicit none
    type(monolis_prm) :: monoPRM_deflated_eq
    type(monolis_com) :: monoCOM_deflated_eq
    type(monolis_mat) :: monoMAT_deflated_eq
    type(monolis_mat) :: monoPRE_deflated_eq
    integer(kint) :: M
    integer(kint) :: M_neib
    integer(kint) :: NNDOF
    real(kdouble) :: W(:,:)
    real(kdouble) :: WtA(:,:)
    real(kdouble) :: P(:)
    real(kdouble) :: Z(:)
    real(kdouble) :: tdemv, time, WtAZ(M_neib), WLinvWtAZ(NNDOF)

    call monolis_dense_matvec_local_R(M_neib, NNDOF, WtA, Z, WtAZ, tdemv)
    !> call monolis_lapack_dense_matvec_local_R(M_neib, NNDOF, WtA, Z, WtAZ, tdemv)

    monoMAT_deflated_eq%R%B(1:M) = WtAZ(1:M)

    call monolis_mpi_update_reverse_R(monoCOM_deflated_eq, M, WtAZ, time)
    monoMAT_deflated_eq%R%B(1:M) = monoMAT_deflated_eq%R%B(1:M) + WtAZ(1:M)

    call interface_monolis_solve_main_R(monoPRM_deflated_eq, monoCOM_deflated_eq, &
      & monoMAT_deflated_eq, monoPRE_deflated_eq)

    call monolis_dense_matvec_local_R(NNDOF, M, W, monoMAT_deflated_eq%R%X, WLinvWtAZ, tdemv)
    !> call monolis_lapack_dense_matvec_local_R(NNDOF, M, W, monoMAT_deflated_eq%R%X, WLinvWtAZ, tdemv)

    call monolis_vec_AXPBY_R(NNDOF, 1, -1.0d0, WLinvWtAZ, 1.0d0, P, P)
  end subroutine deflatedCG_omit

  !> @ingroup dev_solver
  !> Deflated CG 法
  subroutine deflatedCG_residual_replacement(M, N, NDOF, W, R, WtW, IPV_R)
    implicit none
    !> 縮約次数 (local)
    integer(kint), intent(in) :: M
    !> 縮約次数 (local)
    integer(kint), intent(in) :: N
    !> 縮約次数 (local)
    integer(kint), intent(in) :: NDOF
    !> W
    real(kdouble), intent(in) :: W(:,:)
    !> 残差ベクトル
    real(kdouble), intent(inout) :: R(:)
    !> WTW
    real(kdouble), intent(inout) :: WtW(:,:)
    !> LUピボット情報(残差再計算用)
    integer(kint), intent(inout) :: IPV_R(:)
    integer(kint) :: NNDOF
    real(kdouble) :: WTR(M), time
    real(kdouble), allocatable :: WWtWinvWtR(:)

    NNDOF = N*NDOF
    call monolis_alloc_R_1d(WWtWinvWtR, NNDOF)
    call monolis_dense_matvec_local_R(M, NNDOF, transpose(W(1:NNDOF,1:M)), R(1:NNDOF), WTR, time)
    !> call monolis_lapack_dense_matvec_trans_local_R(M, NNDOF, W, R, WTR, time)

    call monolis_lapack_LU_solve_R(M, WTW, WTR, IPV_R, WTR)
    call monolis_dense_matvec_local_R(NNDOF, M, W, WTR, WWtWinvWtR, time)
    !> call monolis_lapack_dense_matvec_local_R(NNDOF, M, W, WTR, WWtWinvWtR, time)

    call monolis_vec_AXPBY_R(N, NDOF, -1.0d0, WWtWinvWtR, 1.0d0, R, R)
  end subroutine deflatedCG_residual_replacement

  !> @ingroup dev_solver
  !> Deflated CG 法
  subroutine deflatedCG_post(monoPRM, monoCOM, monoMAT, &
    & monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, M, M_neib, W, X, X0, g)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    type(monolis_prm) :: monoPRM_deflated_eq
    type(monolis_com) :: monoCOM_deflated_eq
    type(monolis_mat) :: monoMAT_deflated_eq
    type(monolis_mat) :: monoPRE_deflated_eq
    !> 縮約次数 (local)
    integer(kint) :: M
    !> 縮約次数 (local+隣接領域)
    integer(kint) :: M_neib
    !> 縮約基底ベクトル
    real(kdouble) :: W(:,:)
    !> X
    real(kdouble) :: X(:)
    !> X0
    real(kdouble) :: X0(:)
    !> g
    real(kdouble) :: g(:)
    integer(kint) :: N, NP, NDOF, NNDOF, NPNDOF
    real(kdouble) :: time
    real(kdouble), allocatable :: AX(:), h(:), WgMh(:), s(:)

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NNDOF = N*NDOF
    NPNDOF= NP*NDOF

    call monolis_alloc_R_1d(AX, NPNDOF)
    call monolis_alloc_R_1d(s, NPNDOF)
    call monolis_alloc_R_1d(h, M_neib)
    call monolis_alloc_R_1d(WgMh, NNDOF)

    S = X - X0
    call monolis_matvec_product_main_R(monoCOM, monoMAT, s, AX, time, time)
    call monolis_dense_matvec_local_R(M, NNDOF, transpose(W(1:NNDOF,1:M)), AX(1:NNDOF), monoMAT_deflated_eq%R%B(1:M), time)
    !> call monolis_lapack_dense_matvec_trans_local_R(M, NNDOF, W, AX, monoMAT_deflated_eq%B, time)

    call interface_monolis_solve_main_R(monoPRM_deflated_eq, monoCOM_deflated_eq, &
      & monoMAT_deflated_eq, monoPRE_deflated_eq)

    h = monoMAT_deflated_eq%R%X
    call monolis_dense_matvec_local_R(NNDOF, M_neib, W, g - h, WgMh, time)
    call monolis_vec_AXPBY_R(N, NDOF, 1.0d0, WgMh, 1.0d0, X, X)
    call monolis_mpi_update_R(monoCOM, NDOF, X, time)
  end subroutine deflatedCG_post

  !> @ingroup dev_solver
  !> Deflated CG 法
  subroutine deflatedCG_finalize(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, &
    & IPV_R, W, WtA, WtW, g)
    implicit none
    type(monolis_prm) :: monoPRM_deflated_eq
    type(monolis_com) :: monoCOM_deflated_eq
    type(monolis_mat) :: monoMAT_deflated_eq
    integer(kint), allocatable :: IPV_R(:)
    real(kdouble), allocatable :: W(:,:), WtA(:,:), WtW(:,:), g(:)

    call monolis_dealloc_R_2d(W)
    call monolis_dealloc_R_2d(WtA)
    call monolis_dealloc_R_2d(WtW)
    call monolis_dealloc_R_1d(g)
    call monolis_dealloc_I_1d(IPV_R)

    call monolis_prm_finalize(monoPRM_deflated_eq)
    call monolis_mat_finalize(monoMAT_deflated_eq)
    call monolis_com_finalize(monoCOM_deflated_eq)
  end subroutine deflatedCG_finalize
end module mod_monolis_solver_DeflatedCG
