!> Deflated CG 法モジュール
module mod_monolis_solver_DeflatedCG_util
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
  use mod_monolis_spmat_handler
  use mod_monolis_spmat_nonzero_pattern_util

  implicit none

contains

  !> @ingroup dev_solver
  !> Deflated CG 法
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

    NDOF  = monoMAT%NDOF
    W(:,1:M) = monoPRM%deflation_mode(:,1:M)

    do i = 1, M
      call monolis_mpi_update_R(monoCOM, NDOF, W(:,i), tcomm)
    enddo

    id = M
    do i = 1, monoCOM%recv_n_neib
      iS = monoCOM%recv_index(i)
      in = monoCOM%recv_index(i + 1) - iS
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
  !> @detail E = W^t A W
  subroutine deflatedCG_E_initialize(monoCOM, monoPRM, monoMAT, monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, &
    & M, M_neib, NNDOF, W, AW, WtA)
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
    integer(kint), intent(in) :: M
    integer(kint), intent(in) :: M_neib
    integer(kint), intent(in) :: NNDOF
    integer(kint) :: i, NP, n_dof_loc
    real(kdouble) :: tdemv, time
    real(kdouble) :: W(:,:)
    real(kdouble), allocatable :: AW(:,:), WtA(:,:), L(:,:)
    logical :: is_coarse_E = .false.

    !# allocation
    call monolis_alloc_R_2d(AW, NNDOF, M_neib)
    call monolis_alloc_R_2d(WtA, M_neib, NNDOF)
    call monolis_alloc_R_2d(L, M, M_neib)

    !# initizalization
    call monolis_prm_initialize(monoPRM_deflated_eq)
    call monolis_mat_initialize(monoMAT_deflated_eq)
    call monolis_com_initialize_by_self(monoCOM_deflated_eq_self)

    monoPRM_deflated_eq%Iarray(monolis_prm_I_method) = monolis_iter_CG
    monoPRM_deflated_eq%Iarray(monolis_prm_I_precond) = monolis_prec_Diag
    !monoPRM_deflated_eq%Iarray(monolis_prm_I_precond) = monolis_prec_MUMPS
    monoPRM_deflated_eq%Rarray(monolis_prm_R_tol) = 1.0d-10
    monoPRM_deflated_eq%Iarray(monolis_prm_I_show_iterlog) = .false.
    monoPRM_deflated_eq%Iarray(monolis_prm_I_show_time) = .false.
    monoPRM_deflated_eq%Iarray(monolis_prm_I_show_summary) = .false.
    monoPRM_deflated_eq%Iarray(monolis_prm_I_show_time_statistics) = .false.
    monoPRM_deflated_eq%Iarray(monolis_prm_I_max_iter) = 10000
    monoPRM_deflated_eq%Iarray(monolis_prm_I_is_prec_stored) = 1

    !> com section
    NP = monoCOM%recv_n_neib + 1

    monoCOM_deflated_eq%my_rank = monoCOM%my_rank
    monoCOM_deflated_eq%comm = monoCOM%comm
    monoCOM_deflated_eq%comm_size = monoCOM%comm_size
    monoCOM_deflated_eq%n_internal_vertex = 1

    if(monoCOM%recv_n_neib == 0)then
      monoCOM_deflated_eq%recv_n_neib = 0
      monoCOM_deflated_eq%send_n_neib = 0

      call monolis_palloc_I_1d(monoCOM_deflated_eq%recv_neib_pe, 1)
      call monolis_palloc_I_1d(monoCOM_deflated_eq%recv_index, 2)
      call monolis_palloc_I_1d(monoCOM_deflated_eq%recv_item, 1)
      call monolis_palloc_I_1d(monoCOM_deflated_eq%send_neib_pe, 1)
      call monolis_palloc_I_1d(monoCOM_deflated_eq%send_index, 2)
      call monolis_palloc_I_1d(monoCOM_deflated_eq%send_item, 1)
    else
      monoCOM_deflated_eq%recv_n_neib = monoCOM%recv_n_neib
      call monolis_palloc_I_1d(monoCOM_deflated_eq%recv_neib_pe, monoCOM%recv_n_neib)
      call monolis_palloc_I_1d(monoCOM_deflated_eq%recv_index, monoCOM%recv_n_neib + 1)
      call monolis_palloc_I_1d(monoCOM_deflated_eq%recv_item, monoCOM%recv_n_neib)

      monoCOM_deflated_eq%recv_neib_pe = monoCOM%recv_neib_pe

      do i = 1, monoCOM_deflated_eq%recv_n_neib
        monoCOM_deflated_eq%recv_index(i + 1) = i
      enddo

      do i = 1, monoCOM_deflated_eq%recv_n_neib
        monoCOM_deflated_eq%recv_item(i) = 1 + i
      enddo

      monoCOM_deflated_eq%send_n_neib = monoCOM%send_n_neib
      call monolis_palloc_I_1d(monoCOM_deflated_eq%send_neib_pe, monoCOM%send_n_neib)
      call monolis_palloc_I_1d(monoCOM_deflated_eq%send_index, monoCOM%send_n_neib + 1)
      call monolis_palloc_I_1d(monoCOM_deflated_eq%send_item, monoCOM%recv_n_neib)

      monoCOM_deflated_eq%send_neib_pe = monoCOM%send_neib_pe

      do i = 1, monoCOM_deflated_eq%send_n_neib
        monoCOM_deflated_eq%send_index(i + 1) = i
      enddo

      do i = 1, monoCOM_deflated_eq%send_n_neib
        monoCOM_deflated_eq%send_item(i) = 1
      enddo
    endif

    !> mat section
    call monolis_matmat_product_main_local_R(monoCOM_deflated_eq_self, monoMAT, M_neib, W, AW, time, time)
    call monolis_dense_matmul_local_R(M, NNDOF, M_neib, transpose(W), AW, L, tdemv)
    WtA = transpose(AW)

if(is_coarse_E)then
    !> sparse block
    n_dof_loc = monoPRM%Iarray(monolis_prm_I_n_local_block_size_of_AtAW)
    if(n_dof_loc == 0) stop "monolis_get_sparse_matrix_from_dense_matrix_R"
    call monolis_get_sparse_matrix_from_dense_matrix_R(monoMAT_deflated_eq, M/n_dof_loc, M_neib/n_dof_loc, n_dof_loc, L)
else
    !> dense block
    monoMAT_deflated_eq%N = 1
    monoMAT_deflated_eq%NP = NP
    monoMAT_deflated_eq%NDOF = M

    call monolis_palloc_I_1d(monoMAT_deflated_eq%CSR%index, NP + 1)
    call monolis_palloc_I_1d(monoMAT_deflated_eq%CSR%item, NP)

    monoMAT_deflated_eq%CSR%index(1) = 0

    do i = 1, NP
      monoMAT_deflated_eq%CSR%index(1 + i) = NP
      monoMAT_deflated_eq%CSR%item(i) = i
    enddo

    call monolis_alloc_nonzero_pattern_mat_val_R(monoMAT_deflated_eq)

    !# matrix value assign
    call monolis_set_block_to_sparse_matrix_main_R(monoMAT_deflated_eq%CSR%index, monoMAT_deflated_eq%CSR%item, &
      monoMAT_deflated_eq%R%A, monoMAT_deflated_eq%n_dof_index, monoMAT_deflated_eq%n_dof_index2, &
      1, 1, L(1:M,1:M))

    do i = 1, monoCOM%recv_n_neib
      call monolis_set_block_to_sparse_matrix_main_R(monoMAT_deflated_eq%CSR%index, monoMAT_deflated_eq%CSR%item, &
        monoMAT_deflated_eq%R%A, monoMAT_deflated_eq%n_dof_index, monoMAT_deflated_eq%n_dof_index2, &
        1, i + 1, L(1:M, i*M + 1:i*M + M))
    enddo
endif
  end subroutine deflatedCG_E_initialize

  !> @ingroup dev_solver
  !> Deflated CG 法
  subroutine deflatedCG_E(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, monoPRE_deflated_eq, &
    & M, X, B)
    implicit none
    type(monolis_prm) :: monoPRM_deflated_eq
    type(monolis_com) :: monoCOM_deflated_eq
    type(monolis_mat) :: monoMAT_deflated_eq
    type(monolis_mat) :: monoPRE_deflated_eq
    integer(kint) :: M
    real(kdouble) :: X(:)
    real(kdouble) :: B(:)

    monoMAT_deflated_eq%R%B(1:M) = B(1:M)

    call interface_monolis_solve_main_R(monoPRM_deflated_eq, monoCOM_deflated_eq, &
      & monoMAT_deflated_eq, monoPRE_deflated_eq)

    X(1:M) = monoMAT_deflated_eq%R%X(1:M)
  end subroutine deflatedCG_E

  !> @ingroup dev_solver
  !> Deflated CG 法
  subroutine deflatedCG_Q(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, monoPRE_deflated_eq, &
      & M, NNDOF, W, X, Y, tdemv)
    implicit none
    type(monolis_prm) :: monoPRM_deflated_eq
    type(monolis_com) :: monoCOM_deflated_eq
    type(monolis_mat) :: monoMAT_deflated_eq
    type(monolis_mat) :: monoPRE_deflated_eq
    integer(kint) :: M
    integer(kint) :: NNDOF
    real(kdouble) :: W(:,:)
    real(kdouble) :: Y(:)
    real(kdouble) :: X(:)
    real(kdouble) :: tdemv, time, WtZ(M)

    if(M == 0)then
      Y = 0.0d0
      return
    endif

    call monolis_dense_matvec_local_R(M, NNDOF, transpose(W), X, WtZ, tdemv)

    monoMAT_deflated_eq%R%B(1:M) = WtZ(1:M)

    call deflatedCG_E(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, monoPRE_deflated_eq, &
    & M, monoMAT_deflated_eq%R%X, monoMAT_deflated_eq%R%B)

    call monolis_dense_matvec_local_R(NNDOF, M, W, monoMAT_deflated_eq%R%X, Y, tdemv)
 end subroutine deflatedCG_Q

  !> @ingroup dev_solver
  !> Deflated CG 法
  subroutine deflatedCG_P(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, monoPRE_deflated_eq, &
      & M, M_neib, NNDOF, W, AW, Z, P, tdemv)
    implicit none
    type(monolis_prm) :: monoPRM_deflated_eq
    type(monolis_com) :: monoCOM_deflated_eq
    type(monolis_mat) :: monoMAT_deflated_eq
    type(monolis_mat) :: monoPRE_deflated_eq
    integer(kint) :: M
    integer(kint) :: M_neib
    integer(kint) :: NNDOF
    real(kdouble) :: W(:,:)
    real(kdouble) :: AW(:,:)
    real(kdouble) :: P(:)
    real(kdouble) :: Z(:)
    real(kdouble) :: tdemv, time, WtZ(M_neib), AWEinvWtZ(NNDOF)

    if(M == 0)then
      call monolis_vec_copy_R(NNDOF, Z, P)
      return
    endif

    call monolis_dense_matvec_local_R(M, NNDOF, transpose(W), Z, WtZ, tdemv)

    monoMAT_deflated_eq%R%B(1:M) = WtZ(1:M)

    call deflatedCG_E(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, monoPRE_deflated_eq, &
    & M, monoMAT_deflated_eq%R%X, monoMAT_deflated_eq%R%B)

    call monolis_dense_matvec_local_R(NNDOF, M_neib, AW, monoMAT_deflated_eq%R%X, AWEinvWtZ, tdemv)

    call monolis_vec_AXPBY_R(NNDOF, -1.0d0, AWEinvWtZ, 1.0d0, Z, P)
  end subroutine deflatedCG_P

  !> @ingroup dev_solver
  !> Deflated CG 法
  subroutine deflatedCG_P_coarse(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, monoPRE_deflated_eq, &
      & M, M_neib, NNDOF, monoMAT_Wt, monoMAT_AW, monoCOM_self, Z, P, tdemv)
    implicit none
    type(monolis_prm) :: monoPRM_deflated_eq
    type(monolis_com) :: monoCOM_deflated_eq
    type(monolis_mat) :: monoMAT_deflated_eq
    type(monolis_mat) :: monoPRE_deflated_eq
    type(monolis_mat) :: monoMAT_Wt
    type(monolis_mat) :: monoMAT_AW
    type(monolis_com) :: monoCOM_self
    integer(kint) :: M
    integer(kint) :: M_neib
    integer(kint) :: NNDOF
    real(kdouble) :: P(:)
    real(kdouble) :: Z(:)
    real(kdouble) :: tdemv, time, WtZ(M_neib), AWEinvWtZ(NNDOF)

    if(M == 0)then
      call monolis_vec_copy_R(NNDOF, Z, P)
      return
    endif

    !call monolis_dense_matvec_local_R(M, NNDOF, transpose(W), Z, WtZ, tdemv)
    call monolis_matvec_product_main_R(monoCOM_self, monoMAT_Wt, Z, WtZ, time, time)

    monoMAT_deflated_eq%R%B(1:M) = WtZ(1:M)

    call deflatedCG_E(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, monoPRE_deflated_eq, &
    & M, monoMAT_deflated_eq%R%X, monoMAT_deflated_eq%R%B)

    !call monolis_dense_matvec_local_R(NNDOF, M_neib, AW, monoMAT_deflated_eq%R%X, AWEinvWtZ, tdemv)
    call monolis_matvec_product_main_R(monoCOM_self, monoMAT_AW, monoMAT_deflated_eq%R%X, AWEinvWtZ, time, time)

    call monolis_vec_AXPBY_R(NNDOF, -1.0d0, AWEinvWtZ, 1.0d0, Z, P)
  end subroutine deflatedCG_P_coarse

  subroutine deflatedCG_get_coarse_W(NNDOF, M, W, monoMAT_Wt)
    implicit none
    type(monolis_mat) :: monoMAT_Wt
    integer(kint) :: M
    integer(kint) :: NNDOF
    real(kdouble) :: W(:,:)
    integer(kint) :: i, j, in

    monoMAT_Wt%N = M
    monoMAT_Wt%NP = M
    monoMAT_Wt%NDOF = 1

    call monolis_palloc_I_1d(monoMAT_Wt%CSR%index, M + 1)

    in = 0
    do i = 1, M
      do j = 1, NNDOF
        if(W(j,i) /= 0.0d0)then
          in = in + 1
          monoMAT_Wt%CSR%index(i + 1) = monoMAT_Wt%CSR%index(i + 1) + 1
        endif
      enddo
    enddo

    do i = 1, M
      monoMAT_Wt%CSR%index(i + 1) = monoMAT_Wt%CSR%index(i + 1) + monoMAT_Wt%CSR%index(i)
    enddo

    call monolis_palloc_I_1d(monoMAT_Wt%CSR%item, in)
    call monolis_palloc_R_1d(monoMAT_Wt%R%A, in)

    in = 0
    do i = 1, M
      do j = 1, NNDOF
        if(W(j,i) /= 0.0d0)then
          in = in + 1
          monoMAT_Wt%CSR%item(in) = j
          monoMAT_Wt%R%A(in) = W(j,i)
        endif
      enddo
    enddo
  end subroutine deflatedCG_get_coarse_W

  subroutine deflatedCG_get_coarse_AW(NNDOF, M, W, monoMAT_Wt)
    implicit none
    type(monolis_mat) :: monoMAT_Wt
    integer(kint) :: M
    integer(kint) :: NNDOF
    real(kdouble) :: W(:,:)
    integer(kint) :: i, j, in

    monoMAT_Wt%N = NNDOF
    monoMAT_Wt%NP = NNDOF
    monoMAT_Wt%NDOF = 1

    call monolis_palloc_I_1d(monoMAT_Wt%CSR%index, NNDOF + 1)

    in = 0
    do i = 1, NNDOF
      do j = 1, M
        if(W(i,j) /= 0.0d0)then
          in = in + 1
          monoMAT_Wt%CSR%index(i + 1) = monoMAT_Wt%CSR%index(i + 1) + 1
        endif
      enddo
    enddo

    do i = 1, NNDOF
      monoMAT_Wt%CSR%index(i + 1) = monoMAT_Wt%CSR%index(i + 1) + monoMAT_Wt%CSR%index(i)
    enddo

    call monolis_palloc_I_1d(monoMAT_Wt%CSR%item, in)
    call monolis_palloc_R_1d(monoMAT_Wt%R%A, in)

    in = 0
    do i = 1, NNDOF
      do j = 1, M
        if(W(i,j) /= 0.0d0)then
          in = in + 1
          monoMAT_Wt%CSR%item(in) = j
          monoMAT_Wt%R%A(in) = W(i,j)
        endif
      enddo
    enddo
  end subroutine deflatedCG_get_coarse_AW

  !> @ingroup dev_solver
  !> Deflated CG 法
  subroutine deflatedCG_Pt(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, monoPRE_deflated_eq, &
      & M, M_neib, NNDOF, W, WtA, Z, P,tdemv)
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
    real(kdouble) :: tdemv, time, WtAZ(M_neib), WEinvWtAZ(NNDOF)

    if(M == 0)then
      call monolis_vec_copy_R(NNDOF, Z, P)
      return
    endif

    call monolis_dense_matvec_local_R(M_neib, NNDOF, WtA, Z, WtAZ, tdemv)

    monoMAT_deflated_eq%R%B(1:M) = WtAZ(1:M)

    WtAZ(1:M) = 0.0d0

    call monolis_mpi_update_reverse_R(monoCOM_deflated_eq, M, WtAZ, time)

    monoMAT_deflated_eq%R%B(1:M) = monoMAT_deflated_eq%R%B(1:M) + WtAZ(1:M)

    call deflatedCG_E(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, monoPRE_deflated_eq, &
    & M, monoMAT_deflated_eq%R%X, monoMAT_deflated_eq%R%B)

    call monolis_dense_matvec_local_R(NNDOF, M, W, monoMAT_deflated_eq%R%X, WEinvWtAZ, tdemv)

    call monolis_vec_AXPBY_R(NNDOF, -1.0d0, WEinvWtAZ, 1.0d0, Z, P)
  end subroutine deflatedCG_Pt

  !> @ingroup dev_solver
  !> Deflated CG 法
  subroutine deflatedCG_residual_replacement_initialize(M, NNDOF, W, WtW, IPV_R)
    implicit none
    !> 縮約次数 (local)
    integer(kint), intent(in) :: M
    integer(kint), intent(in) :: NNDOF
    !> 縮約基底ベクトル
    real(kdouble), intent(in) :: W(:,:)
    !> WtW
    real(kdouble), allocatable :: WtW(:,:)
    !> LUピボット情報(残差再計算用)
    integer(kint), allocatable :: IPV_R(:)
    !> 計算時間
    real(kdouble) :: tdemv

    call monolis_alloc_R_2d(WtW, M, M)
    call monolis_alloc_I_1d(IPV_R, M)

    call monolis_dense_matmul_local_R(M, NNDOF, M, transpose(W), W, WtW, tdemv)
    call monolis_lapack_LU_fact_R(M, WtW, IPV_R)
  end subroutine deflatedCG_residual_replacement_initialize

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
    !> WtW
    real(kdouble), intent(inout) :: WtW(:,:)
    !> LUピボット情報(残差再計算用)
    integer(kint), intent(inout) :: IPV_R(:)
    integer(kint) :: NNDOF
    real(kdouble) :: WTR(M), time
    real(kdouble), allocatable :: WWtWinvWtR(:)

    NNDOF = N*NDOF
    call monolis_alloc_R_1d(WWtWinvWtR, NNDOF)
    call monolis_dense_matvec_local_R(M, NNDOF, transpose(W(1:NNDOF,1:M)), R(1:NNDOF), WTR, time)

    call monolis_lapack_LU_solve_R(M, WtW, WTR, IPV_R, WTR)

    call monolis_dense_matvec_local_R(NNDOF, M, W, WTR, WWtWinvWtR, time)

    call monolis_vec_AXPBY_R(N*NDOF, -1.0d0, WWtWinvWtR, 1.0d0, R, R)
  end subroutine deflatedCG_residual_replacement

  !> @ingroup dev_solver
  !> Deflated CG 法
  subroutine deflatedCG_finalize(monoPRM_deflated_eq, monoCOM_deflated_eq, monoMAT_deflated_eq, monoPRE_deflated_eq, &
    & IPV_R, W, AW, WtA, WtW)
    implicit none
    type(monolis_prm) :: monoPRM_deflated_eq
    type(monolis_com) :: monoCOM_deflated_eq
    type(monolis_mat) :: monoMAT_deflated_eq
    type(monolis_mat) :: monoPRE_deflated_eq
    integer(kint), allocatable :: IPV_R(:)
    real(kdouble), allocatable :: W(:,:), AW(:,:), WtA(:,:), WtW(:,:)

    call monolis_dealloc_R_2d(W)
    call monolis_dealloc_R_2d(AW)
    call monolis_dealloc_R_2d(WtA)
    call monolis_dealloc_R_2d(WtW)
    call monolis_dealloc_I_1d(IPV_R)

    call monolis_prm_finalize(monoPRM_deflated_eq)
    call monolis_mat_finalize(monoMAT_deflated_eq)
    call monolis_mat_finalize(monoPRE_deflated_eq)
    call monolis_com_finalize(monoCOM_deflated_eq)
  end subroutine deflatedCG_finalize
end module mod_monolis_solver_DeflatedCG_util
