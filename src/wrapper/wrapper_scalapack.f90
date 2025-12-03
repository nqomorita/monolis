!> ScaLAPACK ラッパーモジュール
module mod_monolis_scalapack
  use mod_monolis_utils
  use mod_monolis_vec_util

  implicit none

contains

  !> @ingroup wrapper
  !> Scalapack コミュニケータの初期化
  subroutine monolis_scalapack_comm_initialize(comm, scalapack_comm)
    implicit none
    !> [in] コミュニケータ
    integer(kint), intent(in) :: comm
    !> [in] コミュニケータ
    integer(kint) :: scalapack_comm
    integer(kint), allocatable :: user_map(:,:)
    integer(kint) :: n_col, n_row

    !# scalapack コミュニケータの初期化処理
    n_row = monolis_mpi_get_local_comm_size(comm)
    n_col = 1

    call blacs_get(0, 0, scalapack_comm)

    call monolis_alloc_I_2d(user_map, n_row, 1)
    user_map(monolis_mpi_get_local_my_rank(comm) + 1,1) = monolis_mpi_get_global_my_rank()
    call monolis_allreduce_I(n_row, user_map(:,1), monolis_mpi_max, comm)

    call blacs_gridmap(scalapack_comm, user_map, n_row, n_row, n_col)
  end subroutine monolis_scalapack_comm_initialize

  !> @ingroup wrapper
  !> Scalapack コミュニケータの終了
  subroutine monolis_scalapack_comm_finalize(scalapack_comm)
    implicit none
    !> [in] コミュニケータ
    integer(kint), intent(in) :: scalapack_comm

    !# scalapack コミュニケータの終了処理
    call blacs_gridexit(scalapack_comm)
  end subroutine monolis_scalapack_comm_finalize

  !> @ingroup wrapper
  !> PDGESVD 関数（実数型）
  subroutine monolis_scalapack_gesvd_R(N_loc, M, A, S, V, D, comm, scalapack_comm)
    implicit none
    !> [in] 行列の大きさ（行数 N）
    integer(kint), intent(in) :: N_loc
    !> [in] 行列の大きさ（列数 M）
    integer(kint), intent(in) :: M
    !> [in] 入力行列（N_loc x M）
    real(kdouble), intent(in) :: A(:,:)
    !> [out] 左特異行列（N_loc x P）
    real(kdouble), intent(out) :: S(:,:)
    !> [out] 特異値（P）
    real(kdouble), intent(out) :: V(:)
    !> [out] 右特異行列（P x M）
    real(kdouble), intent(out) :: D(:,:)
    !> [in] コミュニケータ
    integer(kint), intent(in) :: comm
    !> [in] コミュニケータ
    integer(kint), intent(in) :: scalapack_comm
    integer(kint) :: N_loc_max, M_fix, N
    integer(kint) :: comm_size, P
    integer(kint) :: i, j
    real(kdouble), allocatable :: A_temp(:,:)
    real(kdouble), allocatable :: S_temp(:,:)
    real(kdouble), allocatable :: V_temp(:)
    real(kdouble), allocatable :: D_temp(:,:)

    !# N の取得（各領域の行数を統一）
    comm_size = monolis_mpi_get_local_comm_size(comm)
    N_loc_max = N_loc
    call monolis_allreduce_I1(N_loc_max, monolis_mpi_max, comm)

    !# M の取得（各領域の列数を統一）
    if(mod(M, comm_size) == 0)then
      M_fix = M
    else
      M_fix = (M/comm_size + 1)*comm_size
    endif

    P = min(N_loc_max*comm_size, M_fix)

    !# 係数行列のパディング
    call monolis_alloc_R_2d(A_temp, N_loc_max, M_fix)
    call monolis_alloc_R_2d(S_temp, N_loc_max, P)
    call monolis_alloc_R_1d(V_temp, P)
    call monolis_alloc_R_2d(D_temp, P, M_fix)

    do i = 1, M
    do j = 1, N_loc
      A_temp(j,i) = A(j,i)
    enddo
    enddo

    !# 特異値分解
    call monolis_scalapack_gesvd_R_main(N_loc_max, M_fix, A_temp, S_temp, V_temp, D_temp, comm, scalapack_comm)

    !# 出力行列サイズの修正
    N = N_loc
    call monolis_allreduce_I1(N, monolis_mpi_sum, comm)
    P = min(N, M)

    do i = 1, P
    do j = 1, N_loc
      S(j,i) = S_temp(j,i)
    enddo
    enddo

    do i = 1, P
      V(i) = V_temp(i)
    enddo

    do i = 1, M
    do j = 1, P
      D(j,i) = D_temp(j,i)
    enddo
    enddo
  end subroutine monolis_scalapack_gesvd_R

  !> @ingroup wrapper
  !> PDGESVD 関数（実数型、メイン関数）
  subroutine monolis_scalapack_gesvd_R_main(N_loc, M, A, S, V, D, comm, scalapack_comm)
    implicit none
    !> [in] 行列の大きさ（行数 N）
    integer(kint), intent(in) :: N_loc
    !> [in] 行列の大きさ（列数 M）
    integer(kint), intent(in) :: M
    !> [in] 入力行列（N_loc x M）
    real(kdouble), intent(in) :: A(:,:)
    !> [out] 左特異行列（N_loc x P）
    real(kdouble), intent(out) :: S(:,:)
    !> [out] 特異値（P）
    real(kdouble), intent(out) :: V(:)
    !> [out] 右特異行列（P x M）
    real(kdouble), intent(out) :: D(:,:)
    !> [in] コミュニケータ
    integer(kint), intent(in) :: comm
    !> [in] コミュニケータ
    integer(kint), intent(in) :: scalapack_comm
    integer(kint) :: N
    integer(kint) :: NB, P, desc_A(9), desc_S(9), desc_D(9)
    integer(kint) :: lld_A, lld_S, lld_D
    integer(kint) :: NW, info
    integer(kint) :: my_col, my_row, n_col, n_row
    real(kdouble), allocatable :: W(:)
    real(kdouble), allocatable :: A_temp(:,:)

    integer :: numroc
    external :: numroc

    call blacs_gridinfo(scalapack_comm, n_row, n_col, my_row, my_col)

    !# Scalapack 用パラメータの取得
    desc_A = 0
    desc_S = 0
    desc_D = 0

    NB = 1

    !# N の取得
    N = N_loc
    call monolis_allreduce_I1(N, monolis_mpi_sum, comm)
    P = min(N, M)

    lld_A = numroc(N, NB, my_row, 0, n_row)
    lld_S = numroc(N, NB, my_row, 0, n_row)
    lld_D = numroc(P, NB, my_row, 0, n_row)

    call descinit(desc_A, N, M, NB, NB, 0, 0, scalapack_comm, lld_A, info)
    call descinit(desc_S, N, P, NB, NB, 0, 0, scalapack_comm, lld_S, info)
    call descinit(desc_D, P, M, NB, NB, 0, 0, scalapack_comm, lld_D, info)

    !# 一時ベクトルの大きさ取得
    call monolis_alloc_R_1d(W, 1)
    call monolis_alloc_R_2d(A_temp, N_loc, M)
    A_temp = A

    call pdgesvd("V", "V", N, M, &
      & A_temp, 1, 1, desc_A, &
      & V, &
      & S, 1, 1, desc_S, &
      & D, 1, 1, desc_D, &
      & W, -1, info)

    NW = int(W(1))
    call monolis_dealloc_R_1d(W)
    call monolis_alloc_R_1d(W, NW)

    !# 特異値分解
    call pdgesvd("V", "V", N, M, &
      & A_temp, 1, 1, desc_A, &
      & V, &
      & S, 1, 1, desc_S, &
      & D, 1, 1, desc_D, &
      & W, NW, info)

    !# 計算結果 D 行列の通信
    call gesvd_R_update_D(n_row, P, M, lld_D, D, comm)
  end subroutine monolis_scalapack_gesvd_R_main

  !> @ingroup wrapper
  !> 右特異ベクトルのアップデート（実数型、行列情報の更新）
  subroutine gesvd_R_update_D(n_row, P, M, lld_D, D, comm)
    implicit none
    !> [in] 
    integer(kint), intent(in) :: n_row
    !> [in] 
    integer(kint), intent(in) :: P
    !> [in] 
    integer(kint), intent(in) :: M
    !> [in] 
    integer(kint), intent(in) :: lld_D
    !> [in,out] 右特異ベクトル
    real(kdouble), intent(inout) :: D(:,:)
    !> [in] コミュニケータ
    integer(kint), intent(in) :: comm
    integer(kint) :: i, j, size
    integer(kint), allocatable :: counts(:), displs(:)
    real(kdouble), allocatable :: D_temp(:)
    real(kdouble), allocatable :: D_full(:)
    real(kdouble), allocatable :: D_perm(:)

    size = lld_D*M

    call monolis_alloc_R_1d(D_temp, size)
    call monolis_alloc_R_1d(D_full, P*M)
    call monolis_alloc_R_1d(D_perm, P*M)
    call monolis_alloc_I_1d(counts, n_row)
    call monolis_alloc_I_1d(displs, n_row)

    counts = size

    do i = 2, n_row
      displs(i) = size*(i-1)
    enddo

    call monolis_mat_to_vec_R(lld_D, M, D, D_temp)
    call monolis_allgather_V_R(size, D_temp, D_full, counts, displs, comm)

    do i = 1, size
    do j = 1, n_row
      D_perm(n_row*(i-1) + j) = D_full(i + size*(j - 1))
    enddo
    enddo

    call monolis_vec_to_mat_R(P, M, D_perm, D)
  end subroutine gesvd_R_update_D

  !> @ingroup wrapper
  !> PDGETRF 関数（実数型、LU分解）
  subroutine monolis_scalapack_getrf_R(N_loc, N, A, ipiv, scalapack_comm)
    implicit none
    !> [in] 行列の大きさ（ローカル行数）
    integer(kint), intent(in) :: N_loc
    !> [in] 行列の大きさ（全体のサイズ N x N）
    integer(kint), intent(in) :: N
    !> [in,out] 入力行列（N_loc x N）、出力はLU分解後の行列
    real(kdouble), intent(inout) :: A(:,:)
    !> [out] ピボット情報（N_loc）
    integer(kint), intent(out) :: ipiv(:)
    !> [in] Scalapack コミュニケータ
    integer(kint), intent(in) :: scalapack_comm
    integer(kint) :: NB, desc_A(9)
    integer(kint) :: lld_A
    integer(kint) :: info
    integer(kint) :: my_col, my_row, n_col, n_row

    integer :: numroc
    external :: numroc

    call blacs_gridinfo(scalapack_comm, n_row, n_col, my_row, my_col)

    desc_A = 0
    NB = 1
    lld_A = numroc(N, NB, my_row, 0, n_row)

    call descinit(desc_A, N, N, NB, NB, 0, 0, scalapack_comm, lld_A, info)

    call pdgetrf(N, N, A, 1, 1, desc_A, ipiv, info)
  end subroutine monolis_scalapack_getrf_R

  !> @ingroup wrapper
  !> PDGETRS 関数（実数型、LU分解による線形方程式の求解）
  subroutine monolis_scalapack_getrs_R(N_loc, N, NRHS, A, ipiv, B, scalapack_comm)
    implicit none
    !> [in] 行列の大きさ（ローカル行数）
    integer(kint), intent(in) :: N_loc
    !> [in] 行列の大きさ（全体のサイズ N x N）
    integer(kint), intent(in) :: N
    !> [in] 右辺ベクトルの数
    integer(kint), intent(in) :: NRHS
    !> [in] LU分解された行列（N_loc x N）
    real(kdouble), intent(in) :: A(:,:)
    !> [in] ピボット情報（N_loc）
    integer(kint), intent(in) :: ipiv(:)
    !> [in,out] 右辺ベクトル（N_loc x NRHS）、出力は解ベクトル
    real(kdouble), intent(inout) :: B(:,:)
    !> [in] Scalapack コミュニケータ
    integer(kint), intent(in) :: scalapack_comm
    integer(kint) :: NB, desc_A(9), desc_B(9)
    integer(kint) :: lld_A, lld_B
    integer(kint) :: info
    integer(kint) :: my_col, my_row, n_col, n_row

    integer :: numroc
    external :: numroc

    call blacs_gridinfo(scalapack_comm, n_row, n_col, my_row, my_col)

    desc_A = 0
    desc_B = 0
    NB = 1
    lld_A = numroc(N, NB, my_row, 0, n_row)
    lld_B = numroc(N, NB, my_row, 0, n_row)

    call descinit(desc_A, N, N, NB, NB, 0, 0, scalapack_comm, lld_A, info)
    call descinit(desc_B, N, NRHS, NB, NB, 0, 0, scalapack_comm, lld_B, info)

    call pdgetrs("N", N, NRHS, A, 1, 1, desc_A, ipiv, B, 1, 1, desc_B, info)
  end subroutine monolis_scalapack_getrs_R
end module mod_monolis_scalapack
