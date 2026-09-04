!> ScaLAPACK ラッパーモジュール
module mod_monolis_scalapack
  use mod_monolis_utils
  use mod_monolis_vec_util
  use mod_monolis_lapack

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

!$omp parallel default(none) &
!$omp & shared(A, A_temp) &
!$omp & firstprivate(M, N_loc) &
!$omp & private(i, j)
!$omp do collapse(2)
    do i = 1, M
    do j = 1, N_loc
      A_temp(j,i) = A(j,i)
    enddo
    enddo
!$omp end do
!$omp end parallel

    !# 特異値分解
    call monolis_scalapack_gesvd_R_main(N_loc_max, M_fix, A_temp, S_temp, V_temp, D_temp, comm, scalapack_comm)

    !# 出力行列サイズの修正
    N = N_loc
    call monolis_allreduce_I1(N, monolis_mpi_sum, comm)
    P = min(N, M)

!$omp parallel default(none) &
!$omp & shared(S, S_temp) &
!$omp & firstprivate(P, N_loc) &
!$omp & private(i, j)
!$omp do collapse(2)
    do i = 1, P
    do j = 1, N_loc
      S(j,i) = S_temp(j,i)
    enddo
    enddo
!$omp end do
!$omp end parallel

    do i = 1, P
      V(i) = V_temp(i)
    enddo

!$omp parallel default(none) &
!$omp & shared(D, D_temp) &
!$omp & firstprivate(M, P) &
!$omp & private(i, j)
!$omp do collapse(2)
    do i = 1, M
    do j = 1, P
      D(j,i) = D_temp(j,i)
    enddo
    enddo
!$omp end do
!$omp end parallel
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

!$omp parallel default(none) &
!$omp & shared(D_perm, D_full) &
!$omp & firstprivate(size, n_row) &
!$omp & private(i, j)
!$omp do collapse(2)
    do i = 1, size
    do j = 1, n_row
      D_perm(n_row*(i-1) + j) = D_full(i + size*(j - 1))
    enddo
    enddo
!$omp end do
!$omp end parallel

    call monolis_vec_to_mat_R(P, M, D_perm, D)
  end subroutine gesvd_R_update_D

  !> @ingroup wrapper
  !> 正規方程式 (A^T A) w = A^T b の分散求解（実数型、PDGESV 使用）
  !> @details 行列 A は列ブロック分散とする（各ランクが m x p_loc の列ブロックを保持し、
  !> 全体の列順序はランク番号順に連結した順序とする）。右辺ベクトル b は全ランクに複製する。
  !> Gram 行列 A^T A をブロックサイクリック分散（ブロック幅 1）で組み立て、PDGESV で求解する。
  !> Gram 行列と求解は分散実行されるため、全体行列を単一ランクに集約しない。
  subroutine monolis_scalapack_gram_gesv_R(m, p_loc, A_loc, b, w_loc, comm, scalapack_comm)
    implicit none
    !> [in] 行列の行数（全ランク共通）
    integer(kint), intent(in) :: m
    !> [in] 自ランクが保持する列数
    integer(kint), intent(in) :: p_loc
    !> [in] 列ブロック行列（m x p_loc）
    real(kdouble), intent(in) :: A_loc(:,:)
    !> [in] 右辺ベクトル（m、全ランク複製）
    real(kdouble), intent(in) :: b(:)
    !> [out] 解ベクトル（p_loc、自ランク保持列に対応）
    real(kdouble), intent(out) :: w_loc(:)
    !> [in] MPI コミュニケータ
    integer(kint), intent(in) :: comm
    !> [in] scalapack コミュニケータ
    integer(kint), intent(in) :: scalapack_comm
    integer(kint) :: my_col, my_row, n_col, n_row
    integer(kint) :: p, lld, r, i, j, g, lr, info
    integer(kint) :: desc_G(9), desc_c(9)
    integer(kint), allocatable :: p_list(:), displs(:), ipiv(:)
    integer(kint), allocatable :: counts_w(:), displs_w(:)
    real(kdouble), allocatable :: A_cyc(:,:), G_loc(:,:), blk(:)
    real(kdouble), allocatable :: c_loc(:), w_cyc(:), w_all(:)

    integer :: numroc
    external :: numroc

    call blacs_gridinfo(scalapack_comm, n_row, n_col, my_row, my_col)

    !# 各ランクの列数の共有と全体列数の取得
    call monolis_alloc_I_1d(p_list, n_row)
    call monolis_alloc_I_1d(displs, n_row)
    call monolis_allgather_I1(p_loc, p_list, comm)

    p = 0
    do r = 1, n_row
      displs(r) = p
      p = p + p_list(r)
    enddo

    if(p == 0)then
      call monolis_dealloc_I_1d(p_list)
      call monolis_dealloc_I_1d(displs)
      return
    endif

    !# ブロックサイクリック分散（ブロック幅 1）で自ランクが保持する行数
    lld = numroc(p, 1, my_row, 0, n_row)

    call monolis_alloc_R_2d(A_cyc, m, max(1, lld))
    call monolis_alloc_R_1d(blk, m*maxval(p_list))

    !# パス 1: 自ランク保持行に対応する列ベクトルの取得
    do r = 0, n_row - 1
      if(p_list(r+1) == 0) cycle
      if(my_row == r)then
        do i = 1, p_loc
          do j = 1, m
            blk(m*(i-1) + j) = A_loc(j,i)
          enddo
        enddo
      endif
      call monolis_bcast_R(m*p_list(r+1), blk, r, comm)
      do i = 1, p_list(r+1)
        g = displs(r+1) + i
        if(mod(g - 1, n_row) /= my_row) cycle
        lr = (g - 1)/n_row + 1
        do j = 1, m
          A_cyc(j,lr) = blk(m*(i-1) + j)
        enddo
      enddo
    enddo

    !# パス 2: Gram 行列（自ランク保持行 lld x p）の組み立て
    call monolis_alloc_R_2d(G_loc, max(1, lld), p)

    do r = 0, n_row - 1
      if(p_list(r+1) == 0) cycle
      if(my_row == r)then
        do i = 1, p_loc
          do j = 1, m
            blk(m*(i-1) + j) = A_loc(j,i)
          enddo
        enddo
      endif
      call monolis_bcast_R(m*p_list(r+1), blk, r, comm)
      if(lld > 0)then
        call dgemm("T", "N", lld, p_list(r+1), m, 1.0d0, A_cyc, m, blk, m, &
          & 0.0d0, G_loc(1,displs(r+1)+1), max(1, lld))
      endif
    enddo

    !# 右辺ベクトル（自ランク保持行）
    call monolis_alloc_R_1d(c_loc, max(1, lld))
    if(lld > 0)then
      call dgemv("T", m, lld, 1.0d0, A_cyc, m, b, 1, 0.0d0, c_loc, 1)
    endif

    !# PDGESV による分散求解
    call descinit(desc_G, p, p, 1, 1, 0, 0, scalapack_comm, max(1, lld), info)
    call descinit(desc_c, p, 1, 1, 1, 0, 0, scalapack_comm, max(1, lld), info)

    call monolis_alloc_I_1d(ipiv, lld + 1)

    call pdgesv(p, 1, G_loc, 1, 1, desc_G, ipiv, c_loc, 1, 1, desc_c, info)

    if(info /= 0)then
      call monolis_std_error_string("monolis_scalapack_gram_gesv_R: PDGESV")
      call monolis_std_error_stop()
    endif

    !# サイクリック分散された解を列ブロック分散に再配置
    call monolis_alloc_I_1d(counts_w, n_row)
    call monolis_alloc_I_1d(displs_w, n_row)

    do r = 0, n_row - 1
      counts_w(r+1) = numroc(p, 1, r, 0, n_row)
    enddo
    do r = 2, n_row
      displs_w(r) = displs_w(r-1) + counts_w(r-1)
    enddo

    call monolis_alloc_R_1d(w_cyc, max(1, lld))
    call monolis_alloc_R_1d(w_all, p)

    do i = 1, lld
      w_cyc(i) = c_loc(i)
    enddo
    call monolis_allgather_V_R(lld, w_cyc, w_all, counts_w, displs_w, comm)

    do i = 1, p_loc
      g = displs(my_row+1) + i
      r = mod(g - 1, n_row)
      lr = (g - 1)/n_row + 1
      w_loc(i) = w_all(displs_w(r+1) + lr)
    enddo

    call monolis_dealloc_I_1d(p_list)
    call monolis_dealloc_I_1d(displs)
    call monolis_dealloc_I_1d(ipiv)
    call monolis_dealloc_I_1d(counts_w)
    call monolis_dealloc_I_1d(displs_w)
    call monolis_dealloc_R_2d(A_cyc)
    call monolis_dealloc_R_2d(G_loc)
    call monolis_dealloc_R_1d(blk)
    call monolis_dealloc_R_1d(c_loc)
    call monolis_dealloc_R_1d(w_cyc)
    call monolis_dealloc_R_1d(w_all)
  end subroutine monolis_scalapack_gram_gesv_R

!  !> @ingroup wrapper
!  !> PDGETRF 関数（実数型、LU分解）
!  subroutine monolis_scalapack_getrf_R(N_loc, N, A, ipiv, comm, scalapack_comm)
!    implicit none
!    !> [in] 行列の大きさ（ローカル行数）
!    integer(kint), intent(in) :: N_loc
!    !> [in] 行列の大きさ（全体のサイズ N x N）
!    integer(kint), intent(in) :: N
!    !> [in,out] 入力行列（N_loc x N）、出力はLU分解後の行列
!    real(kdouble), intent(inout) :: A(:,:)
!    !> [out] ピボット情報（N_loc）
!    integer(kint), intent(inout) :: ipiv(:)
!    !> [in] Scalapack コミュニケータ
!    integer(kint), intent(in) :: comm
!    integer(kint), intent(in) :: scalapack_comm
!    integer(kint) :: NB, desc_A(9)
!    integer(kint) :: lld_A
!    integer(kint) :: info
!    integer(kint) :: my_col, my_row, n_col, n_row
!    integer, external :: numroc
!
!    call blacs_gridinfo(scalapack_comm, n_row, n_col, my_row, my_col)
!
!    NB = 1
!    desc_A = 0
!    lld_A = max(1, numroc(n, nb, my_row, 0, n_row))
!    call descinit(desc_A, N, N, NB, NB, 0, 0, scalapack_comm, lld_A, info)
!
!    call pdgetrf(N, N, A, 1, 1, desc_A, ipiv, info)
!  end subroutine monolis_scalapack_getrf_R

  !> @ingroup wrapper
  !> PDGETRS 関数（実数型、LU分解による線形方程式の求解）
  subroutine monolis_scalapack_getrs_R(N_loc, N, NRHS, A, ipiv, B, comm, scalapack_comm)
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
    integer(kint), intent(inout) :: ipiv(:)
    !> [in,out] 右辺ベクトル（N_loc x NRHS）、出力は解ベクトル
    real(kdouble), intent(inout) :: B(:,:)
    !> [in] コミュニケータ
    integer(kint), intent(in) :: comm
    !> [in] Scalapack コミュニケータ
    integer(kint), intent(in) :: scalapack_comm
    integer(kint) :: NB, desc_A(9), desc_B(9)
    integer(kint) :: lld_A, lld_B
    integer(kint) :: i, j, info, comm_size, N_loc_max, N_fix
    integer(kint) :: my_col, my_row, n_col, n_row
    integer, external :: numroc
    real(kdouble), allocatable :: A_temp(:,:)
    real(kdouble), allocatable :: B_temp(:,:)

    call blacs_gridinfo(scalapack_comm, n_row, n_col, my_row, my_col)

    !# N の取得（各領域の行数を統一）
    comm_size = monolis_mpi_get_local_comm_size(comm)
    N_loc_max = N_loc
    call monolis_allreduce_I1(N_loc_max, monolis_mpi_max, comm)
    N_fix = comm_size*N_loc

    call monolis_alloc_R_2d(A_temp, N_loc_max, N_fix)
    call monolis_alloc_R_2d(B_temp, N_loc_max, NRHS)

!$omp parallel default(none) &
!$omp & shared(A, A_temp, B, B_temp) &
!$omp & firstprivate(N_fix, N_loc_max, NRHS) &
!$omp & private(i, j)
!$omp do collapse(2)
    do i = 1, N_fix
    do j = 1, N_loc_max
      A_temp(j,i) = A(j,i)
    enddo
    enddo
!$omp end do

!$omp do collapse(2)
    do i = 1, NRHS
    do j = 1, N_loc_max
      B_temp(j,i) = B(j,i)
    enddo
    enddo
!$omp end do
!$omp end parallel

    NB = 1
    desc_A = 0
    desc_B = 0
    lld_A = max(1, numroc(N, NB, my_row, 0, n_row))
    lld_B = max(1, numroc(N, NB, my_row, 0, n_row))

    call descinit(desc_A, N, N, NB, NB, 0, 0, scalapack_comm, lld_A, info)
    call descinit(desc_B, N, NRHS, NB, NB, 0, 0, scalapack_comm, lld_B, info)

    call pdgetrf(N, N, A_temp, 1, 1, desc_A, ipiv, info)
    call pdgetrs("N", N, NRHS, A_temp, 1, 1, desc_A, ipiv, B_temp, 1, 1, desc_B, info)

    call getrs_R_update_X(N_loc_max, NRHS, B_temp, n_row, comm)

    do i = 1, NRHS
    do j = 1, N_loc
      B(j,i) = B_temp(j,i)
    enddo
    enddo
  end subroutine monolis_scalapack_getrs_R

  !> @ingroup wrapper
  !> 解ベクトルのアップデート（実数型、解情報の更新）
  !> ScaLAPACK の block-cyclic 分散から連続した 1 次元分散に再配置
  subroutine getrs_R_update_X(N_loc, NRHS, B, n_row, comm)
    implicit none
    !> [in] 行列の全体サイズ
    integer(kint), intent(in) :: N_loc
    !> [in] 右辺ベクトル数
    integer(kint), intent(in) :: NRHS
    !> [in,out] 解ベクトル
    real(kdouble), intent(inout) :: B(:,:)
    !> [in] プロセス数
    integer(kint), intent(in) :: n_row
    !> [in] コミュニケータ
    integer(kint), intent(in) :: comm
    integer(kint) :: i, j, k, in, comm_size
    integer(kint), allocatable :: scounts(:), sdispls(:)
    integer(kint), allocatable :: rcounts(:), rdispls(:)
    real(kdouble), allocatable :: sendbuf(:), recvbuf(:)

    call monolis_alloc_R_1d(sendbuf, N_loc * NRHS)
    call monolis_alloc_R_1d(recvbuf, N_loc * NRHS)
    call monolis_alloc_I_1d(scounts, n_row)
    call monolis_alloc_I_1d(sdispls, n_row)
    call monolis_alloc_I_1d(rcounts, n_row)
    call monolis_alloc_I_1d(rdispls, n_row)

    comm_size = monolis_mpi_get_local_comm_size(comm)

    scounts = NRHS*N_loc/comm_size
    rcounts = NRHS*N_loc/comm_size
    do i = 1, n_row
      sdispls(i) = (i - 1) * NRHS*N_loc/comm_size
      rdispls(i) = (i - 1) * NRHS*N_loc/comm_size
    enddo

    do i = 1, N_loc
      k = mod(i - 1, comm_size) + (i - 1)/comm_size + 1
      do j = 1, NRHS
        sendbuf(NRHS*(i-1) + j) = B(k, j)
      enddo
    enddo

    call monolis_alltoall_V_R(sendbuf, scounts, sdispls, recvbuf, rcounts, rdispls, comm)

    do i = 1, N_loc
      k = mod(i - 1, comm_size) + (i - 1)/comm_size + 1
      do j = 1, NRHS
        B(k, j) = recvbuf(NRHS*(i-1) + j)
      enddo
    enddo
  end subroutine getrs_R_update_X
end module mod_monolis_scalapack
