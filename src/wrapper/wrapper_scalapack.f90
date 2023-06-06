!> ScaLAPACK ラッパーモジュール
module mod_monolis_scalapack
  use mod_monolis_utils
  use mod_monolis_vec_util

  implicit none

contains

  !> PDGESVD 関数（実数型）
  subroutine monolis_scalapack_gesvd_R(N_loc, M, A, S, V, D, comm)
    implicit none
    !> 行列の大きさ（行数 N）
    integer(kint), intent(in) :: N_loc
    !> 行列の大きさ（列数 M）
    integer(kint), intent(in) :: M
    !> 入力行列（N_loc x M）
    real(kdouble) :: A(:,:)
    !> 左特異行列（N_loc x P）
    real(kdouble) :: S(:,:)
    !> 特異値（P）
    real(kdouble) :: V(:)
    !> 右特異行列（P x M）
    real(kdouble) :: D(:,:)
    !> コミュニケータ
    integer(kint) :: N, comm
    integer(kint) :: scalapack_comm
    integer(kint) :: NB, P, desc_A(9), desc_S(9), desc_D(9)
    integer(kint) :: lld_A, lld_S, lld_D
    integer(kint) :: NW, info
    integer(kint) :: my_col, my_row, n_col, n_row
    real(kdouble), allocatable :: W(:)

    integer :: numroc
    external :: numroc

    !# scalapack コミュニケータの初期化処理
    n_row = monolis_mpi_get_local_comm_size(comm)
    n_col = 1

    call blacs_get(0, 0, scalapack_comm)
    call blacs_gridinit(scalapack_comm, "r", n_row, n_col)
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
    lld_D = numroc(M, NB, my_row, 0, n_row)

    call descinit(desc_A, N, M, NB, NB, 0, 0, scalapack_comm, lld_A, info)
    call descinit(desc_S, N, P, NB, NB, 0, 0, scalapack_comm, lld_S, info)
    call descinit(desc_D, P, M, NB, NB, 0, 0, scalapack_comm, lld_D, info)

    !# 一次ベクトルの大きさ取得
    call monolis_alloc_R_1d(W, 1)

    call pdgesvd("V", "V", N, M, &
      & A, 1, 1, desc_A, &
      & V, &
      & S, 1, 1, desc_S, &
      & D, 1, 1, desc_D, &
      & W, -1, info)

    NW = int(W(1))
    call monolis_dealloc_R_1d(W)
    call monolis_alloc_R_1d(W, NW)

    !# 特異値分解
    call pdgesvd("V", "V", N, M, &
      & A, 1, 1, desc_A, &
      & V, &
      & S, 1, 1, desc_S, &
      & D, 1, 1, desc_D, &
      & W, NW, info)

    !# scalapack コミュニケータの終了処理
    call blacs_gridexit(scalapack_comm)
  end subroutine monolis_scalapack_gesvd_R

end module mod_monolis_scalapack
