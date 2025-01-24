!> 固有値ソルバモジュール
module mod_monolis_eigen_solver
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_eigen_lanczos
  use mod_monolis_spmat_handler_util
  use mod_monolis_scalapack

  implicit none

contains

  !> @ingroup eigen
  !> Lanczos 法（順反復、最大固有値、実数型）
  subroutine monolis_eigen_standard_lanczos_R( &
    & monolis, monoCOM, n_get_eigen, ths, maxiter, val, vec, is_bc)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in,out] 取得固有値数
    integer(kint), intent(inout) :: n_get_eigen
    !> [in] 収束判定閾値
    real(kdouble), intent(in) :: ths
    !> [in] 最大反復回数
    integer(kint), intent(in) :: maxiter
    !> [out] 固有値
    real(kdouble), intent(out) :: val(:)
    !> [out] 固有ベクトル
    real(kdouble), intent(out) :: vec(:,:)
    !> [in] Dirhchlet 境界条件判定フラグ
    logical, intent(in) :: is_bc(:)

    if(monoCOM%comm_size > 1) monolis%MAT%N = monoCOM%n_internal_vertex

    call monolis_eigen_standard_lanczos_R_main( &
      & monolis%PRM, monoCOM, monolis%MAT, n_get_eigen, ths, maxiter, val, vec, is_bc)
  end subroutine monolis_eigen_standard_lanczos_R

  !> @ingroup eigen
  !> Lanczos 法（逆反復、最小固有値、実数型）
  subroutine monolis_eigen_inverted_standard_lanczos_R( &
    & monolis, monoCOM, n_get_eigen, ths, maxiter, val, vec, is_bc)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in,out] 取得固有値数
    integer(kint), intent(inout) :: n_get_eigen
    !> [in] 収束判定閾値
    real(kdouble), intent(in) :: ths
    !> [in] 最大反復回数
    integer(kint), intent(in) :: maxiter
    !> [out] 固有値
    real(kdouble), intent(out) :: val(:)
    !> [out] 固有ベクトル
    real(kdouble), intent(out) :: vec(:,:)
    !> [in] Dirhchlet 境界条件判定フラグ
    logical, intent(in) :: is_bc(:)

    if(monoCOM%comm_size > 1) monolis%MAT%N = monoCOM%n_internal_vertex

    call monolis_eigen_inverted_standard_lanczos_R_main( &
      & monolis%PRM, monoCOM, monolis%MAT, monolis%PREC, n_get_eigen, ths, maxiter, val, vec, is_bc)
  end subroutine monolis_eigen_inverted_standard_lanczos_R

  !> @ingroup eigen
  !> 条件数の推定（Scalapack 利用、実数型、非対称対応）
  subroutine monolis_get_condition_number_R(monolis, monoCOM, condition_number)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [out] 条件数の推定値
    real(kdouble), intent(out) :: condition_number
    integer(kint) :: scalapack_comm
    integer(kint) :: N, NT, NDOF, i
    real(kdouble) :: min, max
    real(kdouble), allocatable :: dense(:,:)
    real(kdouble), allocatable :: S(:,:), V(:), D(:,:)

    call monolis_convert_sparse_matrix_to_dense_matrix_R(monolis%MAT, monoCOM, dense)

    call monolis_scalapack_comm_initialize(monoCOM%comm, scalapack_comm)

    N = monolis%MAT%N
    if(monolis_mpi_get_local_comm_size(monoCOM%comm) > 1) N = monoCOM%n_internal_vertex
    NT = N
    NDOF = monolis%MAT%NDOF
    call monolis_allreduce_I1(NT, monolis_mpi_sum, monoCOM%comm)

    call monolis_alloc_R_2d(S, NDOF*N, NDOF*NT)
    call monolis_alloc_R_1d(V, NDOF*NT)
    call monolis_alloc_R_2d(D, NDOF*NT, NDOF*NT)

    call monolis_scalapack_gesvd_R(NDOF*N, NDOF*NT, dense, S, V, D, monoCOM%comm, scalapack_comm)
    
    min = 1.0d300
    max = 0.0d0
    do i = 1, NDOF*NT
      if(min > dabs(V(i)) .and. dabs(V(i)) > 0.0d0) min = dabs(V(i))
      if(max < dabs(V(i)) .and. dabs(V(i)) > 0.0d0) max = dabs(V(i))
    enddo

    condition_number = max/min

    call monolis_scalapack_comm_finalize(scalapack_comm)
  end subroutine monolis_get_condition_number_R

end module mod_monolis_eigen_solver

