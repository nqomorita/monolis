!> 固有値ソルバモジュール
module mod_monolis_eigen_solver
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_eigen_lanczos

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
end module mod_monolis_eigen_solver

