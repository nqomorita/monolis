!> 固有値ソルバモジュール
module mod_monolis_eigen_solver
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_eigen_lanczos

  implicit none

contains

  !> Lanczos 法（順反復、最大固有値、実数型）
  subroutine monolis_eigen_standard_lanczos_R &
    & (monolis, n_get_eigen, ths, maxiter, val, vec, is_bc)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> 取得固有値数
    integer(kint) :: n_get_eigen
    !> 収束判定閾値
    real(kdouble) :: ths
    !> 最大反復回数
    integer(kint) :: maxiter
    !> 固有値
    real(kdouble) :: val(:)
    !> 固有ベクトル
    real(kdouble) :: vec(:,:)
    !> Dirhchlet 境界条件判定フラグ
    logical :: is_bc(:)

    call monolis_eigen_standard_lanczos_R_main( &
      & monolis%PRM, monolis%COM, monolis%MAT, n_get_eigen, ths, maxiter, val, vec, is_bc)
  end subroutine monolis_eigen_standard_lanczos_R

  !> Lanczos 法（逆反復、最小固有値、実数型）
  subroutine monolis_eigen_inverted_standard_lanczos_R &
    & (monolis, n_get_eigen, ths, maxiter, val, vec, is_bc)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> 取得固有値数
    integer(kint) :: n_get_eigen
    !> 収束判定閾値
    real(kdouble) :: ths
    !> 最大反復回数
    integer(kint) :: maxiter
    !> 固有値
    real(kdouble) :: val(:)
    !> 固有ベクトル
    real(kdouble) :: vec(:,:)
    !> Dirhchlet 境界条件判定フラグ
    logical :: is_bc(:)

    call monolis_eigen_inverted_standard_lanczos_R_main( &
      & monolis%PRM, monolis%COM, monolis%MAT, monolis%PREC, n_get_eigen, ths, maxiter, val, vec, is_bc)
  end subroutine monolis_eigen_inverted_standard_lanczos_R

end module mod_monolis_eigen_solver

