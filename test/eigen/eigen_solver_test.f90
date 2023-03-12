!> 固有値ソルバテストモジュール
module mod_monolis_eigen_solver_test
  use mod_monolis

  implicit none

contains

  subroutine monolis_eigen_solver_test()
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
    real(kdouble) :: val(6)
    !> 固有ベクトル
    real(kdouble) :: vec(6,6)
    !> Dirhchlet 境界条件判定フラグ
    logical :: is_bc(6)

    call monolis_eigen_standard_lanczos_R &
      & (monolis, n_get_eigen, ths, maxiter, val, vec, is_bc)

    call monolis_eigen_inverted_standard_lanczos_R &
      & (monolis, n_get_eigen, ths, maxiter, val, vec, is_bc)
  end subroutine monolis_eigen_solver_test
end module mod_monolis_eigen_solver_test

