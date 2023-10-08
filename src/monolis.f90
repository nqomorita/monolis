!> monolis モジュールファイル
!> @details monolis_solver、gedatsu、monolis_utils モジュールをまとめたモジュールファイル
module mod_monolis
  use mod_monolis_solver
  use mod_monolis_utils
  use mod_gedatsu
end module mod_monolis

  subroutine interface_monolis_solve_main_R(monoPRM, monoCOM, monoMAT, monoPREC)
    use mod_monolis_utils
    use mod_monolis_def_solver
    use mod_monolis_def_mat
    use mod_monolis_solve
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in,out] 行列構造体
    type(monolis_mat), intent(inout) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC

    call monolis_solve_main_R(monoPRM, monoCOM, monoMAT, monoPREC)
  end subroutine interface_monolis_solve_main_R

  subroutine interface_monolis_solve_main_C(monoPRM, monoCOM, monoMAT, monoPREC)
    use mod_monolis_utils
    use mod_monolis_def_solver
    use mod_monolis_def_mat
    use mod_monolis_solve
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in,out] 行列構造体
    type(monolis_mat), intent(inout) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC

    call monolis_solve_main_C(monoPRM, monoCOM, monoMAT, monoPREC)
  end subroutine interface_monolis_solve_main_C
