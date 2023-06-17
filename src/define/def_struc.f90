!> monolis 構造体の定義
module mod_monolis_def_struc
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_solver

  implicit none

  !> monolis 構造体
  type monolis_structure
    !> パラメータ構造体
    type(monolis_prm) :: PRM
    !> 行列構造体
    type(monolis_mat) :: MAT
    !> 前処理構造体
    type(monolis_mat) :: PREC
  end type monolis_structure

contains

  !> @ingroup def_init
  !> monolis ライブラリの初期化処理処理
  subroutine monolis_global_initialize()
    implicit none
    call monolis_mpi_initialize()
  end subroutine monolis_global_initialize

  !> @ingroup def_init
  !> monolis ライブラリの終了処理処理
  subroutine monolis_global_finalize()
    implicit none
    call monolis_mpi_finalize()
  end subroutine monolis_global_finalize

  !> @ingroup def_init
  !> monolis 構造体の初期化処理
  subroutine monolis_initialize(monolis)
    implicit none
    !> [out] monolis 構造体
    type(monolis_structure), intent(out) :: monolis

    call monolis_prm_initialize(monolis%PRM)
    call monolis_mat_initialize(monolis%MAT)
    call monolis_mat_initialize(monolis%PREC)
  end subroutine monolis_initialize

  !> @ingroup def_init
  !> monoils 構造体の終了処理
  subroutine monolis_finalize(monolis)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis

    call monolis_prm_finalize(monolis%PRM)
    call monolis_mat_finalize(monolis%MAT)
    call monolis_mat_finalize(monolis%PREC)
  end subroutine monolis_finalize
end module mod_monolis_def_struc
