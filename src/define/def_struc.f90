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
    !> 通信テーブル構造体
    type(monolis_com) :: COM
    !> 行列構造体
    type(monolis_mat) :: MAT
  end type monolis_structure

contains

  !> monolis ライブラリの初期化処理処理
  subroutine monolis_global_initialize()
    implicit none
    call monolis_mpi_initialize()
  end subroutine monolis_global_initialize

  !> monolis ライブラリの終了処理処理
  subroutine monolis_global_finalize()
    implicit none
    call monolis_mpi_finalize()
  end subroutine monolis_global_finalize

  !> monolis 構造体の初期化処理
  subroutine monolis_initialize(monolis, fname_in)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> 読込ディレクトリ
    character(*) :: fname_in

    call monolis_prm_initialize(monolis%PRM)
    call monolis_com_initialize(monolis%COM)
    call monolis_mat_initialize(monolis%MAT)
    !call monolis_com_input_comm_table(monolis%COM, fname_in)
  end subroutine monolis_initialize

  !> monoils 構造体の終了処理
  subroutine monolis_finalize(monolis)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis

    call monolis_prm_finalize(monolis%PRM)
    call monolis_com_finalize(monolis%COM)
    call monolis_mat_finalize(monolis%MAT)
  end subroutine monolis_finalize
end module mod_monolis_def_struc
