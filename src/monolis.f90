!> monolis モジュールファイル
!> @details 全ての monolis モジュールをまとめたモジュールファイル
module mod_monolis
  use mod_monolis_utils
  use mod_gedatsu
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_def_solver
  use mod_monolis_def_solver_util

  !> @defgroup def_init 初期化・終了処理
  !> 初期化・終了処理に関連する関数グループ

  !> @defgroup dev_def_init 開発者用：初期化・終了処理
  !> 初期化・終了処理に関連する関数グループ（開発者用）
end module mod_monolis
