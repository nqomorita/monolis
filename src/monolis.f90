!> monolis モジュールファイル
!> @details 全ての monolis モジュールをまとめたモジュールファイル
module mod_monolis
  use mod_monolis_utils
  use mod_gedatsu
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_def_solver
  use mod_monolis_def_solver_util
  use mod_monolis_spmat_handler_util
  use mod_monolis_spmat_handler
  use mod_monolis_spmat_nonzero_pattern_util
  use mod_monolis_spmat_nonzero_pattern
  use mod_monolis_linalg
  use mod_monolis_converge
  use mod_monolis_matvec
  use mod_monolis_vec_util
  use mod_monolis_solve

  !> @defgroup def_init 初期化・終了処理関数群
  !> 初期化・終了処理に関連する関数グループ

  !> @defgroup linalg 代数演算関数群
  !> 代数演算に関連する関数グループ

  !> @defgroup dev_def_init 開発者用：初期化・終了処理関数群
  !> 初期化・終了処理に関連する関数グループ（開発者用）

  !> @defgroup dev_linalg 開発者用：代数演算関数群
  !> 代数演算に関連する関数グループ（開発者用）
end module mod_monolis
