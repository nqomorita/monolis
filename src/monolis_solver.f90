!> monolis モジュールファイル
!> @details 全ての monolis モジュールをまとめたモジュールファイル
module mod_monolis_solver
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_def_solver
  use mod_monolis_def_solver_util
  use mod_monolis_spmat_handler_util
  use mod_monolis_spmat_handler
  use mod_monolis_spmat_nonzero_pattern_util
  use mod_monolis_spmat_nonzero_pattern
  use mod_monolis_spmat_copy
  use mod_monolis_inner_product
  use mod_monolis_converge
  use mod_monolis_matvec
  use mod_monolis_vec_util
  use mod_monolis_solve
  use mod_monolis_precond
  use mod_monolis_lapack
  use mod_monolis_scalapack
  use mod_monolis_eigen_solver
  use mod_monolis_eigen_lanczos
  use mod_monolis_eigen_lanczos_util

  !> @defgroup def_init 全体初期化・終了処理関数群
  !> 初期化・終了処理に関連する関数グループ

  !> @defgroup def_mat_init 疎行列構造体の初期化・終了処理関数群
  !> 疎行列構造体の初期化・終了処理に関連する関数グループ

  !> @defgroup param パラメータ設定関数群
  !> パラメータ設定に関連する関数グループ

  !> @defgroup linalg 代数演算関数群
  !> 代数演算に関連する関数グループ

  !> @defgroup matrix 疎行列演算関数群（行列値の設定・取得）
  !> 疎行列演算に関連する関数グループ

  !> @defgroup matrix_copy 疎行列演算関数群（行列情報のコピーとクリア）
  !> 疎行列演算に関連する関数グループ

  !> @defgroup nzpattern 疎行列演算関数群（疎行列パターン決定）
  !> 疎行列演算に関連する関数グループ

  !> @defgroup solver 線形ソルバ関数群
  !> 線形ソルバに関連する関数グループ

  !> @defgroup prec 前処理関数群
  !> 前処理に関連する関数グループ

  !> @defgroup eigen 固有値解析関数群
  !> 固有値計算に関連する関数グループ

  !> @defgroup wrapper 外部ライブラリ呼出関数群
  !> 外部ライブラリ呼出に関連する関数グループ

  !> @defgroup dev_linalg 開発者用：代数演算関数群
  !> 代数演算に関連する関数グループ（開発者用）

  !> @defgroup dev_matrix 開発者用：疎行列演算関数群
  !> 疎行列演算に関連する関数グループ（開発者用）

  !> @defgroup dev_solver 開発者用：線形ソルバ関数群
  !> 線形ソルバに関連する関数グループ（開発者用）
end module mod_monolis_solver
