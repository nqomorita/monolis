!> ソルバパラメータの設定関数群
module mod_monolis_def_solver_util
  use mod_monolis_utils
  use mod_monolis_def_struc

  implicit none

contains

  !> ソルバの設定
  subroutine monolis_set_method(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    integer(kint) :: param
    monolis%PRM%monolis_prm_Iarray(monolis_prm_method) = param
  end subroutine monolis_set_method

  !> 前処理の設定
  subroutine monolis_set_precond(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    integer(kint) :: param
    monolis%PRM%monolis_prm_Iarray(monolis_prm_precond) = param
  end subroutine monolis_set_precond

  !> 最大反復回数の設定
  subroutine monolis_set_maxiter(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    integer(kint) :: param
    monolis%PRM%monolis_prm_Iarray(monolis_prm_max_iter) = param
  end subroutine monolis_set_maxiter

  !> 現在の反復回数の設定
  subroutine monolis_get_converge_iter(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    integer(kint) :: param
    param = monolis%PRM%monolis_prm_Iarray(monolis_prm_cur_iter)
  end subroutine monolis_get_converge_iter

  !> エラー番号の取得
  subroutine monolis_get_error_tag(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    integer(kint) :: param
    param = monolis%PRM%monolis_prm_Iarray(monolis_prm_ierr)
  end subroutine monolis_get_error_tag

  !> スケーリングの有無の設定
  subroutine monolis_set_scaling(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    logical :: param
    integer(kint) :: iparam

    monolis%PRM%monolis_prm_Iarray(monolis_prm_is_scaling) = monolis_conv_L2I(param)
  end subroutine monolis_set_scaling

  !> リオーダリングの有無の設定
  subroutine monolis_set_reordering(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    logical :: param
    integer(kint) :: iparam

    monolis%PRM%monolis_prm_Iarray(monolis_prm_is_reordering) = monolis_conv_L2I(param)
  end subroutine monolis_set_reordering

  !> 解ベクトル初期化の有無の設定
  subroutine monolis_set_init_x(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    logical :: param
    integer(kint) :: iparam

    monolis%PRM%monolis_prm_Iarray(monolis_prm_is_init_x) = monolis_conv_L2I(param)
  end subroutine monolis_set_init_x

  !> 対称行列向け処理の有無の設定
  subroutine monolis_set_sym_matrix(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    logical :: param
    integer(kint) :: iparam

    monolis%PRM%monolis_prm_Iarray(monolis_prm_is_sym_matrix) = monolis_conv_L2I(param)
  end subroutine monolis_set_sym_matrix

  !> デバッグ出力の有無の設定
  subroutine monolis_set_debug(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    logical :: param
    integer(kint) :: iparam

    monolis%PRM%monolis_prm_Iarray(monolis_prm_is_debug) = monolis_conv_L2I(param)
  end subroutine monolis_set_debug

  !> 詳細な計算時間測定の有無の設定
  subroutine monolis_set_performance_measurement(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    logical :: param
    integer(kint) :: iparam

    monolis%PRM%monolis_prm_Iarray(monolis_prm_is_measurement) = monolis_conv_L2I(param)
  end subroutine monolis_set_performance_measurement

  !> 行列対角成分確認の有無の設定
  subroutine monolis_set_check_diag(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    logical :: param
    integer(kint) :: iparam

    monolis%PRM%monolis_prm_Iarray(monolis_prm_is_check_diag) = monolis_conv_L2I(param)
  end subroutine monolis_set_check_diag

  !> 前処理情報保存の有無の設定
  subroutine monolis_set_prec_stored(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    logical :: param
    integer(kint) :: iparam

    monolis%PRM%monolis_prm_Iarray(monolis_prm_is_prec_stored) = monolis_conv_L2I(param)
  end subroutine monolis_set_prec_stored

  !> 反復回数と残差履歴の表示の設定
  subroutine monolis_show_iterlog(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    logical :: param
    integer(kint) :: iparam

    monolis%PRM%monolis_prm_Iarray(monolis_prm_show_iterlog) = monolis_conv_L2I(param)
  end subroutine monolis_show_iterlog

  !> 詳細な計算時間の表示の設定
  subroutine monolis_show_timelog(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    logical :: param
    integer(kint) :: iparam

    monolis%PRM%monolis_prm_Iarray(monolis_prm_show_time) = monolis_conv_L2I(param)
  end subroutine monolis_show_timelog

  !> ソルバ収束後のサマリの表示の設定
  subroutine monolis_show_summary(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    logical :: param
    integer(kint) :: iparam

    monolis%PRM%monolis_prm_Iarray(monolis_prm_show_summary) = monolis_conv_L2I(param)
  end subroutine monolis_show_summary

  !> 計算時間の統計的処理結果の表示の設定
  subroutine monolis_show_timelog_statistics(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    logical :: param
    integer(kint) :: iparam

    monolis%PRM%monolis_prm_Iarray(monolis_prm_show_time_statistics) = monolis_conv_L2I(param)
  end subroutine monolis_show_timelog_statistics

  !> 収束判定閾値の設定
  subroutine monolis_set_tolerance(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    real(kdouble) :: param
    monolis%PRM%monolis_prm_Rarray(monolis_prm_tol) = param
  end subroutine monolis_set_tolerance

  !> 現在の残差の取得
  subroutine monolis_get_converge_residual(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    real(kdouble) :: param
    param = monolis%PRM%monolis_prm_Rarray(monolis_prm_cur_resid)
  end subroutine monolis_get_converge_residual

  !> ソルバの全計算時間の取得
  subroutine monolis_get_time_solver(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    real(kdouble) :: param
    param = monolis%PRM%monolis_prm_Rarray(monolis_time_sol)
  end subroutine monolis_get_time_solver

  !> 前処理時間（生成時間）の取得
  subroutine monolis_get_time_preparing(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    real(kdouble) :: param
    param = monolis%PRM%monolis_prm_Rarray(monolis_time_prep)
  end subroutine monolis_get_time_preparing

  !> 疎行列ベクトル積時間の取得
  subroutine monolis_get_time_spmv(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    real(kdouble) :: param
    param = monolis%PRM%monolis_prm_Rarray(monolis_time_spmv)
  end subroutine monolis_get_time_spmv

  !> ベクトル内積時間の取得
  subroutine monolis_get_time_inner_product(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    real(kdouble) :: param
    param = monolis%PRM%monolis_prm_Rarray(monolis_time_dotp)
  end subroutine monolis_get_time_inner_product

  !> 前処理時間（適用時間）の取得
  subroutine monolis_get_time_precondition(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    real(kdouble) :: param
    param = monolis%PRM%monolis_prm_Rarray(monolis_time_prec)
  end subroutine monolis_get_time_precondition

  !> ベクトル内積の通信時間の取得
  subroutine monolis_get_time_comm_inner_product(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    real(kdouble) :: param
    param = monolis%PRM%monolis_prm_Rarray(monolis_time_comm_dotp)
  end subroutine monolis_get_time_comm_inner_product

  integer(kint), parameter :: monolis_time_comm_dotp = 8

  !> 疎行列ベクトル積の通信時間の取得
  subroutine monolis_get_time_comm_spmv(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    real(kdouble) :: param
    param = monolis%PRM%monolis_prm_Rarray(monolis_time_comm_spmv)
  end subroutine monolis_get_time_comm_spmv
end module mod_monolis_def_solver_util
