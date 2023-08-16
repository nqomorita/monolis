!> ソルバパラメータの設定関数群
!# subroutine monolis_set_method(monolis, param)
!# subroutine monolis_set_precond(monolis, param)
!# subroutine monolis_set_maxiter(monolis, param)
!# subroutine monolis_get_converge_iter(monolis, param)
!# subroutine monolis_get_error_tag(monolis, param)
!# subroutine monolis_set_scaling(monolis, param)
!# subroutine monolis_set_reordering(monolis, param)
!# subroutine monolis_set_init_x(monolis, param)
!# subroutine monolis_set_sym_matrix(monolis, param)
!# subroutine monolis_set_debug(monolis, param)
!# subroutine monolis_set_performance_measurement(monolis, param)
!# subroutine monolis_set_check_diag(monolis, param)
!# subroutine monolis_set_prec_stored(monolis, param)
!# subroutine monolis_show_iterlog(monolis, param)
!# subroutine monolis_show_timelog(monolis, param)
!# subroutine monolis_show_summary(monolis, param)
!# subroutine monolis_show_timelog_statistics(monolis, param)
!# subroutine monolis_set_tolerance(monolis, param)
!# subroutine monolis_get_converge_residual(monolis, param)
!# subroutine monolis_get_time_solver(monolis, param)
!# subroutine monolis_get_time_preparing(monolis, param)
!# subroutine monolis_get_time_spmv(monolis, param)
!# subroutine monolis_get_time_inner_product(monolis, param)
!# subroutine monolis_get_time_precondition(monolis, param)
!# subroutine monolis_get_time_comm_inner_product(monolis, param)
!# subroutine monolis_get_time_comm_spmv(monolis, param)
module mod_monolis_def_solver_util
  use mod_monolis_utils
  use mod_monolis_def_struc

  implicit none

contains

  !> @ingroup param
  !> ソルバの設定
  subroutine monolis_set_method(monolis, param)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] パラメータ
    integer(kint), intent(in) :: param
    monolis%PRM%Iarray(monolis_prm_I_method) = param
  end subroutine monolis_set_method

  !> @ingroup param
  !> 前処理の設定
  subroutine monolis_set_precond(monolis, param)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] パラメータ
    integer(kint), intent(in) :: param
    monolis%PRM%Iarray(monolis_prm_I_precond) = param
  end subroutine monolis_set_precond

  !> @ingroup param
  !> 最大反復回数の設定
  subroutine monolis_set_maxiter(monolis, param)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] パラメータ
    integer(kint), intent(in) :: param
    monolis%PRM%Iarray(monolis_prm_I_max_iter) = param
  end subroutine monolis_set_maxiter

  !> @ingroup param
  !> 現在の反復回数の取得
  subroutine monolis_get_converge_iter(monolis, param)
    implicit none
    !> [in] monolis 構造体
    type(monolis_structure), intent(in) :: monolis
    !> [out] パラメータ
    integer(kint), intent(out) :: param
    param = monolis%PRM%Iarray(monolis_prm_I_cur_iter)
  end subroutine monolis_get_converge_iter

  !> @ingroup param
  !> エラー番号の取得
  subroutine monolis_get_error_tag(monolis, param)
    implicit none
    !> [in] monolis 構造体
    type(monolis_structure), intent(in) :: monolis
    !> [out] パラメータ
    integer(kint), intent(out) :: param
    param = monolis%PRM%Iarray(monolis_prm_I_ierr)
  end subroutine monolis_get_error_tag

  !> スケーリングの有無の設定
  !subroutine monolis_set_scaling(monolis, param)
  !  implicit none
    !> [in,out] monolis 構造体
  !  type(monolis_structure), intent(inout) :: monolis
    !> [in] パラメータ
  !  logical, intent(in) :: param
  !  integer(kint) :: iparam

  !  monolis%PRM%Iarray(monolis_prm_I_is_scaling) = monolis_conv_L2I(param)
  !end subroutine monolis_set_scaling

  !> リオーダリングの有無の設定
  !subroutine monolis_set_reordering(monolis, param)
  !  implicit none
    !> [in,out] monolis 構造体
  !  type(monolis_structure), intent(inout) :: monolis
    !> [in] パラメータ
  !  logical, intent(in) :: param
  !  integer(kint) :: iparam

  !  monolis%PRM%Iarray(monolis_prm_I_is_reordering) = monolis_conv_L2I(param)
  !end subroutine monolis_set_reordering

  !> @ingroup param
  !> 解ベクトル初期化の有無の設定
  subroutine monolis_set_init_x(monolis, param)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] パラメータ
    logical, intent(in) :: param
    integer(kint) :: iparam

    monolis%PRM%Iarray(monolis_prm_I_is_init_x) = monolis_conv_L2I(param)
  end subroutine monolis_set_init_x

  !> @ingroup param
  !> 対称行列向け処理の有無の設定
  subroutine monolis_set_sym_matrix(monolis, param)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] パラメータ
    logical, intent(in) :: param
    integer(kint) :: iparam

    monolis%PRM%Iarray(monolis_prm_I_is_sym_matrix) = monolis_conv_L2I(param)
  end subroutine monolis_set_sym_matrix

  !> @ingroup param
  !> デバッグ出力の有無の設定
  subroutine monolis_set_debug(monolis, param)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] パラメータ
    logical, intent(in) :: param
    integer(kint) :: iparam

    monolis%PRM%Iarray(monolis_prm_I_is_debug) = monolis_conv_L2I(param)
  end subroutine monolis_set_debug

  !> @ingroup param
  !> 詳細な計算時間測定の有無の設定
  subroutine monolis_set_performance_measurement(monolis, param)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] パラメータ
    logical, intent(in) :: param
    integer(kint) :: iparam

    monolis%PRM%Iarray(monolis_prm_I_is_measurement) = monolis_conv_L2I(param)
  end subroutine monolis_set_performance_measurement

  !> @ingroup param
  !> 行列対角成分確認の有無の設定
  subroutine monolis_set_check_diag(monolis, param)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] パラメータ
    logical, intent(in) :: param
    integer(kint) :: iparam

    monolis%PRM%Iarray(monolis_prm_I_is_check_diag) = monolis_conv_L2I(param)
  end subroutine monolis_set_check_diag

  !> @ingroup param
  !> 前処理情報保存の有無の設定
  subroutine monolis_set_prec_stored(monolis, param)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] パラメータ
    logical, intent(in) :: param
    integer(kint) :: iparam

    monolis%PRM%Iarray(monolis_prm_I_is_prec_stored) = monolis_conv_L2I(param)
  end subroutine monolis_set_prec_stored

  !> @ingroup param
  !> 反復回数と残差履歴の表示の設定
  subroutine monolis_show_iterlog(monolis, param)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] パラメータ
    logical, intent(in) :: param
    integer(kint) :: iparam

    monolis%PRM%Iarray(monolis_prm_I_show_iterlog) = monolis_conv_L2I(param)
  end subroutine monolis_show_iterlog

  !> @ingroup param
  !> 詳細な計算時間の表示の設定
  subroutine monolis_show_timelog(monolis, param)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] パラメータ
    logical, intent(in) :: param
    integer(kint) :: iparam

    monolis%PRM%Iarray(monolis_prm_I_show_time) = monolis_conv_L2I(param)
  end subroutine monolis_show_timelog

  !> @ingroup param
  !> ソルバ収束後のサマリの表示の設定
  subroutine monolis_show_summary(monolis, param)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] パラメータ
    logical, intent(in) :: param
    integer(kint) :: iparam

    monolis%PRM%Iarray(monolis_prm_I_show_summary) = monolis_conv_L2I(param)
  end subroutine monolis_show_summary

  !> @ingroup param
  !> 計算時間の統計的処理結果の表示の設定
  subroutine monolis_show_timelog_statistics(monolis, param)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] パラメータ
    logical, intent(in) :: param
    integer(kint) :: iparam

    monolis%PRM%Iarray(monolis_prm_I_show_time_statistics) = monolis_conv_L2I(param)
  end subroutine monolis_show_timelog_statistics

  !> @ingroup param
  !> 収束判定閾値の設定
  subroutine monolis_set_tolerance(monolis, param)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] パラメータ
    real(kdouble), intent(in) :: param
    monolis%PRM%Rarray(monolis_prm_R_tol) = param
  end subroutine monolis_set_tolerance

  !> @ingroup param
  !> 現在の残差の取得
  subroutine monolis_get_converge_residual(monolis, param)
    implicit none
    !> [in] monolis 構造体
    type(monolis_structure), intent(in) :: monolis
    !> [out] パラメータ
    real(kdouble), intent(out) :: param
    param = monolis%PRM%Rarray(monolis_prm_R_cur_resid)
  end subroutine monolis_get_converge_residual

  !> @ingroup param
  !> ソルバの全計算時間の取得
  subroutine monolis_get_time_solver(monolis, param)
    implicit none
    !> [in] monolis 構造体
    type(monolis_structure), intent(in) :: monolis
    !> [out] パラメータ
    real(kdouble), intent(out) :: param
    param = monolis%PRM%Rarray(monolis_R_time_sol)
  end subroutine monolis_get_time_solver

  !> @ingroup param
  !> 前処理時間（生成時間）の取得
  subroutine monolis_get_time_preparing(monolis, param)
    implicit none
    !> [in] monolis 構造体
    type(monolis_structure), intent(in) :: monolis
    !> [out] パラメータ
    real(kdouble), intent(out) :: param
    param = monolis%PRM%Rarray(monolis_R_time_prep)
  end subroutine monolis_get_time_preparing

  !> @ingroup param
  !> 疎行列ベクトル積時間の取得
  subroutine monolis_get_time_spmv(monolis, param)
    implicit none
    !> [in] monolis 構造体
    type(monolis_structure), intent(in) :: monolis
    !> [out] パラメータ
    real(kdouble), intent(out) :: param
    param = monolis%PRM%Rarray(monolis_R_time_spmv)
  end subroutine monolis_get_time_spmv

  !> @ingroup param
  !> ベクトル内積時間の取得
  subroutine monolis_get_time_inner_product(monolis, param)
    implicit none
    !> [in] monolis 構造体
    type(monolis_structure), intent(in) :: monolis
    !> [out] パラメータ
    real(kdouble), intent(out) :: param
    param = monolis%PRM%Rarray(monolis_R_time_dotp)
  end subroutine monolis_get_time_inner_product

  !> @ingroup param
  !> 前処理時間（適用時間）の取得
  subroutine monolis_get_time_precondition(monolis, param)
    implicit none
    !> [in] monolis 構造体
    type(monolis_structure), intent(in) :: monolis
    !> [out] パラメータ
    real(kdouble), intent(out) :: param
    param = monolis%PRM%Rarray(monolis_R_time_prec)
  end subroutine monolis_get_time_precondition

  !> @ingroup param
  !> ベクトル内積の通信時間の取得
  subroutine monolis_get_time_comm_inner_product(monolis, param)
    implicit none
    !> [in] monolis 構造体
    type(monolis_structure), intent(in) :: monolis
    !> [out] パラメータ
    real(kdouble), intent(out) :: param
    param = monolis%PRM%Rarray(monolis_R_time_comm_dotp)
  end subroutine monolis_get_time_comm_inner_product

  !> @ingroup param
  !> 疎行列ベクトル積の通信時間の取得
  subroutine monolis_get_time_comm_spmv(monolis, param)
    implicit none
    !> [in] monolis 構造体
    type(monolis_structure), intent(in) :: monolis
    !> [out] パラメータ
    real(kdouble), intent(out) :: param
    param = monolis%PRM%Rarray(monolis_R_time_comm_spmv)
  end subroutine monolis_get_time_comm_spmv

end module mod_monolis_def_solver_util
