!> ソルバパラメータの定義
module mod_monolis_def_solver
  use mod_monolis_utils
  implicit none

  !> パラメータ：CG 法
  integer(kint), parameter :: monolis_iter_CG       = 1
  !> パラメータ：GropCG 法
  !integer(kint), parameter :: monolis_iter_GropCG   = 2
  !> パラメータ：Pipelined CG 法
  !integer(kint), parameter :: monolis_iter_PipeCG   = 3
  !> パラメータ：Pipelined CR 法
  !integer(kint), parameter :: monolis_iter_PipeCR   = 4
  !> パラメータ：BiCGSTAB 法
  integer(kint), parameter :: monolis_iter_BiCGSTAB = 5
  !> パラメータ：Pipelined BiCGSTAB 法
  !integer(kint), parameter :: monolis_iter_PipeBiCGSTAB = 6
  !> パラメータ：BiCGSTAB 法（前処理なし）
  !integer(kint), parameter :: monolis_iter_BiCGSTAB_noprec = 7
  !> パラメータ：Communication Avoiding BiCGSTAB 法（前処理なし）
  !integer(kint), parameter :: monolis_iter_CABiCGSTAB_noprec = 8
  !> パラメータ：Pipelined BiCGSTAB 法（前処理なし）
  !integer(kint), parameter :: monolis_iter_PipeBiCGSTAB_noprec = 9
  !> パラメータ：SOR 法
  !integer(kint), parameter :: monolis_iter_SOR      = 10
  !> パラメータ：Iterative Refinement 法
  !integer(kint), parameter :: monolis_iter_IR       = 11
  !> パラメータ：GMRES 法
  !integer(kint), parameter :: monolis_iter_GMRES    = 12

  !> パラメータ：前処理なし
  integer(kint), parameter :: monolis_prec_NONE   = 0
  !> パラメータ：対角スケーリング前処理
  integer(kint), parameter :: monolis_prec_DIAG   = 1
  !> パラメータ：ブロック ILU 前処理
  !integer(kint), parameter :: monolis_prec_ILU    = 2
  !> パラメータ：ブロック Jacobi 前処理
  !integer(kint), parameter :: monolis_prec_JACOBI = 3
  !> パラメータ：ブロック SOR 前処理
  integer(kint), parameter :: monolis_prec_SOR    = 4
  !> パラメータ：ブロック SAINV 前処理
  !integer(kint), parameter :: monolis_prec_SAINV  = 5
  !> パラメータ：ブロック RIF 前処理
  !integer(kint), parameter :: monolis_prec_RIF    = 6
  !> パラメータ：SPIKE 前処理
  !integer(kint), parameter :: monolis_prec_SPIKE  = 7
  !> パラメータ：LU 分解
  !integer(kint), parameter :: monolis_prec_LU     = 8
  !> パラメータ：LU 分解（MUMPS）
  integer(kint), parameter :: monolis_prec_MUMPS  = 9
  !> パラメータ：ROM 前処理
  !integer(kint), parameter :: monolis_prec_ROM    = 10
  !> パラメータ：LU 分解（MF 法）
  !integer(kint), parameter :: monolis_prec_MF     = 11
  !> パラメータ：ブロック LU 分解（MUMPS）
  integer(kint), parameter :: monolis_prec_MUMPS_LOCAL = 12

  character*24, dimension(12) :: monolis_str_iter = (/&
  & "CG                 ", &
  & "GropCG             ", &
  & "PipeCG             ", &
  & "PipeCR             ", &
  & "BiCGSTAB           ", &
  & "PipeBiCGSTAB       ", &
  & "BiCGSTAB_noprec    ", &
  & "CABiCGSTAB_noprec  ", &
  & "PipeBiCGSTAB_noprec", &
  & "SOR                ", &
  & "IR                 ", &
  & "GMRES              "/)

  character*24, dimension(0:12)  :: monolis_str_prec = (/&
  & "None  ", &
  & "Diag  ", &
  & "ILU   ", &
  & "Jacobi", &
  & "SOR   ", &
  & "SAINV ", &
  & "RIF   ", &
  & "SPIKE ", &
  & "LU    ", &
  & "MUMPS ", &
  & "ROM   ", &
  & "MF    ", &
  & "MUMPSL"/)

  !> 整数パラメータのサイズ
  integer(kint), parameter :: monolis_prm_Iarray_size = 100

  !> 実数パラメータのサイズ
  integer(kint), parameter :: monolis_prm_Rarray_size = 100

  !> パラメータ：ソルバ
  integer(kint), parameter :: monolis_prm_method = 1
  !> パラメータ：前処理
  integer(kint), parameter :: monolis_prm_precond = 2
  !> パラメータ：最大反復回数
  integer(kint), parameter :: monolis_prm_max_iter = 3
  !> パラメータ：現在の反復回数
  integer(kint), parameter :: monolis_prm_cur_iter = 4
  !> パラメータ：エラー番号
  integer(kint), parameter :: monolis_prm_ierr = 5
  !> パラメータ：スケーリングの有無
  integer(kint), parameter :: monolis_prm_is_scaling = 6
  !> パラメータ：リオーダリングの有無
  integer(kint), parameter :: monolis_prm_is_reordering = 7
  !> パラメータ：解ベクトル初期化の有無
  integer(kint), parameter :: monolis_prm_is_init_x = 8
  !> パラメータ：対称行列向け処理の有無
  integer(kint), parameter :: monolis_prm_is_sym_matrix = 9
  !> パラメータ：デバッグ出力の有無
  integer(kint), parameter :: monolis_prm_is_debug = 10
  !> パラメータ：詳細な計算時間測定の有無
  integer(kint), parameter :: monolis_prm_is_measurement = 11
  !> パラメータ：行列対角成分確認の有無
  integer(kint), parameter :: monolis_prm_is_check_diag = 12
  !> パラメータ：前処理情報保存の有無
  integer(kint), parameter :: monolis_prm_is_prec_stored = 13
  !> パラメータ：反復回数と残差履歴の表示
  integer(kint), parameter :: monolis_prm_show_iterlog = 14
  !> パラメータ：詳細な計算時間の表示
  integer(kint), parameter :: monolis_prm_show_time = 15
  !> パラメータ：ソルバ収束後のサマリの表示
  integer(kint), parameter :: monolis_prm_show_summary = 16
  !> パラメータ：計算時間の統計的処理結果の表示
  integer(kint), parameter :: monolis_prm_show_time_statistics = 17

  !> パラメータ：収束判定閾値
  integer(kint), parameter :: monolis_prm_tol = 1
  !> パラメータ：現在の残差
  integer(kint), parameter :: monolis_prm_cur_resid = 2
  !> パラメータ：ソルバの全計算時間
  integer(kint), parameter :: monolis_time_sol = 3
  !> パラメータ：前処理時間（生成時間）
  integer(kint), parameter :: monolis_time_prep = 4
  !> パラメータ：疎行列ベクトル積時間
  integer(kint), parameter :: monolis_time_spmv = 5
  !> パラメータ：ベクトル内積時間
  integer(kint), parameter :: monolis_time_dotp = 6
  !> パラメータ：前処理時間（適用時間）
  integer(kint), parameter :: monolis_time_prec = 7
  !> パラメータ：ベクトル内積の通信時間
  integer(kint), parameter :: monolis_time_comm_dotp = 8
  !> パラメータ：疎行列ベクトル積の通信時間
  integer(kint), parameter :: monolis_time_comm_spmv = 9

  !> パラメータ 構造体
  type monolis_prm
    !> 整数パラメータ
    integer(kint) :: monolis_prm_Iarray(monolis_prm_Iarray_size) = 0
    !> 実数パラメータ
    real(kdouble) :: monolis_prm_Rarray(monolis_prm_Rarray_size) = 0.0d0
  end type monolis_prm

contains

  !> パラメータ 構造体の初期化処理
  subroutine monolis_prm_initialize(monoPRM)
    implicit none
    !> パラメータ 構造体
    type(monolis_prm) :: monoPRM

    monoPRM%monolis_prm_Iarray(monolis_prm_method) = 1
    monoPRM%monolis_prm_Iarray(monolis_prm_precond) = 1
    monoPRM%monolis_prm_Iarray(monolis_prm_max_iter) = 1000
    monoPRM%monolis_prm_Iarray(monolis_prm_cur_iter) = 0
    monoPRM%monolis_prm_Iarray(monolis_prm_ierr) = -1
    monoPRM%monolis_prm_Iarray(monolis_prm_is_scaling) = monolis_I_false
    monoPRM%monolis_prm_Iarray(monolis_prm_is_reordering) = monolis_I_false
    monoPRM%monolis_prm_Iarray(monolis_prm_is_init_x) = monolis_I_true
    monoPRM%monolis_prm_Iarray(monolis_prm_is_sym_matrix) = monolis_I_false
    monoPRM%monolis_prm_Iarray(monolis_prm_is_debug) = monolis_I_false
    monoPRM%monolis_prm_Iarray(monolis_prm_is_measurement) = monolis_I_false
    monoPRM%monolis_prm_Iarray(monolis_prm_is_check_diag) = monolis_I_false
    monoPRM%monolis_prm_Iarray(monolis_prm_is_prec_stored) = monolis_I_false
    monoPRM%monolis_prm_Iarray(monolis_prm_show_iterlog) = monolis_I_true
    monoPRM%monolis_prm_Iarray(monolis_prm_show_time) = monolis_I_true
    monoPRM%monolis_prm_Iarray(monolis_prm_show_summary) = monolis_I_true
    monoPRM%monolis_prm_Iarray(monolis_prm_show_time_statistics) = monolis_I_false

    monoPRM%monolis_prm_Rarray(monolis_prm_tol) = 0.0d0
    monoPRM%monolis_prm_Rarray(monolis_prm_cur_resid) = 0.0d0
    monoPRM%monolis_prm_Rarray(monolis_time_sol) = 0.0d0
    monoPRM%monolis_prm_Rarray(monolis_time_prep) = 0.0d0
    monoPRM%monolis_prm_Rarray(monolis_time_spmv) = 0.0d0
    monoPRM%monolis_prm_Rarray(monolis_time_dotp) = 0.0d0
    monoPRM%monolis_prm_Rarray(monolis_time_prec) = 0.0d0
    monoPRM%monolis_prm_Rarray(monolis_time_comm_dotp) = 0.0d0
    monoPRM%monolis_prm_Rarray(monolis_time_comm_spmv) = 0.0d0
  end subroutine monolis_prm_initialize

  !> パラメータ 構造体の終了処理
  subroutine monolis_prm_finalize(monoPRM)
    implicit none
    !> パラメータ 構造体
    type(monolis_prm) :: monoPRM

    monoPRM%monolis_prm_Iarray = 0
    monoPRM%monolis_prm_Rarray = 0.0d0
  end subroutine monolis_prm_finalize
end module mod_monolis_def_solver