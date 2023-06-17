!> ソルバパラメータの定義
module mod_monolis_def_solver
  use mod_monolis_utils
  implicit none

  !> パラメータ：CG 法
  integer(kint), parameter :: monolis_iter_CG       = 1
  !> パラメータ：GropCG 法
  integer(kint), parameter :: monolis_iter_GropCG   = 2
  !> パラメータ：Pipelined CG 法
  integer(kint), parameter :: monolis_iter_PipeCG   = 3
  !> パラメータ：Pipelined CR 法
  integer(kint), parameter :: monolis_iter_PipeCR   = 4
  !> パラメータ：BiCGSTAB 法
  integer(kint), parameter :: monolis_iter_BiCGSTAB = 5
  !> パラメータ：Pipelined BiCGSTAB 法
  integer(kint), parameter :: monolis_iter_PipeBiCGSTAB = 6
  !> パラメータ：BiCGSTAB 法（前処理なし）
  integer(kint), parameter :: monolis_iter_BiCGSTAB_noprec = 7
  !> パラメータ：Communication Avoiding BiCGSTAB 法（前処理なし）
  !integer(kint), parameter :: monolis_iter_CABiCGSTAB_noprec = 8
  !> パラメータ：Pipelined BiCGSTAB 法（前処理なし）
  integer(kint), parameter :: monolis_iter_PipeBiCGSTAB_noprec = 8
  !> パラメータ：GMRES 法
  !integer(kint), parameter :: monolis_iter_GMRES    = 10
  !> パラメータ：COCG 法
  integer(kint), parameter :: monolis_iter_COCG       = 9

  !> パラメータ：前処理なし
  integer(kint), parameter :: monolis_prec_NONE   = 0
  !> パラメータ：対角スケーリング前処理
  integer(kint), parameter :: monolis_prec_DIAG   = 1
  !> パラメータ：ブロック ILU 前処理
  !integer(kint), parameter :: monolis_prec_ILU    = 2
  !> パラメータ：ブロック Jacobi 前処理
  !integer(kint), parameter :: monolis_prec_JACOBI = 3
  !> パラメータ：ブロック SOR 前処理
  integer(kint), parameter :: monolis_prec_SOR    = 2
  !> パラメータ：SPIKE 前処理
  !integer(kint), parameter :: monolis_prec_SPIKE  = 5
  !> パラメータ：LU 分解
  !integer(kint), parameter :: monolis_prec_LU     = 6
  !> パラメータ：LU 分解（MUMPS）
  integer(kint), parameter :: monolis_prec_MUMPS  = 3
  !> パラメータ：ROM 前処理
  !integer(kint), parameter :: monolis_prec_ROM    = 8
  !> パラメータ：LU 分解（MF 法）
  !integer(kint), parameter :: monolis_prec_MF     = 9
  !> パラメータ：ブロック LU 分解（MUMPS）
  integer(kint), parameter :: monolis_prec_MUMPS_LOCAL = 4

  character*24, dimension(9) :: monolis_str_iter = (/&
  & "CG                 ", &
  & "GropCG             ", &
  & "PipeCG             ", &
  & "PipeCR             ", &
  & "BiCGSTAB           ", &
  & "PipeBiCGSTAB       ", &
  & "BiCGSTAB_noprec    ", &
  !& "CABiCGSTAB_noprec  ", &
  & "PipeBiCGSTAB_noprec", &
  !& "GMRES              ", &
  & "COCG               "/)

  character*24, dimension(0:4)  :: monolis_str_prec = (/&
  & "None  ", &
  & "Diag  ", &
!  & "ILU   ", &
!  & "Jacobi", &
  & "SOR   ", &
!  & "SPIKE ", &
!  & "LU    ", &
  & "MUMPS ", &
!  & "ROM   ", &
!  & "MF    ", &
  & "MUMPSL"/)

  !> 整数パラメータのサイズ
  integer(kint), parameter :: monolis_prm_Iarray_size = 100

  !> 実数パラメータのサイズ
  integer(kint), parameter :: monolis_prm_Rarray_size = 100

  !> パラメータ：ソルバ
  integer(kint), parameter :: monolis_prm_I_method = 1
  !> パラメータ：前処理
  integer(kint), parameter :: monolis_prm_I_precond = 2
  !> パラメータ：最大反復回数
  integer(kint), parameter :: monolis_prm_I_max_iter = 3
  !> パラメータ：現在の反復回数
  integer(kint), parameter :: monolis_prm_I_cur_iter = 4
  !> パラメータ：エラー番号
  integer(kint), parameter :: monolis_prm_I_ierr = 5
  !> パラメータ：スケーリングの有無
  !integer(kint), parameter :: monolis_prm_I_is_scaling = 6
  !> パラメータ：リオーダリングの有無
  !integer(kint), parameter :: monolis_prm_I_is_reordering = 7
  !> パラメータ：解ベクトル初期化の有無
  integer(kint), parameter :: monolis_prm_I_is_init_x = 8
  !> パラメータ：対称行列向け処理の有無
  integer(kint), parameter :: monolis_prm_I_is_sym_matrix = 9
  !> パラメータ：デバッグ出力の有無
  integer(kint), parameter :: monolis_prm_I_is_debug = 10
  !> パラメータ：詳細な計算時間測定の有無
  integer(kint), parameter :: monolis_prm_I_is_measurement = 11
  !> パラメータ：行列対角成分確認の有無
  integer(kint), parameter :: monolis_prm_I_is_check_diag = 12
  !> パラメータ：前処理情報保存の有無
  integer(kint), parameter :: monolis_prm_I_is_prec_stored = 13
  !> パラメータ：反復回数と残差履歴の表示
  integer(kint), parameter :: monolis_prm_I_show_iterlog = 14
  !> パラメータ：詳細な計算時間の表示
  integer(kint), parameter :: monolis_prm_I_show_time = 15
  !> パラメータ：ソルバ収束後のサマリの表示
  integer(kint), parameter :: monolis_prm_I_show_summary = 16
  !> パラメータ：計算時間の統計的処理結果の表示
  integer(kint), parameter :: monolis_prm_I_show_time_statistics = 17

  !> パラメータ：収束判定閾値
  integer(kint), parameter :: monolis_prm_R_tol = 1
  !> パラメータ：現在の残差
  integer(kint), parameter :: monolis_prm_R_cur_resid = 2
  !> パラメータ：ソルバの全計算時間
  integer(kint), parameter :: monolis_R_time_sol = 3
  !> パラメータ：前処理時間（生成時間）
  integer(kint), parameter :: monolis_R_time_prep = 4
  !> パラメータ：疎行列ベクトル積時間
  integer(kint), parameter :: monolis_R_time_spmv = 5
  !> パラメータ：ベクトル内積時間
  integer(kint), parameter :: monolis_R_time_dotp = 6
  !> パラメータ：前処理時間（適用時間）
  integer(kint), parameter :: monolis_R_time_prec = 7
  !> パラメータ：ベクトル内積の通信時間
  integer(kint), parameter :: monolis_R_time_comm_dotp = 8
  !> パラメータ：疎行列ベクトル積の通信時間
  integer(kint), parameter :: monolis_R_time_comm_spmv = 9

  !> パラメータ構造体
  type monolis_prm
    !> 整数パラメータ
    integer(kint) :: Iarray(monolis_prm_Iarray_size) = 0
    !> 実数パラメータ
    real(kdouble) :: Rarray(monolis_prm_Rarray_size) = 0.0d0
    !> 実数パラメータ
    character(monolis_charlen) :: com_top_dir_name
    !> 実数パラメータ
    character(monolis_charlen) :: com_part_dir_name
    !> 実数パラメータ
    character(monolis_charlen) :: com_file_name
  end type monolis_prm

contains

  !> @ingroup def_init
  !> パラメータ 構造体の初期化処理
  subroutine monolis_prm_initialize(monoPRM)
    implicit none
    !> [out] パラメータ構造体
    type(monolis_prm), intent(out) :: monoPRM

    monoPRM%com_top_dir_name = "."
    monoPRM%com_part_dir_name = "parted.0"
    monoPRM%com_file_name = "node.dat"

    monoPRM%Iarray(monolis_prm_I_method) = 1
    monoPRM%Iarray(monolis_prm_I_precond) = 1
    monoPRM%Iarray(monolis_prm_I_max_iter) = 1000
    monoPRM%Iarray(monolis_prm_I_cur_iter) = 0
    monoPRM%Iarray(monolis_prm_I_ierr) = -1
    !monoPRM%Iarray(monolis_prm_I_is_scaling) = monolis_I_false
    !monoPRM%Iarray(monolis_prm_I_is_reordering) = monolis_I_false
    monoPRM%Iarray(monolis_prm_I_is_init_x) = monolis_I_true
    monoPRM%Iarray(monolis_prm_I_is_sym_matrix) = monolis_I_false
    monoPRM%Iarray(monolis_prm_I_is_debug) = monolis_I_false
    monoPRM%Iarray(monolis_prm_I_is_measurement) = monolis_I_false
    monoPRM%Iarray(monolis_prm_I_is_check_diag) = monolis_I_false
    monoPRM%Iarray(monolis_prm_I_is_prec_stored) = monolis_I_false
    monoPRM%Iarray(monolis_prm_I_show_iterlog) = monolis_I_true
    monoPRM%Iarray(monolis_prm_I_show_time) = monolis_I_true
    monoPRM%Iarray(monolis_prm_I_show_summary) = monolis_I_true
    monoPRM%Iarray(monolis_prm_I_show_time_statistics) = monolis_I_false

    monoPRM%Rarray(monolis_prm_R_tol) = 1.0d-8
    monoPRM%Rarray(monolis_prm_R_cur_resid) = 0.0d0
    monoPRM%Rarray(monolis_R_time_sol) = 0.0d0
    monoPRM%Rarray(monolis_R_time_prep) = 0.0d0
    monoPRM%Rarray(monolis_R_time_spmv) = 0.0d0
    monoPRM%Rarray(monolis_R_time_dotp) = 0.0d0
    monoPRM%Rarray(monolis_R_time_prec) = 0.0d0
    monoPRM%Rarray(monolis_R_time_comm_dotp) = 0.0d0
    monoPRM%Rarray(monolis_R_time_comm_spmv) = 0.0d0
  end subroutine monolis_prm_initialize

  !> @ingroup def_init
  !> パラメータ構造体の終了処理
  subroutine monolis_prm_finalize(monoPRM)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM

    monoPRM%Iarray = 0
    monoPRM%Rarray = 0.0d0
  end subroutine monolis_prm_finalize

  !> @ingroup def_init
  !> 時間計測機能の初期化処理
  subroutine monolis_timer_initialize(monoPRM)
    implicit none
    !> [inout] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM

    monoPRM%Rarray(monolis_R_time_sol) = 0.0d0
    monoPRM%Rarray(monolis_R_time_prep) = 0.0d0
    monoPRM%Rarray(monolis_R_time_spmv) = 0.0d0
    monoPRM%Rarray(monolis_R_time_dotp) = 0.0d0
    monoPRM%Rarray(monolis_R_time_prec) = 0.0d0
    monoPRM%Rarray(monolis_R_time_comm_dotp) = 0.0d0
    monoPRM%Rarray(monolis_R_time_comm_spmv) = 0.0d0

    monoPRM%Rarray(monolis_R_time_sol) = monolis_get_time()
  end subroutine monolis_timer_initialize

  !> @ingroup def_init
  !> 時間計測機能の終了処理
  subroutine monolis_timer_finalize(monoPRM, monoCOM)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    real(kdouble) :: t1, time(6), t_max, t_min, t_avg, t_sd
    logical :: is_output

    call monolis_std_debug_log_header("monolis_timer_finalize")

    call monolis_mpi_local_barrier(monoCOM%comm)

    t1 = monolis_get_time()

    monoPRM%Rarray(monolis_R_time_sol) = t1 - monoPRM%Rarray(monolis_R_time_sol)

    if(monoPRM%Iarray(monolis_prm_I_show_summary) == monolis_I_true .and. monoCOM%my_rank == 0)then
      write(*,"(a,i10)")"** MONOLIS converge iter:", monoPRM%Iarray(monolis_prm_I_cur_iter)
      write(*,"(a,1p4e10.3)")"** MONOLIS rel. residual:", monoPRM%Rarray(monolis_prm_R_cur_resid)
    endif

    if((monoPRM%Iarray(monolis_prm_I_show_summary) == monolis_I_true .or. &
      & monoPRM%Iarray(monolis_prm_I_show_time)    == monolis_I_true) .and. monoCOM%my_rank == 0)then
      write(*,"(a,1p4e10.3)")"** MONOLIS solution time:", monoPRM%Rarray(monolis_R_time_sol)
    endif

    if(monoPRM%Iarray(monolis_prm_I_show_time) == monolis_I_true)then
      time(1) = monoPRM%Rarray(monolis_R_time_prep)
      time(2) = monoPRM%Rarray(monolis_R_time_spmv)
      time(3) = monoPRM%Rarray(monolis_R_time_dotp)
      time(4) = monoPRM%Rarray(monolis_R_time_prec)
      time(5) = monoPRM%Rarray(monolis_R_time_comm_dotp)
      time(6) = monoPRM%Rarray(monolis_R_time_comm_spmv)

      call monolis_allreduce_R(6, time, monolis_mpi_sum, monoCOM%comm)

      time = time/dble(monolis_mpi_get_global_comm_size())

      if(monoCOM%my_rank == 0)then
        write(*,"(a,1p4e10.3)")" - solution/prepost time:", time(1)
        write(*,"(a,1p4e10.3)")" - solution/SpMV    time:", time(2)
        write(*,"(a,1p4e10.3)")" - solution/inner p time:", time(3)
        write(*,"(a,1p4e10.3)")" - solution/precond time:", time(4)
        write(*,"(a,1p4e10.3)")" - (comm time/inner p)  :", time(5)
        write(*,"(a,1p4e10.3)")" - (comm time/spmv)     :", time(6)
      endif
    endif

    if(monoPRM%Iarray(monolis_prm_I_show_time_statistics) == monolis_I_true)then
      if(monoCOM%my_rank == 0) write(*,"(a,1p4e10.3)")"** MONOLIS solution time statistics"
      if(monoCOM%my_rank == 0) write(*,"(a,1p4e10.3)")"                          max       min       average   std"

      time(1) = monoPRM%Rarray(monolis_R_time_prep)
      call monolis_time_statistics (time(1), t_max, t_min, t_avg, t_sd, monoCOM%comm)
      if(monoCOM%my_rank == 0) write(*,"(a,1p4e10.3)")" - solution/prepost time:", t_max, t_min, t_avg, t_sd

      time(2) = monoPRM%Rarray(monolis_R_time_spmv)
      call monolis_time_statistics (time(2), t_max, t_min, t_avg, t_sd, monoCOM%comm)
      if(monoCOM%my_rank == 0) write(*,"(a,1p4e10.3)")" - solution/SpMV    time:", t_max, t_min, t_avg, t_sd

      time(3) = monoPRM%Rarray(monolis_R_time_dotp)
      call monolis_time_statistics (time(3), t_max, t_min, t_avg, t_sd, monoCOM%comm)
      if(monoCOM%my_rank == 0) write(*,"(a,1p4e10.3)")" - solution/inner p time:", t_max, t_min, t_avg, t_sd

      time(4) = monoPRM%Rarray(monolis_R_time_prec)
      call monolis_time_statistics (time(4), t_max, t_min, t_avg, t_sd, monoCOM%comm)
      if(monoCOM%my_rank == 0) write(*,"(a,1p4e10.3)")" - solution/precond time:", t_max, t_min, t_avg, t_sd

      time(5) = monoPRM%Rarray(monolis_R_time_comm_dotp)
      call monolis_time_statistics (time(5), t_max, t_min, t_avg, t_sd, monoCOM%comm)
      if(monoCOM%my_rank == 0) write(*,"(a,1p4e10.3)")" - (comm time/inner p)  :", t_max, t_min, t_avg, t_sd

      time(6) = monoPRM%Rarray(monolis_R_time_comm_spmv)
      call monolis_time_statistics (time(6), t_max, t_min, t_avg, t_sd, monoCOM%comm)
      if(monoCOM%my_rank == 0) write(*,"(a,1p4e10.3)")" - (comm time/spmv)     :", t_max, t_min, t_avg, t_sd
    endif

    time(1) = monoPRM%Rarray(monolis_R_time_prep)
    time(2) = monoPRM%Rarray(monolis_R_time_spmv)
    time(3) = monoPRM%Rarray(monolis_R_time_dotp)
    time(4) = monoPRM%Rarray(monolis_R_time_prec)
    time(5) = monoPRM%Rarray(monolis_R_time_comm_dotp)
    time(6) = monoPRM%Rarray(monolis_R_time_comm_spmv)

    call monolis_allreduce_R(6, time, monolis_mpi_sum, monoCOM%comm)
    time = time/dble(monolis_mpi_get_global_comm_size())

    monoPRM%Rarray(monolis_R_time_prep) = time(1)
    monoPRM%Rarray(monolis_R_time_spmv) = time(2)
    monoPRM%Rarray(monolis_R_time_dotp) = time(3)
    monoPRM%Rarray(monolis_R_time_prec) = time(4)
    monoPRM%Rarray(monolis_R_time_comm_dotp) = time(5)
    monoPRM%Rarray(monolis_R_time_comm_spmv) = time(6)
  end subroutine monolis_timer_finalize

  !> @ingroup def_init
  !> 時間計測機能の統計処理
  subroutine monolis_time_statistics(time, t_max, t_min, t_avg, t_sd, comm)
    implicit none
    !> [in] 時間
    real(kdouble), intent(in) :: time
    !> [out] 最大値
    real(kdouble), intent(out) :: t_max
    !> [out] 最小値
    real(kdouble), intent(out) :: t_min
    !> [out] 平均値
    real(kdouble), intent(out) :: t_avg
    !> [out] 標準偏差
    real(kdouble), intent(out) :: t_sd
    !> [in] MPI コミュニケータ
    integer(kint), intent(in) :: comm
    real(kdouble) :: tmp
    integer(kint) :: np

    np = monolis_mpi_get_global_comm_size()

    t_max = time
    call monolis_allreduce_R1(t_max, monolis_mpi_max, comm)

    t_min = time
    call monolis_allreduce_R1(t_min, monolis_mpi_min, comm)

    t_avg = time
    call monolis_allreduce_R1(t_avg, monolis_mpi_sum, comm)
    t_avg = t_avg / np

    tmp = time*time
    call monolis_allreduce_R1(tmp, monolis_mpi_sum, comm)
    tmp = tmp / np

    if(np == 1)then
      t_sd = 0.0d0
    else
      t_sd = dsqrt(tmp - t_avg*t_avg)
    endif
  end subroutine monolis_time_statistics
end module mod_monolis_def_solver
