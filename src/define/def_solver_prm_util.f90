!> ソルバパラメータの設定関数群
!> comm
!# subroutine monolis_set_input_top_directory_name(monolis, param)
!# subroutine monolis_set_input_part_directory_name(monolis, param)
!# subroutine monolis_set_input_file_name(monolis, param)
!# subroutine monolis_set_communicator(monolis, comm)
!# subroutine monolis_set_my_rank(monolis, my_rank)
!# subroutine monolis_set_comm_size(monolis, comm_size)
!# subroutine monolis_set_n_internal_vertex(monolis, n_internal_vertex)
!# subroutine monolis_get_n_internal_vertex(monolis, n_internal_vertex)
!> solver
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

  !> 読込ファイルのトップディレクトリの設定
  subroutine monolis_set_input_top_directory_name(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    character(*) :: param
    monolis%PRM%com_top_dir_name = trim(param)
  end subroutine monolis_set_input_top_directory_name

  !> 読込ファイルの分割データディレクトリの設定
  subroutine monolis_set_input_part_directory_name(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    character(*) :: param
    monolis%PRM%com_part_dir_name = trim(param)
  end subroutine monolis_set_input_part_directory_name

  !> 読込ファイル名の設定
  subroutine monolis_set_input_file_name(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    character(*) :: param
    monolis%PRM%com_file_name = trim(param)
  end subroutine monolis_set_input_file_name

  !> @ingroup mpi
  !> MPI のグローバルコミュニケータを取得する関数
  function monolis_get_global_comm()
    implicit none
    integer(kint) :: monolis_get_global_comm
    monolis_get_global_comm = monolis_mpi_get_global_comm()
  end function monolis_get_global_comm

  !> @ingroup mpi
  !> MPI のグローバルランクサイズを取得する関数
  function monolis_get_global_comm_size()
    implicit none
    integer(kint) :: monolis_get_global_comm_size
    monolis_get_global_comm_size = monolis_mpi_get_global_comm_size()
  end function monolis_get_global_comm_size

  !> @ingroup mpi
  !> MPI のグローバルランクを取得する関数
  function monolis_get_global_my_rank()
    implicit none
    integer(kint) :: monolis_get_global_my_rank
    monolis_get_global_my_rank = monolis_mpi_get_global_my_rank()
  end function monolis_get_global_my_rank

  !> @ingroup mpi
  !> MPI のローカルコミュニケータのランクサイズを取得する関数
  function monolis_get_local_comm_size(monolis)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> [out] コミュニケータサイズ
    integer(kint) :: monolis_get_local_comm_size
    integer(kint) :: ierr
    monolis_get_local_comm_size = monolis_mpi_get_local_comm_size(monolis%COM%comm)
  end function monolis_get_local_comm_size

  !> @ingroup mpi
  !> MPI のローカルコミュニケータのランクサイズを取得する関数
  function monolis_get_local_my_rank(monolis)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> [out] MPI ランク番号
    integer(kint) :: monolis_get_local_my_rank
    monolis_get_local_my_rank = monolis_mpi_get_local_my_rank(monolis%COM%comm)
  end function monolis_get_local_my_rank

  !> @ingroup mpi
  !> MPI バリア関数（グローバルコミュニケータ）
  subroutine monolis_global_barrier()
    implicit none
    call monolis_mpi_global_barrier()
  end subroutine monolis_global_barrier

  !> @ingroup mpi
  !> MPI バリア関数（ローカルコミュニケータ）
  subroutine monolis_local_barrier(monolis)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    call monolis_mpi_local_barrier(monolis%COM%comm)
  end subroutine monolis_local_barrier

  !> @ingroup mpi
  !> MPI コミュニケータの分割
  subroutine monolis_split_communicator(monolis, group_id, comm_split)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> [in] コミュニケータのグループ id
    integer(kint) :: group_id
    !> [out] 分割後の MPI コミュニケータ
    integer(kint) :: comm_split
    call monolis_mpi_split_comm(monolis%COM%comm, group_id, comm_split)
  end subroutine monolis_split_communicator

  !> @ingroup com
  !> monolis 構造体に MPI コミュニケータを設定
  subroutine monolis_set_communicator(monolis, comm)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> MPI コミュニケータ
    integer(kint) :: comm
    call monolis_com_set_communicator(monolis%COM, comm)
  end subroutine monolis_set_communicator

  !> @ingroup com
  !> monolis 構造体に MPI ランク番号を設定
  subroutine monolis_set_my_rank(monolis, my_rank)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> MPI ランク番号
    integer(kint) :: my_rank
    call monolis_com_set_my_rank(monolis%COM, my_rank)
  end subroutine monolis_set_my_rank

  !> @ingroup com
  !> monolis 構造体に MPI コミュニケータサイズを設定
  subroutine monolis_set_comm_size(monolis, comm_size)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> MPI コミュニケータサイズ
    integer(kint) :: comm_size
    call monolis_com_set_comm_size(monolis%COM, comm_size)
  end subroutine monolis_set_comm_size

  !> @ingroup com
  !> monolis 構造体に内部領域に属する自由度数を設定
  subroutine monolis_set_n_internal_vertex(monolis, n_internal_vertex)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> 内部領域に属する自由度数
    integer(kint) :: n_internal_vertex
    call monolis_com_set_n_internal_vertex(monolis%COM, n_internal_vertex)
  end subroutine monolis_set_n_internal_vertex

  !> @ingroup com
  !> monolis 構造体に内部領域に属する自由度数を取得
  subroutine monolis_get_n_internal_vertex(monolis, n_internal_vertex)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> 内部領域に属する自由度数
    integer(kint) :: n_internal_vertex
    call monolis_com_get_n_internal_vertex(monolis%COM, n_internal_vertex)
  end subroutine monolis_get_n_internal_vertex

  !> @ingroup com
  !> monolis 構造体に内部領域に属する単一メッシュのリストを取得
  subroutine monolis_get_internal_simple_mesh_list(monolis, n_elem, n_base, elem, list)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> 要素数
    integer(kint) :: n_elem
    !> 要素の基底数
    integer(kint) :: n_base
    !> 要素コネクティビティ
    integer(kint) :: elem(:,:)
    !> 内部領域に属する自由度数
    integer(kint) :: list(:)
    integer(kint) :: i, in, j, id(n_base), n_internal_vertex, my_rank
    integer(kint), allocatable :: domain_id(:)

    call monolis_get_n_internal_vertex(monolis, n_internal_vertex)

    my_rank = monolis_get_local_my_rank(monolis)

    call monolis_alloc_I_1d(domain_id, monolis%MAT%NP)

    do i = 1, n_internal_vertex
      domain_id(i) = my_rank
    enddo

    call monolis_mpi_update_I(monolis%COM, 1, domain_id)

    list = 0
    do i = 1, n_elem
      do j = 1, n_base
        id(j) = domain_id(elem(j,i))
      enddo
      in = minval(id)
      if(in == my_rank) list(i) = 1
    enddo
  end subroutine monolis_get_internal_simple_mesh_list

  !> @ingroup com
  !> monolis 構造体に内部領域に属するコネクティビティのリストを取得
  subroutine monolis_get_internal_connectivity_list(monolis, n_elem, index, item, list)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> 要素数
    integer(kint) :: n_elem
    !> 要素コネクティビティの index 配列
    integer(kint) :: index(:)
    !> 要素コネクティビティの item 配列
    integer(kint) :: item(:)
    !> 内部領域に属する自由度数
    integer(kint) :: list(:)
    integer(kint) :: i, in, j, jS, jE, n_internal_vertex, my_rank
    integer(kint), allocatable :: id(:)
    integer(kint), allocatable :: domain_id(:)

    call monolis_get_n_internal_vertex(monolis, n_internal_vertex)

    my_rank = monolis_get_local_my_rank(monolis)

    call monolis_alloc_I_1d(domain_id, monolis%MAT%NP)

    do i = 1, n_internal_vertex
      domain_id(i) = my_rank
    enddo

    call monolis_mpi_update_I(monolis%COM, 1, domain_id)

    list = 0
    do i = 1, n_elem
      jS = index(i) + 1
      jE = index(i + 1)
      call monolis_alloc_I_1d(id, jE - jS + 1)
      do j = jS, jE
        id(j - jS + 1) = domain_id(item(j))
      enddo
      call monolis_dealloc_I_1d(id)
      in = minval(id)
      if(in == my_rank) list(i) = 1
    enddo
  end subroutine monolis_get_internal_connectivity_list

  !> ソルバの設定
  subroutine monolis_set_method(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    integer(kint) :: param
    monolis%PRM%Iarray(monolis_prm_I_method) = param
  end subroutine monolis_set_method

  !> 前処理の設定
  subroutine monolis_set_precond(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    integer(kint) :: param
    monolis%PRM%Iarray(monolis_prm_I_precond) = param
  end subroutine monolis_set_precond

  !> 最大反復回数の設定
  subroutine monolis_set_maxiter(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    integer(kint) :: param
    monolis%PRM%Iarray(monolis_prm_I_max_iter) = param
  end subroutine monolis_set_maxiter

  !> 現在の反復回数の取得
  subroutine monolis_get_converge_iter(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    integer(kint) :: param
    param = monolis%PRM%Iarray(monolis_prm_I_cur_iter)
  end subroutine monolis_get_converge_iter

  !> エラー番号の取得
  subroutine monolis_get_error_tag(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    integer(kint) :: param
    param = monolis%PRM%Iarray(monolis_prm_I_ierr)
  end subroutine monolis_get_error_tag

  !> スケーリングの有無の設定
  !subroutine monolis_set_scaling(monolis, param)
  !  implicit none
  !  !> monolis 構造体
  !  type(monolis_structure) :: monolis
  !  !> パラメータ
  !  logical :: param
  !  integer(kint) :: iparam

  !  monolis%PRM%Iarray(monolis_prm_I_is_scaling) = monolis_conv_L2I(param)
  !end subroutine monolis_set_scaling

  !> リオーダリングの有無の設定
  !subroutine monolis_set_reordering(monolis, param)
  !  implicit none
  !  !> monolis 構造体
  !  type(monolis_structure) :: monolis
  !  !> パラメータ
  !  logical :: param
  !  integer(kint) :: iparam

  !  monolis%PRM%Iarray(monolis_prm_I_is_reordering) = monolis_conv_L2I(param)
  !end subroutine monolis_set_reordering

  !> 解ベクトル初期化の有無の設定
  subroutine monolis_set_init_x(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    logical :: param
    integer(kint) :: iparam

    monolis%PRM%Iarray(monolis_prm_I_is_init_x) = monolis_conv_L2I(param)
  end subroutine monolis_set_init_x

  !> 対称行列向け処理の有無の設定
  subroutine monolis_set_sym_matrix(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    logical :: param
    integer(kint) :: iparam

    monolis%PRM%Iarray(monolis_prm_I_is_sym_matrix) = monolis_conv_L2I(param)
  end subroutine monolis_set_sym_matrix

  !> デバッグ出力の有無の設定
  subroutine monolis_set_debug(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    logical :: param
    integer(kint) :: iparam

    monolis%PRM%Iarray(monolis_prm_I_is_debug) = monolis_conv_L2I(param)
  end subroutine monolis_set_debug

  !> 詳細な計算時間測定の有無の設定
  subroutine monolis_set_performance_measurement(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    logical :: param
    integer(kint) :: iparam

    monolis%PRM%Iarray(monolis_prm_I_is_measurement) = monolis_conv_L2I(param)
  end subroutine monolis_set_performance_measurement

  !> 行列対角成分確認の有無の設定
  subroutine monolis_set_check_diag(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    logical :: param
    integer(kint) :: iparam

    monolis%PRM%Iarray(monolis_prm_I_is_check_diag) = monolis_conv_L2I(param)
  end subroutine monolis_set_check_diag

  !> 前処理情報保存の有無の設定
  subroutine monolis_set_prec_stored(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    logical :: param
    integer(kint) :: iparam

    monolis%PRM%Iarray(monolis_prm_I_is_prec_stored) = monolis_conv_L2I(param)
  end subroutine monolis_set_prec_stored

  !> 反復回数と残差履歴の表示の設定
  subroutine monolis_show_iterlog(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    logical :: param
    integer(kint) :: iparam

    monolis%PRM%Iarray(monolis_prm_I_show_iterlog) = monolis_conv_L2I(param)
  end subroutine monolis_show_iterlog

  !> 詳細な計算時間の表示の設定
  subroutine monolis_show_timelog(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    logical :: param
    integer(kint) :: iparam

    monolis%PRM%Iarray(monolis_prm_I_show_time) = monolis_conv_L2I(param)
  end subroutine monolis_show_timelog

  !> ソルバ収束後のサマリの表示の設定
  subroutine monolis_show_summary(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    logical :: param
    integer(kint) :: iparam

    monolis%PRM%Iarray(monolis_prm_I_show_summary) = monolis_conv_L2I(param)
  end subroutine monolis_show_summary

  !> 計算時間の統計的処理結果の表示の設定
  subroutine monolis_show_timelog_statistics(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    logical :: param
    integer(kint) :: iparam

    monolis%PRM%Iarray(monolis_prm_I_show_time_statistics) = monolis_conv_L2I(param)
  end subroutine monolis_show_timelog_statistics

  !> 収束判定閾値の設定
  subroutine monolis_set_tolerance(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    real(kdouble) :: param
    monolis%PRM%Rarray(monolis_prm_R_tol) = param
  end subroutine monolis_set_tolerance

  !> 現在の残差の取得
  subroutine monolis_get_converge_residual(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    real(kdouble) :: param
    param = monolis%PRM%Rarray(monolis_prm_R_cur_resid)
  end subroutine monolis_get_converge_residual

  !> ソルバの全計算時間の取得
  subroutine monolis_get_time_solver(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    real(kdouble) :: param
    param = monolis%PRM%Rarray(monolis_R_time_sol)
  end subroutine monolis_get_time_solver

  !> 前処理時間（生成時間）の取得
  subroutine monolis_get_time_preparing(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    real(kdouble) :: param
    param = monolis%PRM%Rarray(monolis_R_time_prep)
  end subroutine monolis_get_time_preparing

  !> 疎行列ベクトル積時間の取得
  subroutine monolis_get_time_spmv(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    real(kdouble) :: param
    param = monolis%PRM%Rarray(monolis_R_time_spmv)
  end subroutine monolis_get_time_spmv

  !> ベクトル内積時間の取得
  subroutine monolis_get_time_inner_product(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    real(kdouble) :: param
    param = monolis%PRM%Rarray(monolis_R_time_dotp)
  end subroutine monolis_get_time_inner_product

  !> 前処理時間（適用時間）の取得
  subroutine monolis_get_time_precondition(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    real(kdouble) :: param
    param = monolis%PRM%Rarray(monolis_R_time_prec)
  end subroutine monolis_get_time_precondition

  !> ベクトル内積の通信時間の取得
  subroutine monolis_get_time_comm_inner_product(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    real(kdouble) :: param
    param = monolis%PRM%Rarray(monolis_R_time_comm_dotp)
  end subroutine monolis_get_time_comm_inner_product

  !> 疎行列ベクトル積の通信時間の取得
  subroutine monolis_get_time_comm_spmv(monolis, param)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> パラメータ
    real(kdouble) :: param
    param = monolis%PRM%Rarray(monolis_R_time_comm_spmv)
  end subroutine monolis_get_time_comm_spmv
end module mod_monolis_def_solver_util
