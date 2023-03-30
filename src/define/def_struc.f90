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
    !> 前処理構造体
    type(monolis_mat) :: PREC
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
  subroutine monolis_initialize(monolis)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis

    call monolis_prm_initialize(monolis%PRM)

    call monolis_com_initialize(monolis%COM)

    call monolis_mat_initialize(monolis%MAT)

    call monolis_mat_initialize(monolis%PREC)

    call monolis_com_input_comm_table(monolis%COM, &
      & monolis%PRM%com_top_dir_name, monolis%PRM%com_part_dir_name, monolis%PRM%com_file_name)
  end subroutine monolis_initialize

  !> monoils 構造体の終了処理
  subroutine monolis_finalize(monolis)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis

    call monolis_prm_finalize(monolis%PRM)
    call monolis_com_finalize(monolis%COM)
    call monolis_mat_finalize(monolis%MAT)
    call monolis_mat_finalize(monolis%PREC)
  end subroutine monolis_finalize

  subroutine monolis_com_input_comm_table(monoCOM, top_dir_name, part_dir_name, file_name)
    implicit none
    type(monolis_com) :: monoCOM
    integer(kint) :: i, j, nitem
    character(*) :: top_dir_name
    character(*) :: part_dir_name
    character(*) :: file_name
    character(monolis_charlen) :: header
    character(monolis_charlen) :: fname

    if(monoCOM%comm_size <= 1)then
      monoCOM%comm_size = 1
      return
    endif

    header = trim(top_dir_name)//"/"//trim(part_dir_name)

    fname = monolis_get_output_file_name_by_domain_id(trim(header), trim(file_name)//".recv", monoCOM%my_rank)

    call monolis_input_recv_com_table(fname, monoCOM)

    fname = monolis_get_output_file_name_by_domain_id(trim(header), trim(file_name)//".send", monoCOM%my_rank)

    call monolis_input_send_com_table(fname, monoCOM)

    fname = monolis_get_output_file_name_by_domain_id(trim(header), trim(file_name)//".n_internal", monoCOM%my_rank)

    call monolis_input_internal_vertex_number(fname, monoCOM%n_internal_vertex)
  end subroutine monolis_com_input_comm_table
end module mod_monolis_def_struc
