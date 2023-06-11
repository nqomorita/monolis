!> monolis 構造体の定義テスト
module mod_monolis_def_struc_test
  use mod_monolis

  implicit none

contains

  subroutine monolis_def_struc_test()
    implicit none

    call monolis_std_global_log_string("monolis_global_initialize")
    call monolis_std_global_log_string("monolis_global_finalize")

    call monolis_std_global_log_string("monolis_initialize")

    call monolis_std_global_log_string("monolis_finalize")

    call monolis_std_global_log_string("monolis_com_input_comm_table")
  end subroutine monolis_def_struc_test

end module mod_monolis_def_struc_test
