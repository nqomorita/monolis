!> ソルバパラメータの定義テスト
module mod_monolis_def_solver_test
  use mod_monolis

  implicit none

contains

  subroutine monolis_def_solver_test()
    implicit none

    call monolis_std_global_log_string("monolis_prm_initialize")
    call monolis_std_global_log_string("monolis_prm_finalize")

    call monolis_std_global_log_string("monolis_timer_initialize")
    call monolis_std_global_log_string("monolis_timer_finalize")
    call monolis_std_global_log_string("monolis_time_statistics")
  end subroutine monolis_def_solver_test

end module mod_monolis_def_solver_test