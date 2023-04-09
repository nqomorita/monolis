!> 線形ソルバテストモジュール
module mod_monolis_solve_test
  use mod_monolis

  implicit none

contains

  subroutine monolis_solve_test()
    implicit none

    call monolis_std_global_log_string("monolis_solve_R")
    call monolis_std_global_log_string("monolis_solve_C")
    call monolis_std_global_log_string("monolis_solve_main_R")
    call monolis_std_global_log_string("monolis_solve_main_C")
    call monolis_std_global_log_string("monolis_solver_select_R")
    call monolis_std_global_log_string("monolis_solver_select_C")
  end subroutine monolis_solve_test

end module mod_monolis_solve_test
