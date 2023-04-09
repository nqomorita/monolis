!> 前処理関連モジュール
module mod_monolis_precond_test
  use mod_monolis
  use mod_monolis_precond_diag_test
  use mod_monolis_precond_sor_test

  implicit none

contains

  subroutine monolis_precond_test()
    implicit none

    call monolis_std_log_string("monolis_precond_setup")
    call monolis_std_log_string("monolis_precond_clear")

    call monolis_std_log_string("monolis_precond_setup_R")
    call monolis_std_log_string("monolis_precond_apply_R")
    call monolis_std_log_string("monolis_precond_clear_R")

    call monolis_std_log_string("monolis_precond_setup_C")
    call monolis_std_log_string("monolis_precond_apply_C")
    call monolis_std_log_string("monolis_precond_clear_C")

    call monolis_precond_diag_test()
    call monolis_precond_sor_test()
  end subroutine monolis_precond_test

end module mod_monolis_precond_test
