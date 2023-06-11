!> SOR 前処理関連テストモジュール
module mod_monolis_precond_sor_test
  use mod_monolis
  use mod_monolis_precond_sor_33_test
  use mod_monolis_precond_sor_nn_test

  implicit none

contains

  subroutine monolis_precond_sor_test()
    implicit none

    call monolis_std_log_string("monolis_precond_sor_setup_R")
    call monolis_std_log_string("monolis_precond_sor_apply_R")
    call monolis_std_log_string("monolis_precond_sor_clear_R")

    call monolis_std_log_string("monolis_precond_sor_setup_C")
    call monolis_std_log_string("monolis_precond_sor_apply_C")
    call monolis_std_log_string("monolis_precond_sor_clear_C")

    call  monolis_precond_sor_33_test()
    call  monolis_precond_sor_nn_test()
  end subroutine monolis_precond_sor_test
end module mod_monolis_precond_sor_test
