!> 対角スケーリング前処理関連テストモジュール
module mod_monolis_precond_diag_test
  use mod_monolis
  use mod_monolis_precond_diag_33_test
  use mod_monolis_precond_diag_nn_test

  implicit none

contains

  subroutine monolis_precond_diag_test()
    implicit none

    call monolis_std_log_string("monolis_precond_diag_setup_R")
    call monolis_std_log_string("monolis_precond_diag_apply_R")
    call monolis_std_log_string("monolis_precond_diag_clear_R")

    call monolis_std_log_string("monolis_precond_diag_setup_C")
    call monolis_std_log_string("monolis_precond_diag_apply_C")
    call monolis_std_log_string("monolis_precond_diag_clear_C")

    call  monolis_precond_diag_33_test()
    call  monolis_precond_diag_nn_test()
  end subroutine monolis_precond_diag_test

end module mod_monolis_precond_diag_test
