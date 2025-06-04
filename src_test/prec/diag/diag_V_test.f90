!> 対角スケーリング前処理（nxn ブロック）テストモジュール
module mod_monolis_precond_diag_V_test
  use mod_monolis

  implicit none

contains

  subroutine monolis_precond_diag_V_test
    implicit none

    call monolis_std_log_string("monolis_precond_diag_V_setup_R")
    call monolis_std_log_string("monolis_precond_diag_V_apply_R")
    call monolis_std_log_string("monolis_precond_diag_V_clear_R")

    call monolis_std_log_string("monolis_precond_diag_V_setup_C")
    call monolis_std_log_string("monolis_precond_diag_V_apply_C")
    call monolis_std_log_string("monolis_precond_diag_V_clear_C")
  end subroutine monolis_precond_diag_V_test

end module mod_monolis_precond_diag_V_test
