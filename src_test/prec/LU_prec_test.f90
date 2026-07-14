!> LU 前処理関連テストモジュール
module mod_monolis_precond_LU_test
  use mod_monolis

  implicit none

contains

  subroutine monolis_precond_LU_test()
    implicit none

    call monolis_std_global_log_string("monolis_precond_LU_setup_R")
    call monolis_std_global_log_string("monolis_precond_LU_apply_R")
    call monolis_std_global_log_string("monolis_precond_LU_clear_R")

    call monolis_std_global_log_string("monolis_precond_LU_setup_C")
    call monolis_std_global_log_string("monolis_precond_LU_apply_C")
    call monolis_std_global_log_string("monolis_precond_LU_clear_C")
  end subroutine monolis_precond_LU_test
end module mod_monolis_precond_LU_test
