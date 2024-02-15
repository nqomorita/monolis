!> SOR 前処理関連テストモジュール
module mod_monolis_precond_MUMPS_test
  use mod_monolis

  implicit none

contains

  subroutine monolis_precond_MUMPS_test()
    implicit none

    call monolis_std_log_string("monolis_precond_MUMPS_setup_R")
    call monolis_std_log_string("monolis_precond_MUMPS_apply_R")
    call monolis_std_log_string("monolis_precond_MUMPS_clear_R")

  end subroutine monolis_precond_MUMPS_test
end module mod_monolis_precond_MUMPS_test
