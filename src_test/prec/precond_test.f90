!> 前処理関連モジュール
module mod_monolis_precond_test
  use mod_monolis

  implicit none

contains

  subroutine monolis_precond_test
    implicit none

    !call monolis_precond_setup_R(monoPRM, monoCOM, monoMAT, monoPREC)
    !call monolis_precond_apply_R(monoPRM, monoCOM, monoMAT, monoPREC, X, Y)
    !call monolis_precond_clear_R(monoPRM, monoCOM, monoMAT, monoPREC)
  end subroutine monolis_precond_test

end module mod_monolis_precond_test
