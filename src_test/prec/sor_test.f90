!> SOR 前処理関連テストモジュール
module mod_monolis_precond_sor_test
  use mod_monolis

  implicit none

contains

  subroutine monolis_precond_sor_test
    implicit none

    !call monolis_precond_sor_setup(monoPRM, monoCOM, monoMAT, monoPREC)
    !call monolis_precond_sor_apply(monoPRM, monoCOM, monoMAT, monoPREC, X, Y)
    !call monolis_precond_sor_clear(monoPRM, monoCOM, monoMAT, monoPREC)
  end subroutine monolis_precond_sor_test
end module mod_monolis_precond_sor_test
