!> SOR 前処理（3x3 ブロック）テストモジュール
module mod_monolis_precond_sor_33_test
  use mod_monolis

  implicit none

contains

  subroutine monolis_precond_sor_33_test
    implicit none

    !call monolis_precond_sor_33_setup_R(monoMAT, monoPREC)
    !call monolis_precond_sor_33_apply_R(monoMAT, monoPREC, X, Y)
    !call monolis_precond_sor_33_clear_R(monoPREC)
  end subroutine monolis_precond_sor_33_test

end module mod_monolis_precond_sor_33_test
