!> 対角スケーリング前処理関連テストモジュール
module mod_monolis_precond_diag_test
  use mod_monolis

  implicit none

contains

  subroutine monolis_precond_diag_test
    implicit none
    !call monolis_precond_diag_setup(monoPRM, monoCOM, monoMAT, monoPREC)
    !call monolis_precond_diag_apply(monoPRM, monoCOM, monoMAT, monoPREC, X, Y)
    !call monolis_precond_diag_clear(monoPRM, monoCOM, monoMAT, monoPREC)
  end subroutine monolis_precond_diag_test

end module mod_monolis_precond_diag_test
