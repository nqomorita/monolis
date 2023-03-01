!> 対角スケーリング前処理（nxn ブロック）テストモジュール
module mod_monolis_precond_diag_nn_test
  use mod_monolis

  implicit none

contains

  subroutine monolis_precond_diag_nn_test
    implicit none

    !call monolis_precond_diag_nn_setup(monoMAT, monoPREC)
    !call monolis_precond_diag_nn_apply(monoMAT, monoPREC, X, Y)
    !call monolis_precond_diag_nn_clear(monoPREC)
  end subroutine monolis_precond_diag_nn_test

end module mod_monolis_precond_diag_nn_test
