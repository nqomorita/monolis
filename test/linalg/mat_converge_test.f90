!> 収束判定テストモジュール
module mod_monolis_converge_test
  use mod_monolis

  implicit none

contains

  subroutine monolis_converge_test
    implicit none

    !call monolis_set_converge_R(monoPRM, monoCOM, monoMAT, B, B2, is_converge, tdotp, tcomm)
    !call monolis_check_converge_R(monoPRM, monoCOM, monoMAT, R, B2, iter, is_converge, tdotp, tcomm)

  end subroutine monolis_converge_test

end module mod_monolis_converge_test
