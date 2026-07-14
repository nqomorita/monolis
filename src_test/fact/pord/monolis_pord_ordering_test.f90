!> monolis_pord_ordering テストモジュール (スタブ)
module mod_monolis_pord_ordering_test
  use mod_monolis

  implicit none

contains

  subroutine monolis_pord_ordering_test
    implicit none

    call monolis_std_global_log_string("monolis_pord_ordering")
    call monolis_std_global_log_string("monolis_pord_perm_from_elimtree")
  end subroutine monolis_pord_ordering_test

end module mod_monolis_pord_ordering_test
