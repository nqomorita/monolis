!> monolis_pord_decomp テストモジュール (スタブ)
module mod_monolis_pord_decomp_test
  use mod_monolis

  implicit none

contains

  subroutine monolis_pord_decomp_test
    implicit none

    call monolis_std_global_log_string("newNDnode")
    call monolis_std_global_log_string("freeNDnode")
    call monolis_std_global_log_string("freeNDtree")
    call monolis_std_global_log_string("setupNDroot")
    call monolis_std_global_log_string("splitNDnode")
    call monolis_std_global_log_string("buildNDtree")
    call monolis_std_global_log_string("newMultisector")
    call monolis_std_global_log_string("freeMultisector")
    call monolis_std_global_log_string("trivialMultisector")
    call monolis_std_global_log_string("constructMultisector")
    call monolis_std_global_log_string("extractMS2stage")
    call monolis_std_global_log_string("extractMSmultistage")
  end subroutine monolis_pord_decomp_test

end module mod_monolis_pord_decomp_test
