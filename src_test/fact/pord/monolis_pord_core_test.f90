!> monolis_pord_core テストモジュール (スタブ)
module mod_monolis_pord_core_test
  use mod_monolis

  implicit none

contains

  subroutine monolis_pord_core_test
    implicit none

    call monolis_std_global_log_string("eval_sep")
    call monolis_std_global_log_string("pord_random")
    call monolis_std_global_log_string("minBucket")
    call monolis_std_global_log_string("pord_resettimer")
    call monolis_std_global_log_string("pord_starttimer")
    call monolis_std_global_log_string("pord_stoptimer")
    call monolis_std_global_log_string("distributionCounting")
    call monolis_std_global_log_string("insertUpInts")
    call monolis_std_global_log_string("insertUpIntsWithStaticIntKeys")
    call monolis_std_global_log_string("qsortUpInts")
    call monolis_std_global_log_string("setupBucket")
    call monolis_std_global_log_string("freeBucket")
    call monolis_std_global_log_string("insertBucket")
    call monolis_std_global_log_string("removeBucket")
  end subroutine monolis_pord_core_test

end module mod_monolis_pord_core_test
