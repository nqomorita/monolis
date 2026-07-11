!> monolis_pord_bisect テストモジュール (スタブ)
module mod_monolis_pord_bisect_test
  use mod_monolis

  implicit none

contains

  subroutine monolis_pord_bisect_test
    implicit none

    call monolis_std_global_log_string("freeBipartiteGraph")
    call monolis_std_global_log_string("setupBipartiteGraph")
    call monolis_std_global_log_string("maximumMatching")
    call monolis_std_global_log_string("maximumFlow")
    call monolis_std_global_log_string("DMviaMatching")
    call monolis_std_global_log_string("DMviaFlow")
    call monolis_std_global_log_string("findPseudoPeripheralDomain")
    call monolis_std_global_log_string("constructLevelSep")
    call monolis_std_global_log_string("initialDDSep")
    call monolis_std_global_log_string("updateB2W")
    call monolis_std_global_log_string("updateW2B")
    call monolis_std_global_log_string("improveDDSep")
    call monolis_std_global_log_string("newDomainDecomposition")
    call monolis_std_global_log_string("freeDomainDecomposition")
    call monolis_std_global_log_string("buildInitialDomains")
    call monolis_std_global_log_string("mergeMultisecs")
    call monolis_std_global_log_string("initialDomainDecomposition")
    call monolis_std_global_log_string("constructDomainDecomposition")
    call monolis_std_global_log_string("computePriorities")
    call monolis_std_global_log_string("random_number_int")
    call monolis_std_global_log_string("eliminateMultisecs")
    call monolis_std_global_log_string("findIndMultisecs")
    call monolis_std_global_log_string("coarserDomainDecomposition")
    call monolis_std_global_log_string("shrinkDomainDecomposition")
    call monolis_std_global_log_string("sort_by_key")
    call monolis_std_global_log_string("newGbisect")
    call monolis_std_global_log_string("freeGbisect")
    call monolis_std_global_log_string("constructSeparator")
    call monolis_std_global_log_string("smoothBy2Layers")
    call monolis_std_global_log_string("smoothSeparator")
  end subroutine monolis_pord_bisect_test

end module mod_monolis_pord_bisect_test
