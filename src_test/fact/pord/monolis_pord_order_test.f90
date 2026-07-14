!> monolis_pord_order テストモジュール (スタブ)
module mod_monolis_pord_order_test
  use mod_monolis

  implicit none

contains

  subroutine monolis_pord_order_test
    implicit none

    call monolis_std_global_log_string("freeElimGraph")
    call monolis_std_global_log_string("setupElimGraph")
    call monolis_std_global_log_string("crunchElimGraph")
    call monolis_std_global_log_string("buildElement")
    call monolis_std_global_log_string("grow_adjncy")
    call monolis_std_global_log_string("updateAdjncy")
    call monolis_std_global_log_string("findIndNodes")
    call monolis_std_global_log_string("updateDegree")
    call monolis_std_global_log_string("updateScore")
    call monolis_std_global_log_string("extractElimTree")
    call monolis_std_global_log_string("freeMinPriority")
    call monolis_std_global_log_string("setupMinPriority")
    call monolis_std_global_log_string("orderMinPriority")
    call monolis_std_global_log_string("eliminateStage")
    call monolis_std_global_log_string("eliminateStep")
    call monolis_std_global_log_string("stageinfo_sum")
  end subroutine monolis_pord_order_test

end module mod_monolis_pord_order_test
