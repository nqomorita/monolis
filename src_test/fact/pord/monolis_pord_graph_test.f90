!> monolis_pord_graph テストモジュール (スタブ)
module mod_monolis_pord_graph_test
  use mod_monolis

  implicit none

contains

  subroutine monolis_pord_graph_test
    implicit none

    call monolis_std_global_log_string("newGraph")
    call monolis_std_global_log_string("freeGraph")
    call monolis_std_global_log_string("setupSubgraph")
    call monolis_std_global_log_string("compressGraph")
    call monolis_std_global_log_string("indNodes")
    call monolis_std_global_log_string("setupGridGraph")
    call monolis_std_global_log_string("newElimTree")
    call monolis_std_global_log_string("freeElimTree")
    call monolis_std_global_log_string("initFchSilbRoot")
    call monolis_std_global_log_string("firstPostorder")
    call monolis_std_global_log_string("nextPostorder")
    call monolis_std_global_log_string("permFromElimTree")
    call monolis_std_global_log_string("expandElimTree")
  end subroutine monolis_pord_graph_test

end module mod_monolis_pord_graph_test
