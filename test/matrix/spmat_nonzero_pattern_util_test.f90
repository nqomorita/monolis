!> 非零構造決定テストモジュール
module mod_monolis_spmat_nonzero_pattern_util_test
  use mod_monolis

  implicit none

contains

  subroutine monolis_spmat_nonzero_pattern_util_test
    implicit none

    !call monolis_get_nonzero_pattern_by_nodal_graph_main(monolis, n_node, ndof, index, item)
    !call monolis_get_nonzero_pattern_by_nodal_graph_with_arbit_main &
    !call monolis_get_n_dof_index(n_node, n_dof_list, n_dof_index)
    !call monolis_alloc_nonzero_pattern_mat_val_R(monolis)
    !call monolis_alloc_nonzero_pattern_mat_val_C(monolis)
    !call monolis_get_CSC_format(NC, NR, NZ, index, item, indexR, itemR, permR)
  end subroutine monolis_spmat_nonzero_pattern_util_test

end module mod_monolis_spmat_nonzero_pattern_util_test
