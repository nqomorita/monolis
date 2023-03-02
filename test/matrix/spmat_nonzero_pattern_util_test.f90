!> 非零構造決定テストモジュール
module mod_monolis_spmat_nonzero_pattern_util_test
  use mod_monolis

  implicit none

contains

  subroutine monolis_spmat_nonzero_pattern_util_test
    implicit none

    call monolis_get_nonzero_pattern_by_nodal_graph_main_test()
    call monolis_get_nonzero_pattern_by_nodal_graph_with_arbit_main_test()
    call monolis_get_n_dof_index_test()
    call monolis_alloc_nonzero_pattern_mat_val_R_test()
    call monolis_alloc_nonzero_pattern_mat_val_C_test()
    call monolis_get_CSC_format_test()
  end subroutine monolis_spmat_nonzero_pattern_util_test

  subroutine monolis_get_nonzero_pattern_by_nodal_graph_main_test()
    implicit none

    !call monolis_get_nonzero_pattern_by_nodal_graph_main(monolis, n_node, ndof, index, item)
  end subroutine monolis_get_nonzero_pattern_by_nodal_graph_main_test

  subroutine monolis_get_nonzero_pattern_by_nodal_graph_with_arbit_main_test()
    implicit none

    !call monolis_get_nonzero_pattern_by_nodal_graph_with_arbit_main &
  end subroutine monolis_get_nonzero_pattern_by_nodal_graph_with_arbit_main_test

  subroutine monolis_get_n_dof_index_test()
    implicit none

    !call monolis_get_n_dof_index(n_node, n_dof_list, n_dof_index)
  end subroutine monolis_get_n_dof_index_test

  subroutine monolis_alloc_nonzero_pattern_mat_val_R_test()
    implicit none

    !call monolis_alloc_nonzero_pattern_mat_val_R(monolis)
  end subroutine monolis_alloc_nonzero_pattern_mat_val_R_test

  subroutine monolis_alloc_nonzero_pattern_mat_val_C_test()
    implicit none

    !call monolis_alloc_nonzero_pattern_mat_val_C(monolis)
  end subroutine monolis_alloc_nonzero_pattern_mat_val_C_test

  subroutine monolis_get_CSC_format_test()
    implicit none

    !call monolis_get_CSC_format(NC, NR, NZ, index, item, indexR, itemR, permR)
  end subroutine monolis_get_CSC_format_test
end module mod_monolis_spmat_nonzero_pattern_util_test
