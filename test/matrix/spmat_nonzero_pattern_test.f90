!> 非零構造決定テストモジュール
module mod_monolis_spmat_nonzero_pattern_test
  use mod_monolis

  implicit none

contains

  subroutine monolis_spmat_nonzero_pattern_test
    implicit none
    !call monolis_get_nonzero_pattern_by_simple_mesh_R(monolis, n_node, n_base, ndof, n_elem, elem)
    !call monolis_get_nonzero_pattern_by_simple_mesh_C(monolis, n_node, n_base, ndof, n_elem, elem)
    !call monolis_get_nonzero_pattern_by_simple_mesh_with_arbitrary_dof_R &
    !call monolis_get_nonzero_pattern_by_simple_mesh_with_arbitrary_dof_C &
    !call monolis_get_nonzero_pattern_by_connectivity_R &
    !call monolis_get_nonzero_pattern_by_connectivity_C &
    !call monolis_get_nonzero_pattern_by_nodal_graph_R(monolis, n_node, ndof, index, item)
    !call monolis_get_nonzero_pattern_by_nodal_graph_C(monolis, n_node, ndof, index, item)
    !call monolis_get_nonzero_pattern_by_nodal_graph_with_arbitrary_dof_R &
    !call monolis_get_nonzero_pattern_by_nodal_graph_with_arbitrary_dof_C &
  end subroutine monolis_spmat_nonzero_pattern_test

end module mod_monolis_spmat_nonzero_pattern_test
