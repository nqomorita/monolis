!> 非零構造決定テストモジュール
module mod_monolis_spmat_nonzero_pattern_test
  use mod_monolis
  use mod_monolis_spmat_nonzero_pattern

  implicit none

contains

  subroutine monolis_spmat_nonzero_pattern_test()
    implicit none

    call monolis_get_nonzero_pattern_by_simple_mesh_R_test()
    call monolis_get_nonzero_pattern_by_simple_mesh_C_test()

    call monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test()
    call monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test()
  end subroutine monolis_spmat_nonzero_pattern_test

  subroutine monolis_get_nonzero_pattern_by_simple_mesh_R_test()
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: n_node
    integer(kint) :: n_base
    integer(kint) :: ndof
    integer(kint) :: n_elem
    integer(kint) :: elem(2,3)

    call monolis_std_global_log_string("monolis_get_nonzero_pattern_by_simple_mesh_R")
    call monolis_std_global_log_string("monolis_get_nonzero_pattern_by_connectivity_R")
    call monolis_std_global_log_string("monolis_get_nonzero_pattern_by_nodal_graph_R")
    call monolis_std_global_log_string("monolis_get_nonzero_pattern_by_nodal_graph_with_arbitrary_dof_R")

    n_node = 4

    n_base = 2

    ndof = 2

    n_elem = 3

    elem(1,1) = 1; elem(2,1) = 2
    elem(1,2) = 2; elem(2,2) = 3
    elem(1,3) = 3; elem(2,3) = 4

    call monolis_initialize(monolis)

    call monolis_get_nonzero_pattern_by_simple_mesh_R(monolis, n_node, n_base, ndof, n_elem, elem)

    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 1", monolis%MAT%N, 4)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 2", monolis%MAT%NP, 4)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 3", monolis%MAT%NDOF, 2)

    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 4", monolis%MAT%CSR%index(1), 0)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 5", monolis%MAT%CSR%index(2), 2)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 6", monolis%MAT%CSR%index(3), 5)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 7", monolis%MAT%CSR%index(4), 8)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 8", monolis%MAT%CSR%index(5), 10)

    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 9", monolis%MAT%CSR%item(1), 1)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 10", monolis%MAT%CSR%item(2), 2)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 11", monolis%MAT%CSR%item(3), 1)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 12", monolis%MAT%CSR%item(4), 2)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 13", monolis%MAT%CSR%item(5), 3)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 14", monolis%MAT%CSR%item(6), 2)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 15", monolis%MAT%CSR%item(7), 3)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 16", monolis%MAT%CSR%item(8), 4)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 17", monolis%MAT%CSR%item(9), 3)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 18", monolis%MAT%CSR%item(10), 4)

    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 19", monolis%MAT%CSC%index(1), 0)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 20", monolis%MAT%CSC%index(2), 2)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 21", monolis%MAT%CSC%index(3), 5)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 22", monolis%MAT%CSC%index(4), 8)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 23", monolis%MAT%CSC%index(5), 10)

    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 24", monolis%MAT%CSC%item(1), 1)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 25", monolis%MAT%CSC%item(2), 2)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 26", monolis%MAT%CSC%item(3), 1)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 27", monolis%MAT%CSC%item(4), 2)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 28", monolis%MAT%CSC%item(5), 3)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 29", monolis%MAT%CSC%item(6), 2)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 30", monolis%MAT%CSC%item(7), 3)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 31", monolis%MAT%CSC%item(8), 4)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 32", monolis%MAT%CSC%item(9), 3)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 33", monolis%MAT%CSC%item(10), 4)

    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 34", monolis%MAT%CSC%perm(1), 1)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 35", monolis%MAT%CSC%perm(2), 3)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 36", monolis%MAT%CSC%perm(3), 2)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 37", monolis%MAT%CSC%perm(4), 4)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 38", monolis%MAT%CSC%perm(5), 6)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 39", monolis%MAT%CSC%perm(6), 5)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 40", monolis%MAT%CSC%perm(7), 7)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 41", monolis%MAT%CSC%perm(8), 9)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 42", monolis%MAT%CSC%perm(9), 8)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 43", monolis%MAT%CSC%perm(10), 10)

    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 44", size(monolis%MAT%R%A), 40)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 45", size(monolis%MAT%R%B), 8)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 46", size(monolis%MAT%R%X), 8)
  end subroutine monolis_get_nonzero_pattern_by_simple_mesh_R_test

  subroutine monolis_get_nonzero_pattern_by_simple_mesh_C_test()
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: n_node
    integer(kint) :: n_base
    integer(kint) :: ndof
    integer(kint) :: n_elem
    integer(kint) :: elem(2,3)

    call monolis_std_global_log_string("monolis_get_nonzero_pattern_by_simple_mesh_C")
    call monolis_std_global_log_string("monolis_get_nonzero_pattern_by_connectivity_C")
    call monolis_std_global_log_string("monolis_get_nonzero_pattern_by_nodal_graph_C")
    call monolis_std_global_log_string("monolis_get_nonzero_pattern_by_nodal_graph_with_arbitrary_dof_C")

    n_node = 4

    n_base = 2

    ndof = 2

    n_elem = 3

    elem(1,1) = 1; elem(2,1) = 2
    elem(1,2) = 2; elem(2,2) = 3
    elem(1,3) = 3; elem(2,3) = 4

    call monolis_initialize(monolis)

    call monolis_get_nonzero_pattern_by_simple_mesh_C(monolis, n_node, n_base, ndof, n_elem, elem)

    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 1", monolis%MAT%N, 4)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 2", monolis%MAT%NP, 4)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 3", monolis%MAT%NDOF, 2)

    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 4", monolis%MAT%CSR%index(1), 0)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 5", monolis%MAT%CSR%index(2), 2)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 6", monolis%MAT%CSR%index(3), 5)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 7", monolis%MAT%CSR%index(4), 8)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 8", monolis%MAT%CSR%index(5), 10)

    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 9", monolis%MAT%CSR%item(1), 1)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 10", monolis%MAT%CSR%item(2), 2)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 11", monolis%MAT%CSR%item(3), 1)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 12", monolis%MAT%CSR%item(4), 2)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 13", monolis%MAT%CSR%item(5), 3)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 14", monolis%MAT%CSR%item(6), 2)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 15", monolis%MAT%CSR%item(7), 3)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 16", monolis%MAT%CSR%item(8), 4)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 17", monolis%MAT%CSR%item(9), 3)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 18", monolis%MAT%CSR%item(10), 4)

    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 19", monolis%MAT%CSC%index(1), 0)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 20", monolis%MAT%CSC%index(2), 2)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 21", monolis%MAT%CSC%index(3), 5)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 22", monolis%MAT%CSC%index(4), 8)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 23", monolis%MAT%CSC%index(5), 10)

    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 24", monolis%MAT%CSC%item(1), 1)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 25", monolis%MAT%CSC%item(2), 2)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 26", monolis%MAT%CSC%item(3), 1)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 27", monolis%MAT%CSC%item(4), 2)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 28", monolis%MAT%CSC%item(5), 3)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 29", monolis%MAT%CSC%item(6), 2)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 30", monolis%MAT%CSC%item(7), 3)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 31", monolis%MAT%CSC%item(8), 4)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 32", monolis%MAT%CSC%item(9), 3)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 33", monolis%MAT%CSC%item(10), 4)

    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 34", monolis%MAT%CSC%perm(1), 1)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 35", monolis%MAT%CSC%perm(2), 3)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 36", monolis%MAT%CSC%perm(3), 2)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 37", monolis%MAT%CSC%perm(4), 4)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 38", monolis%MAT%CSC%perm(5), 6)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 39", monolis%MAT%CSC%perm(6), 5)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 40", monolis%MAT%CSC%perm(7), 7)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 41", monolis%MAT%CSC%perm(8), 9)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 42", monolis%MAT%CSC%perm(9), 8)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 43", monolis%MAT%CSC%perm(10), 10)

    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 44", size(monolis%MAT%C%A), 40)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_C_test 45", size(monolis%MAT%C%B), 8)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_simple_mesh_R_test 46", size(monolis%MAT%C%X), 8)
  end subroutine monolis_get_nonzero_pattern_by_simple_mesh_C_test

  subroutine monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test()
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: n_node
    integer(kint) :: n_base
    integer(kint) :: n_dof_list(4)
    integer(kint) :: n_elem
    integer(kint) :: elem(2,3)

    call monolis_std_global_log_string("monolis_get_nonzero_pattern_by_simple_mesh_with_arbitrary_dof_R")

    n_node = 4

    n_base = 2

    n_elem = 3

    elem(1,1) = 1; elem(2,1) = 2
    elem(1,2) = 2; elem(2,2) = 3
    elem(1,3) = 3; elem(2,3) = 4

    n_dof_list(1) = 1
    n_dof_list(2) = 2
    n_dof_list(3) = 1
    n_dof_list(4) = 2

    call monolis_initialize(monolis)

    call monolis_get_nonzero_pattern_by_simple_mesh_V_R &
      & (monolis, n_node, n_base, n_dof_list, n_elem, elem)

    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 1", monolis%MAT%N, 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 2", monolis%MAT%NP, 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 3", monolis%MAT%NDOF, -1)

    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 4", monolis%MAT%CSR%index(1), 0)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 5", monolis%MAT%CSR%index(2), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 6", monolis%MAT%CSR%index(3), 5)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 7", monolis%MAT%CSR%index(4), 8)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 8", monolis%MAT%CSR%index(5), 10)

    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 9", monolis%MAT%CSR%item(1), 1)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 10", monolis%MAT%CSR%item(2), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 11", monolis%MAT%CSR%item(3), 1)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 12", monolis%MAT%CSR%item(4), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 13", monolis%MAT%CSR%item(5), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 14", monolis%MAT%CSR%item(6), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 15", monolis%MAT%CSR%item(7), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 16", monolis%MAT%CSR%item(8), 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 17", monolis%MAT%CSR%item(9), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 18", monolis%MAT%CSR%item(10), 4)

    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 19", monolis%MAT%CSC%index(1), 0)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 20", monolis%MAT%CSC%index(2), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 21", monolis%MAT%CSC%index(3), 5)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 22", monolis%MAT%CSC%index(4), 8)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 23", monolis%MAT%CSC%index(5), 10)

    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 24", monolis%MAT%CSC%item(1), 1)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 25", monolis%MAT%CSC%item(2), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 26", monolis%MAT%CSC%item(3), 1)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 27", monolis%MAT%CSC%item(4), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 28", monolis%MAT%CSC%item(5), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 29", monolis%MAT%CSC%item(6), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 30", monolis%MAT%CSC%item(7), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 31", monolis%MAT%CSC%item(8), 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 32", monolis%MAT%CSC%item(9), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 33", monolis%MAT%CSC%item(10), 4)

    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 34", monolis%MAT%CSC%perm(1), 1)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 35", monolis%MAT%CSC%perm(2), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 36", monolis%MAT%CSC%perm(3), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 37", monolis%MAT%CSC%perm(4), 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 38", monolis%MAT%CSC%perm(5), 6)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 39", monolis%MAT%CSC%perm(6), 5)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 40", monolis%MAT%CSC%perm(7), 7)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 41", monolis%MAT%CSC%perm(8), 9)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 42", monolis%MAT%CSC%perm(9), 8)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 43", monolis%MAT%CSC%perm(10), 10)

    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test i1", monolis%MAT%n_dof_list(1), 1)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test i2", monolis%MAT%n_dof_list(2), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test i3", monolis%MAT%n_dof_list(3), 1)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test i4", monolis%MAT%n_dof_list(4), 2)

    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test j1", monolis%MAT%n_dof_index(1), 0)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test j2", monolis%MAT%n_dof_index(2), 1)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test j3", monolis%MAT%n_dof_index(3), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test j4", monolis%MAT%n_dof_index(4), 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test j5", monolis%MAT%n_dof_index(5), 6)

    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test k1",  monolis%MAT%n_dof_index2(1), 0)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test k2",  monolis%MAT%n_dof_index2(2), 1)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test k3",  monolis%MAT%n_dof_index2(3), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test k4",  monolis%MAT%n_dof_index2(4), 5)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test k5",  monolis%MAT%n_dof_index2(5), 9)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test k6",  monolis%MAT%n_dof_index2(6), 11)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test k7",  monolis%MAT%n_dof_index2(7), 13)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test k8",  monolis%MAT%n_dof_index2(8), 14)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test k9",  monolis%MAT%n_dof_index2(9), 16)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test k10", monolis%MAT%n_dof_index2(10),18)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test k11", monolis%MAT%n_dof_index2(11),22)

    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 44", size(monolis%MAT%R%A), 22)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 45", size(monolis%MAT%R%B), 6)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 46", size(monolis%MAT%R%X), 6)
  end subroutine monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test

  subroutine monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test()
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: n_node
    integer(kint) :: n_base
    integer(kint) :: n_dof_list(4)
    integer(kint) :: n_elem
    integer(kint) :: elem(2,3)

    call monolis_std_global_log_string("monolis_get_nonzero_pattern_by_simple_mesh_with_arbitrary_dof_C")

    n_node = 4

    n_base = 2

    n_elem = 3

    elem(1,1) = 1; elem(2,1) = 2
    elem(1,2) = 2; elem(2,2) = 3
    elem(1,3) = 3; elem(2,3) = 4

    n_dof_list(1) = 1
    n_dof_list(2) = 2
    n_dof_list(3) = 1
    n_dof_list(4) = 2

    call monolis_initialize(monolis)

    call monolis_get_nonzero_pattern_by_simple_mesh_V_C &
      & (monolis, n_node, n_base, n_dof_list, n_elem, elem)

    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 1", monolis%MAT%N, 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 2", monolis%MAT%NP, 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 3", monolis%MAT%NDOF, -1)

    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 4", monolis%MAT%CSR%index(1), 0)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 5", monolis%MAT%CSR%index(2), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 6", monolis%MAT%CSR%index(3), 5)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 7", monolis%MAT%CSR%index(4), 8)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 8", monolis%MAT%CSR%index(5), 10)

    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 9", monolis%MAT%CSR%item(1), 1)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 10", monolis%MAT%CSR%item(2), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 11", monolis%MAT%CSR%item(3), 1)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 12", monolis%MAT%CSR%item(4), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 13", monolis%MAT%CSR%item(5), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 14", monolis%MAT%CSR%item(6), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 15", monolis%MAT%CSR%item(7), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 16", monolis%MAT%CSR%item(8), 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 17", monolis%MAT%CSR%item(9), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 18", monolis%MAT%CSR%item(10), 4)

    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 19", monolis%MAT%CSC%index(1), 0)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 20", monolis%MAT%CSC%index(2), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 21", monolis%MAT%CSC%index(3), 5)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 22", monolis%MAT%CSC%index(4), 8)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 23", monolis%MAT%CSC%index(5), 10)

    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 24", monolis%MAT%CSC%item(1), 1)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 25", monolis%MAT%CSC%item(2), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 26", monolis%MAT%CSC%item(3), 1)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 27", monolis%MAT%CSC%item(4), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 28", monolis%MAT%CSC%item(5), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 29", monolis%MAT%CSC%item(6), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 30", monolis%MAT%CSC%item(7), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 31", monolis%MAT%CSC%item(8), 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 32", monolis%MAT%CSC%item(9), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 33", monolis%MAT%CSC%item(10), 4)

    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 34", monolis%MAT%CSC%perm(1), 1)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 35", monolis%MAT%CSC%perm(2), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 36", monolis%MAT%CSC%perm(3), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 37", monolis%MAT%CSC%perm(4), 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 38", monolis%MAT%CSC%perm(5), 6)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 39", monolis%MAT%CSC%perm(6), 5)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 40", monolis%MAT%CSC%perm(7), 7)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 41", monolis%MAT%CSC%perm(8), 9)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 42", monolis%MAT%CSC%perm(9), 8)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 43", monolis%MAT%CSC%perm(10), 10)

    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test i1", monolis%MAT%n_dof_list(1), 1)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test i2", monolis%MAT%n_dof_list(2), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test i3", monolis%MAT%n_dof_list(3), 1)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test i4", monolis%MAT%n_dof_list(4), 2)

    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test j1", monolis%MAT%n_dof_index(1), 0)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test j2", monolis%MAT%n_dof_index(2), 1)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test j3", monolis%MAT%n_dof_index(3), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test j4", monolis%MAT%n_dof_index(4), 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test j5", monolis%MAT%n_dof_index(5), 6)

    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test k1",  monolis%MAT%n_dof_index2(1), 0)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test k2",  monolis%MAT%n_dof_index2(2), 1)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test k3",  monolis%MAT%n_dof_index2(3), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test k4",  monolis%MAT%n_dof_index2(4), 5)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test k5",  monolis%MAT%n_dof_index2(5), 9)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test k6",  monolis%MAT%n_dof_index2(6), 11)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test k7",  monolis%MAT%n_dof_index2(7), 13)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test k8",  monolis%MAT%n_dof_index2(8), 14)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test k9",  monolis%MAT%n_dof_index2(9), 16)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test k10", monolis%MAT%n_dof_index2(10),18)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test k11", monolis%MAT%n_dof_index2(11),22)

    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 44", size(monolis%MAT%C%A), 22)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 45", size(monolis%MAT%C%B), 6)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 46", size(monolis%MAT%C%X), 6)
  end subroutine monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test
end module mod_monolis_spmat_nonzero_pattern_test
