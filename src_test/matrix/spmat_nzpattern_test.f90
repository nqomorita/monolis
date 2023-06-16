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

    call monolis_get_nonzero_pattern_by_simple_mesh_with_arbitrary_dof_R &
      & (monolis, n_node, n_base, n_dof_list, n_elem, elem)

    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 1", monolis%MAT%N, 6)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 2", monolis%MAT%NP, 6)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 3", monolis%MAT%NDOF, 1)

    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 4", monolis%MAT%CSR%index(1), 0)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 5", monolis%MAT%CSR%index(2), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 6", monolis%MAT%CSR%index(3), 7)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 7", monolis%MAT%CSR%index(4), 11)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 8", monolis%MAT%CSR%index(5), 16)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 9", monolis%MAT%CSR%index(6), 19)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 10", monolis%MAT%CSR%index(7), 22)

    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 11", monolis%MAT%CSR%item(1), 1)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 12", monolis%MAT%CSR%item(2), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 13", monolis%MAT%CSR%item(3), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 14", monolis%MAT%CSR%item(4), 1)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 15", monolis%MAT%CSR%item(5), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 16", monolis%MAT%CSR%item(6), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 17", monolis%MAT%CSR%item(7), 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 18", monolis%MAT%CSR%item(8), 1)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 19", monolis%MAT%CSR%item(9), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 20", monolis%MAT%CSR%item(10), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 21", monolis%MAT%CSR%item(11), 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 22", monolis%MAT%CSR%item(12), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 23", monolis%MAT%CSR%item(13), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 24", monolis%MAT%CSR%item(14), 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 25", monolis%MAT%CSR%item(15), 5)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 26", monolis%MAT%CSR%item(16), 6)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 27", monolis%MAT%CSR%item(17), 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 28", monolis%MAT%CSR%item(18), 5)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 29", monolis%MAT%CSR%item(19), 6)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 30", monolis%MAT%CSR%item(20), 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 31", monolis%MAT%CSR%item(21), 5)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 32", monolis%MAT%CSR%item(22), 6)

    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 33", monolis%MAT%CSC%index(1), 0)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 34", monolis%MAT%CSC%index(2), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 35", monolis%MAT%CSC%index(3), 7)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 36", monolis%MAT%CSC%index(4), 11)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 37", monolis%MAT%CSC%index(5), 16)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 38", monolis%MAT%CSC%index(6), 19)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 39", monolis%MAT%CSC%index(7), 22)

    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 40", monolis%MAT%CSC%item(1), 1)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 41", monolis%MAT%CSC%item(2), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 42", monolis%MAT%CSC%item(3), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 43", monolis%MAT%CSC%item(4), 1)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 44", monolis%MAT%CSC%item(5), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 45", monolis%MAT%CSC%item(6), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 46", monolis%MAT%CSC%item(7), 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 47", monolis%MAT%CSC%item(8), 1)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 48", monolis%MAT%CSC%item(9), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 49", monolis%MAT%CSC%item(10), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 50", monolis%MAT%CSC%item(11), 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 51", monolis%MAT%CSC%item(12), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 52", monolis%MAT%CSC%item(13), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 53", monolis%MAT%CSC%item(14), 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 54", monolis%MAT%CSC%item(15), 5)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 55", monolis%MAT%CSC%item(16), 6)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 56", monolis%MAT%CSC%item(17), 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 57", monolis%MAT%CSC%item(18), 5)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 58", monolis%MAT%CSC%item(19), 6)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 59", monolis%MAT%CSC%item(20), 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 60", monolis%MAT%CSC%item(21), 5)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 61", monolis%MAT%CSC%item(22), 6)

    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 62", monolis%MAT%CSC%perm(1), 1)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 63", monolis%MAT%CSC%perm(2), 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 64", monolis%MAT%CSC%perm(3), 8)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 65", monolis%MAT%CSC%perm(4), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 66", monolis%MAT%CSC%perm(5), 5)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 67", monolis%MAT%CSC%perm(6), 9)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 68", monolis%MAT%CSC%perm(7), 12)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 68", monolis%MAT%CSC%perm(8), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 70", monolis%MAT%CSC%perm(9), 6)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 71", monolis%MAT%CSC%perm(10), 10)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 72", monolis%MAT%CSC%perm(11), 13)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 73", monolis%MAT%CSC%perm(12), 7)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 74", monolis%MAT%CSC%perm(13), 11)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 75", monolis%MAT%CSC%perm(14), 14)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 76", monolis%MAT%CSC%perm(15), 17)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 77", monolis%MAT%CSC%perm(16), 20)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 78", monolis%MAT%CSC%perm(17), 15)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 79", monolis%MAT%CSC%perm(18), 18)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 80", monolis%MAT%CSC%perm(19), 21)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 81", monolis%MAT%CSC%perm(20), 16)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 82", monolis%MAT%CSC%perm(21), 19)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test 83", monolis%MAT%CSC%perm(22), 22)
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

    call monolis_get_nonzero_pattern_by_simple_mesh_with_arbitrary_dof_C &
      & (monolis, n_node, n_base, n_dof_list, n_elem, elem)

    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 1", monolis%MAT%N, 6)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 2", monolis%MAT%NP, 6)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 3", monolis%MAT%NDOF, 1)

    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 4", monolis%MAT%CSR%index(1), 0)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 5", monolis%MAT%CSR%index(2), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 6", monolis%MAT%CSR%index(3), 7)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 7", monolis%MAT%CSR%index(4), 11)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 8", monolis%MAT%CSR%index(5), 16)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 9", monolis%MAT%CSR%index(6), 19)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 10", monolis%MAT%CSR%index(7), 22)

    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 11", monolis%MAT%CSR%item(1), 1)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 12", monolis%MAT%CSR%item(2), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 13", monolis%MAT%CSR%item(3), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 14", monolis%MAT%CSR%item(4), 1)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 15", monolis%MAT%CSR%item(5), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 16", monolis%MAT%CSR%item(6), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 17", monolis%MAT%CSR%item(7), 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 18", monolis%MAT%CSR%item(8), 1)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 19", monolis%MAT%CSR%item(9), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 20", monolis%MAT%CSR%item(10), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 21", monolis%MAT%CSR%item(11), 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 22", monolis%MAT%CSR%item(12), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 23", monolis%MAT%CSR%item(13), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 24", monolis%MAT%CSR%item(14), 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 25", monolis%MAT%CSR%item(15), 5)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 26", monolis%MAT%CSR%item(16), 6)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 27", monolis%MAT%CSR%item(17), 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 28", monolis%MAT%CSR%item(18), 5)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 29", monolis%MAT%CSR%item(19), 6)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 30", monolis%MAT%CSR%item(20), 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 31", monolis%MAT%CSR%item(21), 5)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 32", monolis%MAT%CSR%item(22), 6)

    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 33", monolis%MAT%CSC%index(1), 0)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 34", monolis%MAT%CSC%index(2), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 35", monolis%MAT%CSC%index(3), 7)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 36", monolis%MAT%CSC%index(4), 11)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 37", monolis%MAT%CSC%index(5), 16)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 38", monolis%MAT%CSC%index(6), 19)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 39", monolis%MAT%CSC%index(7), 22)

    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 40", monolis%MAT%CSC%item(1), 1)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 41", monolis%MAT%CSC%item(2), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 42", monolis%MAT%CSC%item(3), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 43", monolis%MAT%CSC%item(4), 1)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 44", monolis%MAT%CSC%item(5), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 45", monolis%MAT%CSC%item(6), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 46", monolis%MAT%CSC%item(7), 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 47", monolis%MAT%CSC%item(8), 1)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 48", monolis%MAT%CSC%item(9), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 49", monolis%MAT%CSC%item(10), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 50", monolis%MAT%CSC%item(11), 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 51", monolis%MAT%CSC%item(12), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 52", monolis%MAT%CSC%item(13), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 53", monolis%MAT%CSC%item(14), 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 54", monolis%MAT%CSC%item(15), 5)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 55", monolis%MAT%CSC%item(16), 6)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 56", monolis%MAT%CSC%item(17), 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 57", monolis%MAT%CSC%item(18), 5)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 58", monolis%MAT%CSC%item(19), 6)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 59", monolis%MAT%CSC%item(20), 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 60", monolis%MAT%CSC%item(21), 5)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 61", monolis%MAT%CSC%item(22), 6)

    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 62", monolis%MAT%CSC%perm(1), 1)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 63", monolis%MAT%CSC%perm(2), 4)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 64", monolis%MAT%CSC%perm(3), 8)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 65", monolis%MAT%CSC%perm(4), 2)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 66", monolis%MAT%CSC%perm(5), 5)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 67", monolis%MAT%CSC%perm(6), 9)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 68", monolis%MAT%CSC%perm(7), 12)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 68", monolis%MAT%CSC%perm(8), 3)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 70", monolis%MAT%CSC%perm(9), 6)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 71", monolis%MAT%CSC%perm(10), 10)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 72", monolis%MAT%CSC%perm(11), 13)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 73", monolis%MAT%CSC%perm(12), 7)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 74", monolis%MAT%CSC%perm(13), 11)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 75", monolis%MAT%CSC%perm(14), 14)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 76", monolis%MAT%CSC%perm(15), 17)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 77", monolis%MAT%CSC%perm(16), 20)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 78", monolis%MAT%CSC%perm(17), 15)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 79", monolis%MAT%CSC%perm(18), 18)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 80", monolis%MAT%CSC%perm(19), 21)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 81", monolis%MAT%CSC%perm(20), 16)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 82", monolis%MAT%CSC%perm(21), 19)
    call monolis_test_check_eq_I1("monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test 83", monolis%MAT%CSC%perm(22), 22)
  end subroutine monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test
end module mod_monolis_spmat_nonzero_pattern_test
