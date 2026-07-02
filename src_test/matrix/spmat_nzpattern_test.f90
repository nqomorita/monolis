!> 非零構造決定テストモジュール
module mod_monolis_spmat_nonzero_pattern_test
  use mod_monolis
  use mod_monolis_spmat_nonzero_pattern
  use mod_monolis_spmat_handler

  implicit none

contains

  subroutine monolis_spmat_nonzero_pattern_test()
    implicit none

    call monolis_get_nonzero_pattern_by_simple_mesh_R_test()
    call monolis_get_nonzero_pattern_by_simple_mesh_C_test()

    call monolis_get_nzp_by_simple_mesh_with_arbit_dof_R_test()
    call monolis_get_nzp_by_simple_mesh_with_arbit_dof_C_test()

    call monolis_nonzero_pattern_with_margin_R_test()
    call monolis_nonzero_pattern_with_margin_C_test()
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

  subroutine monolis_nonzero_pattern_with_margin_R_test()
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: n_node, ndof, n_bandwidth, i
    integer(kint) :: index(5), item(6)
    integer(kint) :: index2(5), item2(6)
    real(kdouble) :: A(12)

    call monolis_std_global_log_string("monolis_alloc_nonzero_pattern_with_margin_R")
    call monolis_std_global_log_string("monolis_alloc_nonzero_pattern_with_margin_main")
    call monolis_std_global_log_string("monolis_update_nonzero_pattern_with_margin_by_nodal_graph_R")
    call monolis_std_global_log_string("monolis_update_nonzero_pattern_fixed_width_main")
    call monolis_std_global_log_string("monolis_set_array_to_sparse_matrix_R")

    n_node = 4
    ndof = 1
    n_bandwidth = 3

    call monolis_initialize(monolis)

    !> 固定バンド幅で確保
    call monolis_alloc_nonzero_pattern_with_margin_R(monolis, n_node, ndof, n_bandwidth)

    call monolis_test_check_eq_I1("margin_R 1", monolis%MAT%N, 4)
    call monolis_test_check_eq_I1("margin_R 2", monolis%MAT%NP, 4)
    call monolis_test_check_eq_I1("margin_R 3", monolis%MAT%NDOF, 1)

    !> 確保直後はサイズのみ（index/item の値は未設定）
    call monolis_test_check_eq_I1("margin_R 4", size(monolis%MAT%CSR%index), 5)
    call monolis_test_check_eq_I1("margin_R 5", size(monolis%MAT%CSR%item), 12)
    call monolis_test_check_eq_I1("margin_R 6", size(monolis%MAT%R%A), 12)
    call monolis_test_check_eq_I1("margin_R 7", size(monolis%MAT%R%B), 4)
    call monolis_test_check_eq_I1("margin_R 8", size(monolis%MAT%R%X), 4)

    !> 鎖状グラフ 1-2-3-4 で構造更新（バンド幅内なので再確保なし）
    index(1) = 0; index(2) = 1; index(3) = 3; index(4) = 5; index(5) = 6
    item(1) = 2
    item(2) = 1; item(3) = 3
    item(4) = 2; item(5) = 4
    item(6) = 3

    call monolis_update_nonzero_pattern_with_margin_by_nodal_graph_R(monolis, n_node, index, item)

    !> index（固定ストライド）
    call monolis_test_check_eq_I1("margin_R 9",  monolis%MAT%CSR%index(1), 0)
    call monolis_test_check_eq_I1("margin_R 10", monolis%MAT%CSR%index(2), 3)
    call monolis_test_check_eq_I1("margin_R 11", monolis%MAT%CSR%index(3), 6)
    call monolis_test_check_eq_I1("margin_R 12", monolis%MAT%CSR%index(4), 9)
    call monolis_test_check_eq_I1("margin_R 13", monolis%MAT%CSR%index(5), 12)

    !> item（対角＋近傍＋連続パディング、行内昇順ソート）
    call monolis_test_check_eq_I1("margin_R 14", monolis%MAT%CSR%item(1), 1)
    call monolis_test_check_eq_I1("margin_R 15", monolis%MAT%CSR%item(2), 2)
    call monolis_test_check_eq_I1("margin_R 16", monolis%MAT%CSR%item(3), 3)
    call monolis_test_check_eq_I1("margin_R 17", monolis%MAT%CSR%item(4), 1)
    call monolis_test_check_eq_I1("margin_R 18", monolis%MAT%CSR%item(5), 2)
    call monolis_test_check_eq_I1("margin_R 19", monolis%MAT%CSR%item(6), 3)
    call monolis_test_check_eq_I1("margin_R 20", monolis%MAT%CSR%item(7), 2)
    call monolis_test_check_eq_I1("margin_R 21", monolis%MAT%CSR%item(8), 3)
    call monolis_test_check_eq_I1("margin_R 22", monolis%MAT%CSR%item(9), 4)
    call monolis_test_check_eq_I1("margin_R 23", monolis%MAT%CSR%item(10), 1)
    call monolis_test_check_eq_I1("margin_R 24", monolis%MAT%CSR%item(11), 3)
    call monolis_test_check_eq_I1("margin_R 25", monolis%MAT%CSR%item(12), 4)

    !> CSC が再構築されていること
    call monolis_test_check_eq_I1("margin_R 26", monolis%MAT%CSC%index(5), 12)

    !> 値配列の一括設定
    do i = 1, 12
      A(i) = dble(i)
    enddo
    call monolis_set_array_to_sparse_matrix_R(monolis, A)
    call monolis_test_check_eq_R1("margin_R 27", monolis%MAT%R%A(1), 1.0d0)
    call monolis_test_check_eq_R1("margin_R 28", monolis%MAT%R%A(5), 5.0d0)
    call monolis_test_check_eq_R1("margin_R 29", monolis%MAT%R%A(12), 12.0d0)

    !> バンド幅を超える更新で再確保されること（星状グラフ、節点 2 の次数 3）
    index2(1) = 0; index2(2) = 1; index2(3) = 4; index2(4) = 5; index2(5) = 6
    item2(1) = 2
    item2(2) = 1; item2(3) = 3; item2(4) = 4
    item2(5) = 2
    item2(6) = 2

    call monolis_update_nonzero_pattern_with_margin_by_nodal_graph_R(monolis, n_node, index2, item2)

    call monolis_test_check_eq_I1("margin_R 30", size(monolis%MAT%CSR%item), 16)
    call monolis_test_check_eq_I1("margin_R 31", size(monolis%MAT%R%A), 16)
    call monolis_test_check_eq_I1("margin_R 32", monolis%MAT%CSR%index(2), 4)
    call monolis_test_check_eq_I1("margin_R 33", monolis%MAT%CSR%item(1), 1)
    call monolis_test_check_eq_I1("margin_R 34", monolis%MAT%CSR%item(4), 4)
    call monolis_test_check_eq_I1("margin_R 35", monolis%MAT%CSR%item(5), 1)
    call monolis_test_check_eq_I1("margin_R 36", monolis%MAT%CSR%item(8), 4)
    call monolis_test_check_eq_I1("margin_R 37", monolis%MAT%CSR%item(16), 4)

    call monolis_finalize(monolis)
  end subroutine monolis_nonzero_pattern_with_margin_R_test

  subroutine monolis_nonzero_pattern_with_margin_C_test()
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: n_node, ndof, n_bandwidth, i
    integer(kint) :: index(4), item(4)
    complex(kdouble) :: A(9)

    call monolis_std_global_log_string("monolis_alloc_nonzero_pattern_with_margin_C")
    call monolis_std_global_log_string("monolis_update_nonzero_pattern_with_margin_by_nodal_graph_C")
    call monolis_std_global_log_string("monolis_set_array_to_sparse_matrix_C")

    n_node = 3
    ndof = 1
    n_bandwidth = 3

    call monolis_initialize(monolis)

    call monolis_alloc_nonzero_pattern_with_margin_C(monolis, n_node, ndof, n_bandwidth)

    call monolis_test_check_eq_I1("margin_C 1", monolis%MAT%N, 3)
    call monolis_test_check_eq_I1("margin_C 2", monolis%MAT%NP, 3)
    call monolis_test_check_eq_I1("margin_C 3", size(monolis%MAT%CSR%item), 9)
    call monolis_test_check_eq_I1("margin_C 4", size(monolis%MAT%C%A), 9)

    !> 鎖状グラフ 1-2-3 で構造更新
    index(1) = 0; index(2) = 1; index(3) = 3; index(4) = 4
    item(1) = 2
    item(2) = 1; item(3) = 3
    item(4) = 2

    call monolis_update_nonzero_pattern_with_margin_by_nodal_graph_C(monolis, n_node, index, item)

    call monolis_test_check_eq_I1("margin_C 5", monolis%MAT%CSR%index(4), 9)
    call monolis_test_check_eq_I1("margin_C 6", monolis%MAT%CSR%item(1), 1)
    call monolis_test_check_eq_I1("margin_C 7", monolis%MAT%CSR%item(2), 2)
    call monolis_test_check_eq_I1("margin_C 8", monolis%MAT%CSR%item(3), 3)
    call monolis_test_check_eq_I1("margin_C 9", monolis%MAT%CSR%item(4), 1)
    call monolis_test_check_eq_I1("margin_C 10", monolis%MAT%CSR%item(5), 2)
    call monolis_test_check_eq_I1("margin_C 11", monolis%MAT%CSR%item(6), 3)
    call monolis_test_check_eq_I1("margin_C 12", monolis%MAT%CSR%item(7), 1)
    call monolis_test_check_eq_I1("margin_C 13", monolis%MAT%CSR%item(8), 2)
    call monolis_test_check_eq_I1("margin_C 14", monolis%MAT%CSR%item(9), 3)

    do i = 1, 9
      A(i) = cmplx(dble(i), -dble(i), kdouble)
    enddo
    call monolis_set_array_to_sparse_matrix_C(monolis, A)
    call monolis_test_check_eq_R1("margin_C 15", dble(monolis%MAT%C%A(1)),  1.0d0)
    call monolis_test_check_eq_R1("margin_C 16", aimag(monolis%MAT%C%A(1)), -1.0d0)
    call monolis_test_check_eq_R1("margin_C 17", dble(monolis%MAT%C%A(9)),  9.0d0)

    call monolis_finalize(monolis)
  end subroutine monolis_nonzero_pattern_with_margin_C_test
end module mod_monolis_spmat_nonzero_pattern_test
