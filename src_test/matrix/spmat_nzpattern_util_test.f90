!> 非零構造決定テストモジュール
module mod_monolis_spmat_nonzero_pattern_util_test
  use mod_monolis
  use mod_monolis_spmat_nonzero_pattern_util

  implicit none

contains

  subroutine monolis_spmat_nonzero_pattern_util_test()
    implicit none

    call monolis_get_nonzero_pattern_by_nodal_graph_main_test()
    call monolis_get_n_dof_index_test()
    call monolis_alloc_nonzero_pattern_mat_val_R_test()
    call monolis_alloc_nonzero_pattern_mat_val_C_test()
    call monolis_get_CSC_format_test()

    call monolis_std_global_log_string("monolis_stop_by_matrix_assemble")
    call monolis_std_global_log_string("monolis_stop_by_set_DBC")
    call monolis_std_global_log_string("monolis_stop_by_set_zero_diag_component")
    call monolis_std_global_log_string("monolis_stop_by_submatrix_access")
  end subroutine monolis_spmat_nonzero_pattern_util_test

  subroutine monolis_get_nonzero_pattern_by_nodal_graph_main_test()
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: n_node
    integer(kint) :: ndof
    integer(kint) :: index(5)
    integer(kint) :: item(6)

    call monolis_std_global_log_string("monolis_get_nonzero_pattern_by_nodal_graph_main")

    n_node = 4

    ndof = 2

    index(1) = 0
    index(2) = 1
    index(3) = 3
    index(4) = 5
    index(5) = 6

    item(1) = 2
    item(2) = 1
    item(3) = 3
    item(4) = 2
    item(5) = 4
    item(6) = 3

    call monolis_initialize(monolis)

    call monolis_get_nonzero_pattern_by_nodal_graph_main(monolis%MAT, n_node, ndof, index, item)

    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 1", monolis%MAT%N, 4)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 2", monolis%MAT%NP, 4)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 3", monolis%MAT%NDOF, 2)

    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 4", monolis%MAT%CSR%index(1), 0)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 5", monolis%MAT%CSR%index(2), 2)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 6", monolis%MAT%CSR%index(3), 5)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 7", monolis%MAT%CSR%index(4), 8)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 8", monolis%MAT%CSR%index(5), 10)

    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 9", monolis%MAT%CSR%item(1), 1)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 10", monolis%MAT%CSR%item(2), 2)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 11", monolis%MAT%CSR%item(3), 1)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 12", monolis%MAT%CSR%item(4), 2)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 13", monolis%MAT%CSR%item(5), 3)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 14", monolis%MAT%CSR%item(6), 2)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 15", monolis%MAT%CSR%item(7), 3)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 16", monolis%MAT%CSR%item(8), 4)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 17", monolis%MAT%CSR%item(9), 3)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 18", monolis%MAT%CSR%item(10), 4)

    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 19", monolis%MAT%CSC%index(1), 0)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 20", monolis%MAT%CSC%index(2), 2)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 21", monolis%MAT%CSC%index(3), 5)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 22", monolis%MAT%CSC%index(4), 8)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 23", monolis%MAT%CSC%index(5), 10)

    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 24", monolis%MAT%CSC%item(1), 1)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 25", monolis%MAT%CSC%item(2), 2)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 26", monolis%MAT%CSC%item(3), 1)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 27", monolis%MAT%CSC%item(4), 2)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 28", monolis%MAT%CSC%item(5), 3)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 29", monolis%MAT%CSC%item(6), 2)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 30", monolis%MAT%CSC%item(7), 3)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 31", monolis%MAT%CSC%item(8), 4)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 32", monolis%MAT%CSC%item(9), 3)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 33", monolis%MAT%CSC%item(10), 4)

    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 34", monolis%MAT%CSC%perm(1), 1)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 35", monolis%MAT%CSC%perm(2), 3)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 36", monolis%MAT%CSC%perm(3), 2)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 37", monolis%MAT%CSC%perm(4), 4)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 38", monolis%MAT%CSC%perm(5), 6)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 39", monolis%MAT%CSC%perm(6), 5)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 40", monolis%MAT%CSC%perm(7), 7)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 41", monolis%MAT%CSC%perm(8), 9)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 42", monolis%MAT%CSC%perm(9), 8)
    call monolis_test_check_eq_I1("monolis_get_nonzero_pattern_by_nodal_graph_main_test 43", monolis%MAT%CSC%perm(10), 10)
  end subroutine monolis_get_nonzero_pattern_by_nodal_graph_main_test

  subroutine monolis_get_n_dof_index_test()
    implicit none
    integer(kint) :: n_node
    integer(kint) :: n_dof_list(3)
    integer(kint) :: n_dof_index(3)

    call monolis_std_global_log_string("monolis_get_n_dof_index")

    n_node = 3

    n_dof_list(1) = 2
    n_dof_list(2) = 3
    n_dof_list(3) = 4

    n_dof_index(1) = 0
    n_dof_index(2) = 2
    n_dof_index(3) = 5

    !call monolis_get_n_dof_index(n_node, n_dof_list, n_dof_index)

    call monolis_test_check_eq_I1("monolis_get_n_dof_index_test 1", n_dof_index(1), 0)
    call monolis_test_check_eq_I1("monolis_get_n_dof_index_test 2", n_dof_index(2), 2)
    call monolis_test_check_eq_I1("monolis_get_n_dof_index_test 3", n_dof_index(3), 5)
  end subroutine monolis_get_n_dof_index_test

  subroutine monolis_alloc_nonzero_pattern_mat_val_R_test()
    implicit none
    type(monolis_structure) :: monolis

    call monolis_std_global_log_string("monolis_alloc_nonzero_pattern_mat_val_R")

    call monolis_initialize(monolis)

    monolis%MAT%NDOF = 2

    monolis%MAT%NP = 5

    call monolis_palloc_I_1d(monolis%MAT%CSR%index, 6)

    monolis%MAT%CSR%index(6) = 10

    call monolis_alloc_nonzero_pattern_mat_val_R(monolis%MAT)

    call monolis_test_check_eq_I1("monolis_alloc_nonzero_pattern_mat_val_R_test 1", size(monolis%MAT%R%A), 40)
    call monolis_test_check_eq_I1("monolis_alloc_nonzero_pattern_mat_val_R_test 2", size(monolis%MAT%R%B), 10)
    call monolis_test_check_eq_I1("monolis_alloc_nonzero_pattern_mat_val_R_test 3", size(monolis%MAT%R%X), 10)
  end subroutine monolis_alloc_nonzero_pattern_mat_val_R_test

  subroutine monolis_alloc_nonzero_pattern_mat_val_C_test()
    implicit none
    type(monolis_structure) :: monolis

    call monolis_std_global_log_string("monolis_alloc_nonzero_pattern_mat_val_C")

    call monolis_initialize(monolis)

    monolis%MAT%NDOF = 2

    monolis%MAT%NP = 5

    call monolis_palloc_I_1d(monolis%MAT%CSR%index, 6)

    monolis%MAT%CSR%index(6) = 10

    call monolis_alloc_nonzero_pattern_mat_val_C(monolis%MAT)

    call monolis_test_check_eq_I1("monolis_alloc_nonzero_pattern_mat_val_C_test 1", size(monolis%MAT%C%A), 40)
    call monolis_test_check_eq_I1("monolis_alloc_nonzero_pattern_mat_val_C_test 2", size(monolis%MAT%C%B), 10)
    call monolis_test_check_eq_I1("monolis_alloc_nonzero_pattern_mat_val_C_test 3", size(monolis%MAT%C%X), 10)
  end subroutine monolis_alloc_nonzero_pattern_mat_val_C_test

  subroutine monolis_get_CSC_format_test()
    implicit none
    integer(kint) :: NC
    integer(kint) :: NR
    integer(kint) :: NZ
    integer(kint) :: index(5)
    integer(kint) :: item(10)
    integer(kint) :: indexR(5)
    integer(kint) :: itemR(10)
    integer(kint) :: permR(10)

    call monolis_std_global_log_string("monolis_get_CSC_format")

    NC = 4

    NR = 4

    NZ = 10

    index(1) = 0
    index(2) = 2
    index(3) = 5
    index(4) = 8
    index(5) = 10

    item(1) = 1
    item(2) = 2
    item(3) = 1
    item(4) = 2
    item(5) = 3
    item(6) = 2
    item(7) = 3
    item(8) = 4
    item(9) = 3
    item(10) = 4

    call monolis_get_CSC_format(NC, NR, NZ, index, item, indexR, itemR, permR)

    call monolis_test_check_eq_I1("monolis_get_CSC_format_test 1", indexR(1), 0)
    call monolis_test_check_eq_I1("monolis_get_CSC_format_test 2", indexR(2), 2)
    call monolis_test_check_eq_I1("monolis_get_CSC_format_test 3", indexR(3), 5)
    call monolis_test_check_eq_I1("monolis_get_CSC_format_test 4", indexR(4), 8)
    call monolis_test_check_eq_I1("monolis_get_CSC_format_test 5", indexR(5), 10)

    call monolis_test_check_eq_I1("monolis_get_CSC_format_test 6", itemR(1), 1)
    call monolis_test_check_eq_I1("monolis_get_CSC_format_test 7", itemR(2), 2)
    call monolis_test_check_eq_I1("monolis_get_CSC_format_test 8", itemR(3), 1)
    call monolis_test_check_eq_I1("monolis_get_CSC_format_test 9", itemR(4), 2)
    call monolis_test_check_eq_I1("monolis_get_CSC_format_test 10", itemR(5), 3)
    call monolis_test_check_eq_I1("monolis_get_CSC_format_test 11", itemR(6), 2)
    call monolis_test_check_eq_I1("monolis_get_CSC_format_test 12", itemR(7), 3)
    call monolis_test_check_eq_I1("monolis_get_CSC_format_test 13", itemR(8), 4)
    call monolis_test_check_eq_I1("monolis_get_CSC_format_test 14", itemR(9), 3)
    call monolis_test_check_eq_I1("monolis_get_CSC_format_test 15", itemR(10), 4)

    call monolis_test_check_eq_I1("monolis_get_CSC_format_test 16", permR(1), 1)
    call monolis_test_check_eq_I1("monolis_get_CSC_format_test 17", permR(2), 3)
    call monolis_test_check_eq_I1("monolis_get_CSC_format_test 18", permR(3), 2)
    call monolis_test_check_eq_I1("monolis_get_CSC_format_test 19", permR(4), 4)
    call monolis_test_check_eq_I1("monolis_get_CSC_format_test 20", permR(5), 6)
    call monolis_test_check_eq_I1("monolis_get_CSC_format_test 21", permR(6), 5)
    call monolis_test_check_eq_I1("monolis_get_CSC_format_test 22", permR(7), 7)
    call monolis_test_check_eq_I1("monolis_get_CSC_format_test 23", permR(8), 9)
    call monolis_test_check_eq_I1("monolis_get_CSC_format_test 24", permR(9), 8)
    call monolis_test_check_eq_I1("monolis_get_CSC_format_test 25", permR(10), 10)
  end subroutine monolis_get_CSC_format_test
end module mod_monolis_spmat_nonzero_pattern_util_test
