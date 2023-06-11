!> 疎行列操作関数群
module mod_monolis_spmat_handler_test
  use mod_monolis
  use mod_monolis_spmat_handler
  use mod_monolis_spmat_nonzero_pattern

  implicit none

contains

  subroutine monolis_spmat_handler_test()
    implicit none

    call monolis_set_scalar_to_sparse_matrix_R_test()
    call monolis_set_scalar_to_sparse_matrix_C_test()
  end subroutine monolis_spmat_handler_test

  subroutine monolis_set_scalar_to_sparse_matrix_R_test()
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: i
    integer(kint) :: j
    integer(kint) :: sub_i
    integer(kint) :: sub_j
    real(kdouble) :: val
    real(kdouble) :: bval(2,2)
    real(kdouble) :: B(8)
    integer(kint) :: n_node
    integer(kint) :: n_base
    integer(kint) :: ndof
    integer(kint) :: n_elem
    integer(kint) :: e(1), f(1)
    integer(kint) :: elem(2,3)
    logical :: is_find

    call monolis_std_global_log_string("monolis_set_scalar_to_sparse_matrix_R")
    call monolis_std_global_log_string("monolis_add_scalar_to_sparse_matrix_R")
    call monolis_std_global_log_string("monolis_set_block_to_sparse_matrix_R")
    call monolis_std_global_log_string("monolis_add_matrix_to_sparse_matrix_R")
    call monolis_std_global_log_string("monolis_add_matrix_to_sparse_matrix_main_R")
    call monolis_std_global_log_string("monolis_add_matrix_to_sparse_matrix_offdiag_R")
    call monolis_std_global_log_string("monolis_set_Dirichlet_bc_R")
    call monolis_std_global_log_string("monolis_set_Dirichlet_bc_main_R")

    n_node = 4

    n_base = 2

    ndof = 2

    n_elem = 3

    elem(1,1) = 1; elem(2,1) = 2
    elem(1,2) = 2; elem(2,2) = 3
    elem(1,3) = 3; elem(2,3) = 4

    call monolis_initialize(monolis)

    call monolis_get_nonzero_pattern_by_simple_mesh_R(monolis, n_node, n_base, ndof, n_elem, elem)

    val = 1.0d0

    call monolis_set_scalar_to_sparse_matrix_R(monolis, 2, 2, 1, 1, val)

    call monolis_test_check_eq_R1("monolis_set_scalar_to_sparse_matrix_R_test 1", monolis%MAT%R%A(13), 1.0d0)

    val = 0.0d0

    bval(1,1) = 1.0d0
    bval(1,2) = 2.0d0
    bval(2,1) = 3.0d0
    bval(2,2) = 4.0d0

    call monolis_set_block_to_sparse_matrix_R(monolis, 2, 2, bval)

    call monolis_test_check_eq_R1("monolis_set_scalar_to_sparse_matrix_R_test 2a", monolis%MAT%R%A(13), 1.0d0)
    call monolis_test_check_eq_R1("monolis_set_scalar_to_sparse_matrix_R_test 2b", monolis%MAT%R%A(14), 2.0d0)
    call monolis_test_check_eq_R1("monolis_set_scalar_to_sparse_matrix_R_test 2c", monolis%MAT%R%A(15), 3.0d0)
    call monolis_test_check_eq_R1("monolis_set_scalar_to_sparse_matrix_R_test 2d", monolis%MAT%R%A(16), 4.0d0)

    call monolis_get_scalar_from_sparse_matrix_R(monolis, 2, 2, 1, 1, val, is_find)

    call monolis_test_check_eq_R1("monolis_set_scalar_to_sparse_matrix_R_test 3", val, 1.0d0)

    val = 1.0d0

    call monolis_add_scalar_to_sparse_matrix_R(monolis, 2, 2, 1, 1, val)

    call monolis_test_check_eq_R1("monolis_set_scalar_to_sparse_matrix_R_test 3", monolis%MAT%R%A(13), 2.0d0)

    bval(1,1) = 1.0d0
    bval(1,2) = 2.0d0
    bval(2,1) = 3.0d0
    bval(2,2) = 4.0d0

    e(1) = 2

    call monolis_add_matrix_to_sparse_matrix_R(monolis, 1, e, bval)

    call monolis_test_check_eq_R1("monolis_set_scalar_to_sparse_matrix_R_test 4a", monolis%MAT%R%A(13), 3.0d0)
    call monolis_test_check_eq_R1("monolis_set_scalar_to_sparse_matrix_R_test 4b", monolis%MAT%R%A(14), 4.0d0)
    call monolis_test_check_eq_R1("monolis_set_scalar_to_sparse_matrix_R_test 4c", monolis%MAT%R%A(15), 6.0d0)
    call monolis_test_check_eq_R1("monolis_set_scalar_to_sparse_matrix_R_test 4d", monolis%MAT%R%A(16), 8.0d0)

    bval(1,1) = 1.0d0
    bval(1,2) = 2.0d0
    bval(2,1) = 3.0d0
    bval(2,2) = 4.0d0

    e(1) = 2
    f(1) = 3

    call monolis_add_matrix_to_sparse_matrix_offdiag_R(monolis, 1, 1, e, f, bval)

    call monolis_test_check_eq_R1("monolis_set_scalar_to_sparse_matrix_R_test 5a", monolis%MAT%R%A(17), 1.0d0)
    call monolis_test_check_eq_R1("monolis_set_scalar_to_sparse_matrix_R_test 5b", monolis%MAT%R%A(18), 2.0d0)
    call monolis_test_check_eq_R1("monolis_set_scalar_to_sparse_matrix_R_test 5c", monolis%MAT%R%A(19), 3.0d0)
    call monolis_test_check_eq_R1("monolis_set_scalar_to_sparse_matrix_R_test 5d", monolis%MAT%R%A(20), 4.0d0)

    monolis%MAT%R%A = 2.0d0

    B = 0.0d0

    call monolis_set_Dirichlet_bc_R(monolis, B, 2, 1, 1.0d0)

    call monolis_test_check_eq_R1("monolis_set_scalar_to_sparse_matrix_R_test 6a", B(1), -2.0d0)
    call monolis_test_check_eq_R1("monolis_set_scalar_to_sparse_matrix_R_test 6b", B(2), -2.0d0)
    call monolis_test_check_eq_R1("monolis_set_scalar_to_sparse_matrix_R_test 6c", B(3),  1.0d0)
    call monolis_test_check_eq_R1("monolis_set_scalar_to_sparse_matrix_R_test 6d", B(4), -2.0d0)
    call monolis_test_check_eq_R1("monolis_set_scalar_to_sparse_matrix_R_test 6e", B(5), -2.0d0)
    call monolis_test_check_eq_R1("monolis_set_scalar_to_sparse_matrix_R_test 6f", B(6), -2.0d0)
    call monolis_test_check_eq_R1("monolis_set_scalar_to_sparse_matrix_R_test 6g", B(7),  0.0d0)
    call monolis_test_check_eq_R1("monolis_set_scalar_to_sparse_matrix_R_test 6h", B(8),  0.0d0)
  end subroutine monolis_set_scalar_to_sparse_matrix_R_test

  subroutine monolis_set_scalar_to_sparse_matrix_C_test()
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: i
    integer(kint) :: j
    integer(kint) :: sub_i
    integer(kint) :: sub_j
    complex(kdouble) :: val
    complex(kdouble) :: bval(2,2)
    complex(kdouble) :: B(8)
    integer(kint) :: n_node
    integer(kint) :: n_base
    integer(kint) :: ndof
    integer(kint) :: n_elem
    integer(kint) :: e(1), f(1)
    integer(kint) :: elem(2,3)
    logical :: is_find

    call monolis_std_global_log_string("monolis_set_scalar_to_sparse_matrix_C")
    call monolis_std_global_log_string("monolis_add_scalar_to_sparse_matrix_C")
    call monolis_std_global_log_string("monolis_set_block_to_sparse_matrix_C")
    call monolis_std_global_log_string("monolis_add_matrix_to_sparse_matrix_C")
    call monolis_std_global_log_string("monolis_add_matrix_to_sparse_matrix_main_C")
    call monolis_std_global_log_string("monolis_add_matrix_to_sparse_matrix_offdiag_C")
    call monolis_std_global_log_string("monolis_set_Dirichlet_bc_C")
    call monolis_std_global_log_string("monolis_set_Dirichlet_bc_main_C")

    n_node = 4

    n_base = 2

    ndof = 2

    n_elem = 3

    elem(1,1) = 1; elem(2,1) = 2
    elem(1,2) = 2; elem(2,2) = 3
    elem(1,3) = 3; elem(2,3) = 4

    call monolis_initialize(monolis)

    call monolis_get_nonzero_pattern_by_simple_mesh_C(monolis, n_node, n_base, ndof, n_elem, elem)

    val = (1.0d0, 1.0d0)

    call monolis_set_scalar_to_sparse_matrix_C(monolis, 2, 2, 1, 1, val)

    call monolis_test_check_eq_C1("monolis_set_scalar_to_sparse_matrix_C_test 1", monolis%MAT%C%A(13), (1.0d0, 1.0d0))

    val = 0.0d0

    bval(1,1) = (1.0d0, 1.0d0)
    bval(1,2) = (2.0d0, 2.0d0)
    bval(2,1) = (3.0d0, 3.0d0)
    bval(2,2) = (4.0d0, 4.0d0)

    call monolis_set_block_to_sparse_matrix_C(monolis, 2, 2, bval)

    call monolis_test_check_eq_C1("monolis_set_scalar_to_sparse_matrix_C_test 2a", monolis%MAT%C%A(13), (1.0d0, 1.0d0))
    call monolis_test_check_eq_C1("monolis_set_scalar_to_sparse_matrix_C_test 2b", monolis%MAT%C%A(14), (2.0d0, 2.0d0))
    call monolis_test_check_eq_C1("monolis_set_scalar_to_sparse_matrix_C_test 2c", monolis%MAT%C%A(15), (3.0d0, 3.0d0))
    call monolis_test_check_eq_C1("monolis_set_scalar_to_sparse_matrix_C_test 2d", monolis%MAT%C%A(16), (4.0d0, 4.0d0))

    call monolis_get_scalar_from_sparse_matrix_C(monolis, 2, 2, 1, 1, val, is_find)

    call monolis_test_check_eq_C1("monolis_set_scalar_to_sparse_matrix_C_test 3", val, (1.0d0, 1.0d0))

    val = (1.0d0, 1.0d0)

    call monolis_add_scalar_to_sparse_matrix_C(monolis, 2, 2, 1, 1, val)

    call monolis_test_check_eq_C1("monolis_set_scalar_to_sparse_matrix_C_test 3", monolis%MAT%C%A(13), (2.0d0, 2.0d0))

    bval(1,1) = (1.0d0, 1.0d0)
    bval(1,2) = (2.0d0, 2.0d0)
    bval(2,1) = (3.0d0, 3.0d0)
    bval(2,2) = (4.0d0, 4.0d0)

    e(1) = 2

    call monolis_add_matrix_to_sparse_matrix_C(monolis, 1, e, bval)

    call monolis_test_check_eq_C1("monolis_set_scalar_to_sparse_matrix_C_test 4a", monolis%MAT%C%A(13), (3.0d0, 3.0d0))
    call monolis_test_check_eq_C1("monolis_set_scalar_to_sparse_matrix_C_test 4b", monolis%MAT%C%A(14), (4.0d0, 4.0d0))
    call monolis_test_check_eq_C1("monolis_set_scalar_to_sparse_matrix_C_test 4c", monolis%MAT%C%A(15), (6.0d0, 6.0d0))
    call monolis_test_check_eq_C1("monolis_set_scalar_to_sparse_matrix_C_test 4d", monolis%MAT%C%A(16), (8.0d0, 8.0d0))

    bval(1,1) = (1.0d0, 1.0d0)
    bval(1,2) = (2.0d0, 2.0d0)
    bval(2,1) = (3.0d0, 3.0d0)
    bval(2,2) = (4.0d0, 4.0d0)

    e(1) = 2
    f(1) = 3

    call monolis_add_matrix_to_sparse_matrix_offdiag_C(monolis, 1, 1, e, f, bval)

    call monolis_test_check_eq_C1("monolis_set_scalar_to_sparse_matrix_C_test 5a", monolis%MAT%C%A(17), (1.0d0, 1.0d0))
    call monolis_test_check_eq_C1("monolis_set_scalar_to_sparse_matrix_C_test 5b", monolis%MAT%C%A(18), (2.0d0, 2.0d0))
    call monolis_test_check_eq_C1("monolis_set_scalar_to_sparse_matrix_C_test 5c", monolis%MAT%C%A(19), (3.0d0, 3.0d0))
    call monolis_test_check_eq_C1("monolis_set_scalar_to_sparse_matrix_C_test 5d", monolis%MAT%C%A(20), (4.0d0, 4.0d0))

    monolis%MAT%C%A = (2.0d0, 2.0d0)

    B = (0.0d0, 0.0d0)

    call monolis_set_Dirichlet_bc_C(monolis, B, 2, 1, (1.0d0, 0.0d0))

    call monolis_test_check_eq_C1("monolis_set_scalar_to_sparse_matrix_C_test 6a", B(1), (-2.0d0, -2.0d0))
    call monolis_test_check_eq_C1("monolis_set_scalar_to_sparse_matrix_C_test 6b", B(2), (-2.0d0, -2.0d0))
    call monolis_test_check_eq_C1("monolis_set_scalar_to_sparse_matrix_C_test 6c", B(3), ( 1.0d0,  0.0d0))
    call monolis_test_check_eq_C1("monolis_set_scalar_to_sparse_matrix_C_test 6d", B(4), (-2.0d0, -2.0d0))
    call monolis_test_check_eq_C1("monolis_set_scalar_to_sparse_matrix_C_test 6e", B(5), (-2.0d0, -2.0d0))
    call monolis_test_check_eq_C1("monolis_set_scalar_to_sparse_matrix_C_test 6f", B(6), (-2.0d0, -2.0d0))
    call monolis_test_check_eq_C1("monolis_set_scalar_to_sparse_matrix_C_test 6g", B(7), ( 0.0d0,  0.0d0))
    call monolis_test_check_eq_C1("monolis_set_scalar_to_sparse_matrix_C_test 6h", B(8), ( 0.0d0,  0.0d0))
  end subroutine monolis_set_scalar_to_sparse_matrix_C_test

end module mod_monolis_spmat_handler_test
