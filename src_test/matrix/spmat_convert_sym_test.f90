!> 疎行列操作関数群
module mod_monolis_spmat_convert_sym_test
  use mod_monolis
  use mod_monolis_spmat_convert_sym
  use mod_monolis_spmat_nonzero_pattern

  implicit none

contains

  subroutine monolis_spmat_convert_sym_test()
    implicit none
    type(monolis_structure) :: mat
    type(monolis_com) :: com
    integer(kint) :: n_node, n_elem, elem(2,4)

    call monolis_initialize(mat)
    call monolis_com_initialize_by_self(com)

    n_node = 5
    n_elem = 4

    if(monolis_mpi_get_global_comm_size() == 2) return

    elem(1,1) = 1; elem(2,1) = 2;
    elem(1,2) = 2; elem(2,2) = 3;
    elem(1,3) = 3; elem(2,3) = 4;
    elem(1,4) = 4; elem(2,4) = 5;

    call monolis_get_nonzero_pattern_by_simple_mesh_R(mat, n_node, 2, 1, n_elem, elem)

    call monolis_add_scalar_to_sparse_matrix_R(mat, 1, 1, 1, 1, 2.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 1, 2, 1, 1, 1.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 2, 1, 1, 1, 1.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 2, 2, 1, 1, 2.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 2, 3, 1, 1, 3.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 3, 2, 1, 1, 1.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 3, 3, 1, 1, 2.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 3, 4, 1, 1, 4.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 4, 3, 1, 1, 1.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 4, 4, 1, 1, 2.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 4, 5, 1, 1, 5.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 5, 4, 1, 1, 1.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 5, 5, 1, 1, 2.0d0)

    call monolis_matrix_convert_to_symmetric_inner_R(mat%MAT, com)

    call monolis_matrix_convert_to_symmetric_outer_R(mat%MAT, com)

    call monolis_test_check_eq_R1("monolis_spmat_convert_sym_test 1", mat%MAT%R%A( 1), 2.0d0)
    call monolis_test_check_eq_R1("monolis_spmat_convert_sym_test 1", mat%MAT%R%A( 2), 1.0d0)
    call monolis_test_check_eq_R1("monolis_spmat_convert_sym_test 1", mat%MAT%R%A( 3), 1.0d0)
    call monolis_test_check_eq_R1("monolis_spmat_convert_sym_test 1", mat%MAT%R%A( 4), 2.0d0)
    call monolis_test_check_eq_R1("monolis_spmat_convert_sym_test 1", mat%MAT%R%A( 5), 2.0d0)
    call monolis_test_check_eq_R1("monolis_spmat_convert_sym_test 1", mat%MAT%R%A( 6), 2.0d0)
    call monolis_test_check_eq_R1("monolis_spmat_convert_sym_test 1", mat%MAT%R%A( 7), 2.0d0)
    call monolis_test_check_eq_R1("monolis_spmat_convert_sym_test 1", mat%MAT%R%A( 8), 2.5d0)
    call monolis_test_check_eq_R1("monolis_spmat_convert_sym_test 1", mat%MAT%R%A( 9), 2.5d0)
    call monolis_test_check_eq_R1("monolis_spmat_convert_sym_test 1", mat%MAT%R%A(10), 2.0d0)
    call monolis_test_check_eq_R1("monolis_spmat_convert_sym_test 1", mat%MAT%R%A(11), 3.0d0)
    call monolis_test_check_eq_R1("monolis_spmat_convert_sym_test 1", mat%MAT%R%A(12), 3.0d0)
    call monolis_test_check_eq_R1("monolis_spmat_convert_sym_test 1", mat%MAT%R%A(13), 2.0d0)
  end subroutine monolis_spmat_convert_sym_test

end module mod_monolis_spmat_convert_sym_test
