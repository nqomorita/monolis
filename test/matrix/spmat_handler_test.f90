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
    !call monolis_set_scalar_to_sparse_matrix_C(monolis, i, j, sub_i, sub_j, val)
    !call monolis_set_block_to_sparse_matrix_R(monolis, i, j, val)
    !call monolis_set_block_to_sparse_matrix_C(monolis, i, j, val)
    !call monolis_get_scalar_from_sparse_matrix_R(monolis, i, j, sub_i, sub_j, val, is_find)
    !call monolis_get_scalar_from_sparse_matrix_C(monolis, i, j, sub_i, sub_j, val, is_find)
    !call monolis_add_scalar_to_sparse_matrix_R(monolis, i, j, sub_i, sub_j, val)
    !call monolis_add_scalar_to_sparse_matrix_C(monolis, i, j, sub_i, sub_j, val)
    !call monolis_add_matrix_to_sparse_matrix_R(monolis, n_base, connectivity, mat)
    !call monolis_add_matrix_to_sparse_matrix_C(monolis, n_base, connectivity, mat)
    !call monolis_add_matrix_to_sparse_matrix_offdiag_R(monolis, n_base_i, n_base_j, &
    !call monolis_add_matrix_to_sparse_matrix_offdiag_C(monolis, n_base_i, n_base_j, &
    !call monolis_set_matrix_BCSR_R(monolis, N, NP, NDOF, NZ, A, index, item)
    !call monolis_set_matrix_BCSR_mat_val_R(monolis, NDOF, NZ, A)
    !call monolis_set_Dirichlet_bc_R(monolis, B, node_id, ndof_bc, val)
    !call monolis_set_Dirichlet_bc_C(monolis, B, node_id, ndof_bc, val)
  end subroutine monolis_spmat_handler_test

  subroutine monolis_set_scalar_to_sparse_matrix_R_test()
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> 行番号
    integer(kint) :: i
    !> 列番号
    integer(kint) :: j
    !> ブロック中の行番号
    integer(kint) :: sub_i
    !> ブロック中の列番号
    integer(kint) :: sub_j
    !> 行列値
    real(kdouble) :: val
    !> 節点数
    integer(kint) :: n_node
    !> 要素を構成する節点数
    integer(kint) :: n_base
    !> 自由度数
    integer(kint) :: ndof
    !> 要素数
    integer(kint) :: n_elem
    !> 要素
    integer(kint) :: elem(2,3)

    call monolis_std_log_string("monolis_set_scalar_to_sparse_matrix_R_test")

    n_node = 4

    n_base = 2

    ndof = 2

    n_elem = 3

    elem(1,1) = 1; elem(2,1) = 2
    elem(1,2) = 2; elem(2,2) = 3
    elem(1,3) = 3; elem(2,3) = 4

    call monolis_get_nonzero_pattern_by_simple_mesh_R(monolis, n_node, n_base, ndof, n_elem, elem)

    val = 1.0d0

    call monolis_set_scalar_to_sparse_matrix_R(monolis, 2, 2, 1, 1, val)

    call monolis_test_check_eq_R1("monolis_set_scalar_to_sparse_matrix_R_test 1", monolis%MAT%R%A(13), 1.0d0)
  end subroutine monolis_set_scalar_to_sparse_matrix_R_test

end module mod_monolis_spmat_handler_test
