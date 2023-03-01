!> 疎行列操作関数群
module mod_monolis_spmat_handler_test
  use mod_monolis

  implicit none

contains

  subroutine monolis_spmat_handler_test
    implicit none

    !call monolis_set_scalar_to_sparse_matrix_R(monolis, i, j, sub_i, sub_j, val)
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

end module mod_monolis_spmat_handler_test
