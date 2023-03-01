!> 疎行列操作関数群テストモジュール
module mod_monolis_spmat_handler_util_test
  use mod_monolis

  implicit none

contains

  subroutine monolis_spmat_handler_util_test
    implicit none

    !call monolis_set_scalar_to_sparse_matrix_main_R(index, item, A, ndof, ci, cj, csub_i, csub_j, val)
    !call monolis_set_scalar_to_sparse_matrix_main_C(index, item, A, ndof, ci, cj, csub_i, csub_j, val)
    !call monolis_set_block_to_sparse_matrix_main_R(index, item, A, ndof, ci, cj, val)
    !call monolis_set_block_to_sparse_matrix_main_C(index, item, A, ndof, ci, cj, val)
    !call monolis_get_scalar_from_sparse_matrix_main_R(index, item, A, ndof, ci, cj, csub_i, csub_j, val, is_find)
    !call monolis_get_scalar_from_sparse_matrix_main_C(index, item, A, ndof, ci, cj, csub_i, csub_j, val, is_find)
    !call monolis_add_scalar_to_sparse_matrix_main_R(index, item, A, ndof, ci, cj, csub_i, csub_j, val)
    !call monolis_add_scalar_to_sparse_matrix_main_C(index, item, A, ndof, ci, cj, csub_i, csub_j, val)
    !call monolis_add_matrix_to_sparse_matrix_main_R(index, item, A, n1, n2, ndof, e1, e2, val)
    !call monolis_add_matrix_to_sparse_matrix_main_C(index, item, A, n1, n2, ndof, e1, e2, val)
    !call monolis_set_Dirichlet_bc_main_R(index, item, A, B, indexR, itemR, permA, &
    !call monolis_set_Dirichlet_bc_main_C(index, item, A, B, indexR, itemR, permA, &
    !call monolis_stop_by_matrix_assemble(ci, cj)
    !call monolis_stop_by_submatrix_access(ndof, sub_dof)
    !call monolis_get_max_matrix_component_main_R(monoMAT, monoCOM, max_val)
    !call monolis_check_diagonal_zero_component_main_R(monoPRM, monoMAT)
  end subroutine monolis_spmat_handler_util_test

end module mod_monolis_spmat_handler_util_test
