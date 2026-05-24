!> 多重フロント法 LU 前進・後退代入のテストモジュール
module mod_monolis_fact_solve_test
  use mod_monolis
  use mod_monolis_fact_analysis
  use mod_monolis_fact_factorize
  use mod_monolis_fact_solve

  implicit none

contains

  subroutine monolis_fact_solve_test
    implicit none
    type(monolis_structure) :: mat
    type(monolis_com) :: com
    type(monolis_mat_lu) :: lu
    integer(kint) :: n_node, nelem, elem(2,4), i
    real(kdouble), allocatable :: x_true(:), b(:)

    call monolis_std_global_log_string("monolis_fact_solve")

    call monolis_initialize(mat)
    call monolis_com_initialize_by_self(com)

    n_node = 5
    nelem  = 4
    do i = 1, nelem
      elem(1,i) = i
      elem(2,i) = i + 1
    end do

    call monolis_get_nonzero_pattern_by_simple_mesh_R(mat, n_node, 2, 1, nelem, elem)

    do i = 1, n_node
      call monolis_add_scalar_to_sparse_matrix_R(mat, i, i, 1, 1, 4.0d0)
    end do
    do i = 1, nelem
      call monolis_add_scalar_to_sparse_matrix_R(mat, elem(1,i), elem(2,i), 1, 1, -1.0d0)
      call monolis_add_scalar_to_sparse_matrix_R(mat, elem(2,i), elem(1,i), 1, 1, -1.0d0)
    end do

    call monolis_alloc_R_1d(x_true, n_node)
    call monolis_alloc_R_1d(b, n_node)
    do i = 1, n_node
      x_true(i) = dble(i)
    end do

    call monolis_matvec_product_R(mat, com, x_true, b)

    call monolis_fact_analysis(mat%MAT, lu)
    call monolis_fact_factorize(mat%MAT, lu)
    call monolis_fact_solve(lu, b)

    call monolis_test_check_eq_R("monolis_fact_solve_test", b, x_true)

    call monolis_dealloc_R_1d(x_true)
    call monolis_dealloc_R_1d(b)
    call monolis_mat_finalize_LU(lu)
    call monolis_finalize(mat)
  end subroutine monolis_fact_solve_test

end module mod_monolis_fact_solve_test
