!> LU 前処理（nxn ブロック）テストモジュール
module mod_monolis_precond_LU_nn_test
  use mod_monolis
  use mod_monolis_fact_LU_nn

  implicit none

contains

  subroutine monolis_precond_LU_nn_test
    implicit none
    type(monolis_structure) :: mat
    type(monolis_structure) :: prec
    type(monolis_com) :: com
    integer(kint) :: n_node, nelem, elem(2,4), i
    real(kdouble), allocatable :: x_true(:), b(:), x(:)

    call monolis_std_global_log_string("monolis_fact_LU_nn_setup_R")
    call monolis_std_global_log_string("monolis_fact_LU_nn_apply_R")
    call monolis_std_global_log_string("monolis_fact_LU_nn_clear_R")

    call monolis_initialize(mat)
    call monolis_initialize(prec)
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
    call monolis_alloc_R_1d(x, n_node)

    call monolis_fact_LU_nn_setup_R(mat%MAT, prec%MAT)

    !> 1 回目: x_true = [1, 2, 3, 4, 5]
    do i = 1, n_node
      x_true(i) = dble(i)
    end do
    call monolis_matvec_product_R(mat, com, x_true, b)
    call monolis_fact_LU_nn_apply_R(mat%MAT, prec%MAT, b, x)
    call monolis_test_check_eq_R("monolis_precond_LU_nn_test 1", x, x_true)

    !> 2 回目: 右辺を変えて、同じ LU データで解けることを確認
    do i = 1, n_node
      x_true(i) = dble(n_node - i + 1)
    end do
    call monolis_matvec_product_R(mat, com, x_true, b)
    call monolis_fact_LU_nn_apply_R(mat%MAT, prec%MAT, b, x)
    call monolis_test_check_eq_R("monolis_precond_LU_nn_test 2", x, x_true)

    !> 3 回目: さらに別の右辺
    do i = 1, n_node
      x_true(i) = sin(dble(i))
    end do
    call monolis_matvec_product_R(mat, com, x_true, b)
    call monolis_fact_LU_nn_apply_R(mat%MAT, prec%MAT, b, x)
    call monolis_test_check_eq_R("monolis_precond_LU_nn_test 3", x, x_true)

    call monolis_fact_LU_nn_clear_R(prec%MAT)

    call monolis_dealloc_R_1d(x_true)
    call monolis_dealloc_R_1d(b)
    call monolis_dealloc_R_1d(x)
    call monolis_finalize(prec)
    call monolis_finalize(mat)
  end subroutine monolis_precond_LU_nn_test

end module mod_monolis_precond_LU_nn_test
