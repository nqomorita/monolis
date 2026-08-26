!> LU 前処理関連テストモジュール
module mod_monolis_precond_LU_test
  use mod_monolis
  use mod_monolis_precond_LU

  implicit none

contains

  subroutine monolis_precond_LU_test()
    implicit none
    type(monolis_structure) :: mat
    type(monolis_structure) :: prec
    type(monolis_com) :: com
    integer(kint) :: n_node, nelem, elem(2,4), i
    real(kdouble), allocatable :: x_true(:), b(:), x(:)

    call monolis_std_global_log_string("monolis_precond_LU_setup_R")
    call monolis_std_global_log_string("monolis_precond_LU_apply_R")
    call monolis_std_global_log_string("monolis_precond_LU_clear_R")

    call monolis_std_global_log_string("monolis_precond_LU_setup_C")
    call monolis_std_global_log_string("monolis_precond_LU_apply_C")
    call monolis_std_global_log_string("monolis_precond_LU_clear_C")

    call monolis_initialize(mat)
    call monolis_initialize(prec)
    call monolis_com_initialize_by_self(com)

    n_node = 5
    nelem = 4
    do i = 1, nelem
      elem(1,i) = i
      elem(2,i) = i + 1
    enddo

    call monolis_get_nonzero_pattern_by_simple_mesh_R(mat, n_node, 2, 1, nelem, elem)

    do i = 1, n_node
      call monolis_add_scalar_to_sparse_matrix_R(mat, i, i, 1, 1, 4.0d0)
    enddo
    do i = 1, nelem
      call monolis_add_scalar_to_sparse_matrix_R(mat, elem(1,i), elem(2,i), 1, 1, -1.0d0)
      call monolis_add_scalar_to_sparse_matrix_R(mat, elem(2,i), elem(1,i), 1, 1, -1.0d0)
    enddo

    call monolis_alloc_R_1d(x_true, n_node)
    call monolis_alloc_R_1d(b, n_node)
    call monolis_alloc_R_1d(x, n_node)

    do i = 1, n_node
      x_true(i) = dble(i)
    enddo
    call monolis_matvec_product_R(mat, com, x_true, b)

    mat%PRM%Iarray(monolis_prm_I_is_prec_stored) = monolis_I_true
    call monolis_precond_LU_setup_R(mat%PRM, com, mat%MAT, prec%MAT)

    !> 保存済み LU があれば、行列値が変わっても再分解しないことを確認
    call monolis_add_scalar_to_sparse_matrix_R(mat, 1, 1, 1, 1, 1.0d0)
    call monolis_precond_LU_setup_R(mat%PRM, com, mat%MAT, prec%MAT)
    call monolis_precond_LU_clear_R(mat%PRM, com, mat%MAT, prec%MAT)

    if(prec%MAT%LU%factorized)then
      call monolis_test_assert_pass("monolis_precond_LU_test stored factor")
      call monolis_precond_LU_apply_R(mat%PRM, com, mat%MAT, prec%MAT, b, x)
      call monolis_test_check_eq_R("monolis_precond_LU_test reused factor", x, x_true)
    else
      call monolis_test_assert_fail("monolis_precond_LU_test stored factor", &
        & "LU factor was cleared despite is_prec_stored")
    endif

    !> 保存指定を解除した後は、従来どおり LU データを解放することを確認
    mat%PRM%Iarray(monolis_prm_I_is_prec_stored) = monolis_I_false
    call monolis_precond_LU_clear_R(mat%PRM, com, mat%MAT, prec%MAT)
    if(.not. prec%MAT%LU%factorized)then
      call monolis_test_assert_pass("monolis_precond_LU_test released factor")
    else
      call monolis_test_assert_fail("monolis_precond_LU_test released factor", &
        & "LU factor remains after is_prec_stored is disabled")
    endif

    call monolis_dealloc_R_1d(x_true)
    call monolis_dealloc_R_1d(b)
    call monolis_dealloc_R_1d(x)
    call monolis_finalize(prec)
    call monolis_finalize(mat)
  end subroutine monolis_precond_LU_test
end module mod_monolis_precond_LU_test
