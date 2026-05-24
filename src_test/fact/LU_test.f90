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

    call monolis_precond_LU_nn_block_test(2)
    call monolis_precond_LU_nn_block_test(3)
  end subroutine monolis_precond_LU_nn_test

  !> NDOF ブロックケースのテスト（NDOF >= 2）
  subroutine monolis_precond_LU_nn_block_test(ndof)
    implicit none
    integer(kint), intent(in) :: ndof
    type(monolis_structure) :: mat
    type(monolis_structure) :: prec
    type(monolis_com) :: com
    integer(kint) :: n_node, nelem, elem(2,4), i, j, k, nb
    real(kdouble), allocatable :: x_true(:), b(:), x(:)
    real(kdouble) :: diag_val, off_val
    character(len=64) :: tag

    call monolis_initialize(mat)
    call monolis_initialize(prec)
    call monolis_com_initialize_by_self(com)

    n_node = 5
    nelem  = 4
    do i = 1, nelem
      elem(1,i) = i
      elem(2,i) = i + 1
    end do

    call monolis_get_nonzero_pattern_by_simple_mesh_R(mat, n_node, 2, ndof, nelem, elem)

    !> ブロック対角：対角 NDOF*2 + (j+k)、非対角 0.1
    do i = 1, n_node
      do j = 1, ndof
        do k = 1, ndof
          if (j == k) then
            diag_val = dble(ndof) * 2.0d0 + dble(j)
          else
            diag_val = 0.1d0
          end if
          call monolis_add_scalar_to_sparse_matrix_R(mat, i, i, j, k, diag_val)
        end do
      end do
    end do

    !> オフダイアゴナルブロック：-1.0/(j+k) 程度の値
    do i = 1, nelem
      do j = 1, ndof
        do k = 1, ndof
          off_val = -1.0d0 / dble(j + k)
          call monolis_add_scalar_to_sparse_matrix_R(mat, elem(1,i), elem(2,i), j, k, off_val)
          call monolis_add_scalar_to_sparse_matrix_R(mat, elem(2,i), elem(1,i), j, k, off_val)
        end do
      end do
    end do

    nb = n_node * ndof
    call monolis_alloc_R_1d(x_true, nb)
    call monolis_alloc_R_1d(b, nb)
    call monolis_alloc_R_1d(x, nb)

    call monolis_fact_LU_nn_setup_R(mat%MAT, prec%MAT)

    do i = 1, nb
      x_true(i) = dble(i) * 0.5d0 - 1.0d0
    end do
    call monolis_matvec_product_R(mat, com, x_true, b)
    call monolis_fact_LU_nn_apply_R(mat%MAT, prec%MAT, b, x)
    write(tag,'(a,i0,a)') "monolis_precond_LU_nn_block_test ndof=", ndof, " 1"
    call monolis_test_check_eq_R(trim(tag), x, x_true)

    do i = 1, nb
      x_true(i) = sin(dble(i))
    end do
    call monolis_matvec_product_R(mat, com, x_true, b)
    call monolis_fact_LU_nn_apply_R(mat%MAT, prec%MAT, b, x)
    write(tag,'(a,i0,a)') "monolis_precond_LU_nn_block_test ndof=", ndof, " 2"
    call monolis_test_check_eq_R(trim(tag), x, x_true)

    call monolis_fact_LU_nn_clear_R(prec%MAT)
    call monolis_dealloc_R_1d(x_true)
    call monolis_dealloc_R_1d(b)
    call monolis_dealloc_R_1d(x)
    call monolis_finalize(prec)
    call monolis_finalize(mat)
  end subroutine monolis_precond_LU_nn_block_test

end module mod_monolis_precond_LU_nn_test
