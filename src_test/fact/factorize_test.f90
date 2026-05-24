!> 多重フロント法 LU 数値分解フェーズのテストモジュール
module mod_monolis_fact_factorize_test
  use mod_monolis
  use mod_monolis_fact_analysis
  use mod_monolis_fact_factorize

  implicit none

contains

  subroutine monolis_fact_factorize_test
    implicit none
    type(monolis_structure) :: mat
    type(monolis_com) :: com
    type(monolis_mat_lu) :: lu
    integer(kint) :: n_node, nelem, elem(2,4), i

    call monolis_std_global_log_string("monolis_fact_factorize")

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

    call monolis_fact_analysis(mat%MAT, lu)
    call monolis_fact_factorize(mat%MAT, lu)

    if (lu%factorized) then
      call monolis_test_assert_pass("monolis_fact_factorize_test factorized")
    else
      call monolis_test_assert_fail("monolis_fact_factorize_test factorized", "lu%factorized is false")
    end if

    call monolis_mat_finalize_LU(lu)
    call monolis_finalize(mat)
  end subroutine monolis_fact_factorize_test

end module mod_monolis_fact_factorize_test
