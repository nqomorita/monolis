!> 多重フロント法 LU 解析フェーズのテストモジュール
module mod_monolis_fact_analysis_test
  use mod_monolis
  use mod_monolis_fact_analysis

  implicit none

contains

  subroutine monolis_fact_analysis_test
    implicit none
    type(monolis_structure) :: mat
    type(monolis_com) :: com
    type(monolis_mat_lu) :: lu
    integer(kint) :: n_node, nelem, elem(2,4), i

    call monolis_std_global_log_string("monolis_fact_analysis")

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

    call monolis_test_check_eq_I1("monolis_fact_analysis_test N", lu%N, n_node)
    call monolis_fact_test_check("monolis_fact_analysis_test analyzed", lu%analyzed)
    call monolis_fact_test_check("monolis_fact_analysis_test nfronts>0", lu%nfronts > 0)
    call monolis_fact_test_check("monolis_fact_analysis_test nfronts<=N", lu%nfronts <= n_node)

    call monolis_mat_finalize_LU(lu)
    call monolis_finalize(mat)
  end subroutine monolis_fact_analysis_test

  subroutine monolis_fact_test_check(header, cond)
    implicit none
    character(*), intent(in) :: header
    logical,      intent(in) :: cond
    if (cond) then
      call monolis_test_assert_pass(header)
    else
      call monolis_test_assert_fail(header, "condition false")
    end if
  end subroutine monolis_fact_test_check

end module mod_monolis_fact_analysis_test
