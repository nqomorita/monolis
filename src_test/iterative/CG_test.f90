!> CG 法テストモジュール
module mod_monolis_solver_CG_test
  use mod_monolis

  implicit none

contains

  subroutine monolis_solver_CG_test()
    implicit none
    integer(kint) :: n_dof

    do n_dof = 1, 1
      !call monolis_solver_CG_test_main(n_dof, monolis_prec_NONE)
      !call monolis_solver_CG_test_main(n_dof, monolis_prec_DIAG)
      !call monolis_solver_CG_test_main(n_dof, monolis_prec_SOR)
      !call monolis_solver_CG_test_main(n_dof, monolis_prec_AMG)
      call monolis_solver_CG_test_main(n_dof, monolis_prec_LU)
    enddo
  end subroutine monolis_solver_CG_test

  subroutine monolis_solver_CG_test_main(n_dof, prec)
    implicit none
    type(monolis_structure) :: mat
    type(monolis_com) :: com
    integer(kint) :: n_node, nelem
    integer(kint) :: i1, i2, j1, j2
    integer(kint) :: n_dof, prec
    real(kdouble) :: val
    integer(kint), allocatable :: elem(:,:)
    real(kdouble), allocatable :: a(:), b(:)

    call monolis_std_global_log_string("monolis_solver_CG")
    call monolis_std_log_I1("DOF", n_dof)
    call monolis_std_log_I1("PRECOND", prec)

    call monolis_initialize(mat)
    call monolis_com_initialize_by_self(com)

    n_node = 3
    nelem = n_node - 1

    call monolis_alloc_I_2d(elem, 2, nelem)
    call monolis_alloc_R_1d(a, n_dof*n_node)
    call monolis_alloc_R_1d(b, n_dof*n_node)

    do i1 = 1, nelem
      elem(1,i1) = i1
      elem(2,i1) = i1 + 1
    enddo

    call monolis_get_nonzero_pattern_by_simple_mesh_R(mat, n_node, 2, n_dof, nelem, elem)

    do i1 = 1, n_node
      do i2 = 1, n_dof
        !call random_number(val)
        val = 4.0d0
        !val = val + 2.0d0
        call monolis_add_scalar_to_sparse_matrix_R(mat, i1, i1, i2, i2, val)
      enddo
    enddo

    do i1 = 1, nelem
      do i2 = 1, n_dof
      do j2 = 1, n_dof
        !call random_number(val)
        val = 2.0d0
        call monolis_add_scalar_to_sparse_matrix_R(mat, elem(1,i1), elem(2,i1), i2, j2, val)
        call monolis_add_scalar_to_sparse_matrix_R(mat, elem(2,i1), elem(1,i1), j2, i2, val)
      enddo
      enddo
    enddo

    a = 1.0d0

    call monolis_matvec_product_R(mat, com, a, b)

    call monolis_set_method(mat, monolis_iter_CG)
    call monolis_set_precond(mat, prec)
    call monolis_set_maxiter(mat, 30)
    call monolis_set_tolerance(mat, 1.0d-10)
    call monolis_show_timelog_statistics(mat, .true.)

    a = 0.0d0

    call monolis_solve_R(mat, com, b, a)

    b = 1.0d0

    call monolis_test_check_eq_R("monolis_solver_CG_test_main", a, b)

    call monolis_finalize(mat)
  end subroutine monolis_solver_CG_test_main

end module mod_monolis_solver_CG_test
