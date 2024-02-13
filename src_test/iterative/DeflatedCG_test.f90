!> DeflatedCG 法テストモジュール
module mod_monolis_solver_DeflatedCG_test
  use mod_monolis

  implicit none

contains

  subroutine monolis_solver_DeflatedCG_test()
    implicit none
    integer(kint) :: n_dof

    do n_dof = 1, 3
      call monolis_solver_DeflatedCG_test_main(n_dof, monolis_prec_NONE)
      !call monolis_solver_DeflatedCG_test_main(n_dof, monolis_prec_DIAG)
      !call monolis_solver_DeflatedCG_test_main(n_dof, monolis_prec_SOR)
    enddo
  end subroutine monolis_solver_DeflatedCG_test

  subroutine monolis_solver_DeflatedCG_test_main(n_dof, prec)
    implicit none
    type(monolis_structure) :: mat
    type(monolis_com) :: com
    integer(kint) :: n_node, nelem, elem(2,9)
    integer(kint) :: i1, i2, j1, j2
    integer(kint) :: n_dof, prec
    integer(kint) :: n_get_eigen
    real(kdouble) :: val
    real(kdouble) :: a(n_dof*10), b(n_dof*10)
    real(kdouble) :: eig_val(n_dof*10), eig_mode(n_dof*10,4*n_dof)
    logical :: is_bc(n_dof*10)

    call monolis_std_global_log_string("monolis_solver_DeflatedCG")
    call monolis_std_log_I1("DOF", n_dof)
    call monolis_std_log_I1("PRECOND", prec)

    call monolis_initialize(mat)
    call monolis_com_initialize_by_self(com)

    n_node = 10

    nelem = 9

    elem(1,1) = 1; elem(2,1) = 2;
    elem(1,2) = 2; elem(2,2) = 3;
    elem(1,3) = 3; elem(2,3) = 4;
    elem(1,4) = 4; elem(2,4) = 5;
    elem(1,5) = 5; elem(2,5) = 6;
    elem(1,6) = 6; elem(2,6) = 7;
    elem(1,7) = 7; elem(2,7) = 8;
    elem(1,8) = 8; elem(2,8) = 9;
    elem(1,9) = 9; elem(2,9) =10;

    call monolis_get_nonzero_pattern_by_simple_mesh_R(mat, n_node, 2, n_dof, nelem, elem)

    do i1 = 1, 10
      do i2 = 1, n_dof
        call random_number(val)
        val = val + 2.0d0
        call monolis_add_scalar_to_sparse_matrix_R(mat, i1, i1, i2, i2, val)
      enddo
    enddo

    do i1 = 1, 9
      do i2 = 1, n_dof
      do j2 = 1, n_dof
        call random_number(val)
        call monolis_add_scalar_to_sparse_matrix_R(mat, elem(1,i1), elem(2,i1), i2, j2, val)
        call monolis_add_scalar_to_sparse_matrix_R(mat, elem(2,i1), elem(1,i1), j2, i2, val)
      enddo
      enddo
    enddo

    a = 1.0d0

    call monolis_matvec_product_R(mat, com, a, b)

    !> get deflation mode
    call monolis_std_global_log_string("get deflation mode")
    call monolis_set_method(mat, monolis_iter_CG)
    call monolis_set_precond(mat, monolis_prec_NONE)
    call monolis_show_iterlog(mat, .false.)
    call monolis_show_timelog(mat, .false.)
    call monolis_show_summary(mat, .false.)
    call monolis_show_timelog_statistics(mat, .false.)

    n_get_eigen = 4*n_dof
    is_bc = .false.
    call monolis_eigen_inverted_standard_lanczos_R &
      & (mat, com, n_get_eigen, 1.0d-10, n_dof*20, eig_val, eig_mode, is_bc)

    call monolis_set_deflation_mode(mat, n_get_eigen, eig_mode)

    !> solve
    !> monolis_iter_DeflatedCG1
    call monolis_std_global_log_string("Deflated CG main")
    call monolis_set_method(mat, monolis_iter_DeflatedCG1)
    call monolis_set_precond(mat, prec)
    call monolis_set_maxiter(mat, 2000)
    call monolis_set_tolerance(mat, 1.0d-9)
    call monolis_show_iterlog(mat, .true.)
    call monolis_show_timelog(mat, .true.)
    call monolis_show_summary(mat, .true.)
    call monolis_show_timelog_statistics(mat, .true.)

    a = 0.0d0

    call monolis_solve_R(mat, com, b, a)

    b = 1.0d0

    call monolis_test_check_eq_R("monolis_solver_DeflatedCG1_test_main", a, b)

    !> monolis_iter_DeflatedCG2
    call monolis_set_method(mat, monolis_iter_DeflatedCG2)

    a = 0.0d0

    call monolis_solve_R(mat, com, b, a)

    b = 1.0d0

    call monolis_test_check_eq_R("monolis_solver_DeflatedCG2_test_main", a, b)

    !> monolis_iter_ADeflatedCG2
    call monolis_set_method(mat, monolis_iter_ADeflatedCG2)

    a = 0.0d0

    call monolis_solve_R(mat, com, b, a)

    b = 1.0d0

    call monolis_test_check_eq_R("monolis_solver_ADeflatedCG_test_main", a, b)

    call monolis_finalize(mat)
  end subroutine monolis_solver_DeflatedCG_test_main

end module mod_monolis_solver_DeflatedCG_test
