!> CG 法テストモジュール
module mod_monolis_solver_CG_test
  use mod_monolis

  implicit none

contains

  subroutine monolis_solver_CG_test()
    implicit none

    call monolis_solver_CG_33_test(monolis_prec_NONE)
    call monolis_solver_CG_33_test(monolis_prec_DIAG)
    call monolis_solver_CG_33_test(monolis_prec_SOR)
    !call monolis_solver_CG_33_test(monolis_prec_MUMPS)
  end subroutine monolis_solver_CG_test

  subroutine monolis_solver_CG_33_test(prec)
    implicit none
    type(monolis_structure) :: mat
    integer(kint) :: nnode, nelem, elem(2,9)
    integer(kint) :: i1, i2, j1, j2
    integer(kint) :: prec
    real(kdouble) :: val
    real(kdouble) :: a(30), b(30)

    call monolis_std_log_string("monolis_solver_CG_33_test")

    call monolis_initialize(mat, "./")

    nnode = 10

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

    call monolis_get_nonzero_pattern_by_simple_mesh_R(mat, nnode, 2, 3, nelem, elem)

    do i1 = 1, 10
      do i2 = 1, 3
        call random_number(val)
        val = val + 1.0d0
        call monolis_add_scalar_to_sparse_matrix_R(mat, i1, i1, i2, i2, val)
      enddo
    enddo

    do i1 = 1, 9
      do i2 = 1, 3
      do j2 = 1, 3
        call random_number(val)
        call monolis_add_scalar_to_sparse_matrix_R(mat, elem(1,i1), elem(2,i1), i2, j2, val)
        call monolis_add_scalar_to_sparse_matrix_R(mat, elem(2,i1), elem(1,i1), j2, i2, val)
      enddo
      enddo
    enddo

    a = 1.0d0

    call monolis_matvec_product_R(mat, a, b)

    call monolis_set_precond(mat, prec)
    call monolis_set_tolerance(mat, 1.0d-10)

    call monolis_solve(mat, b, a)

    b = 1.0d0

    call monolis_test_check_eq_R("monolis_solver_CG_33_test", a, b)

    call monolis_finalize(mat)
  end subroutine monolis_solver_CG_33_test
end module mod_monolis_solver_CG_test
