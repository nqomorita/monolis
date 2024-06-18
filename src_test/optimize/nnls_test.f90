!> NNLS 関数
module mod_monolis_opt_nnls_test
  use mod_monolis
  implicit none

contains

  !> @ingroup linalg
  !> Non-Negative Least Squares
  subroutine monolis_optimize_nnls_test()
    implicit none
    integer(kint) :: max_iter
    real(kdouble) :: tol, residual
    real(kdouble) :: A(4,2)
    real(kdouble) :: b(4)
    real(kdouble) :: x(2)

    call monolis_std_global_log_string("monolis_optimize_nnls")
    
    max_iter = 100
    tol = 1.0d-6

    A(1,1) = 1.0d0; A(1,2) = 1.0d0; 
    A(2,1) = 1.0d0; A(2,2) = 1.0d0; 
    A(3,1) = 1.0d0; A(3,2) = 1.0d0; 
    A(4,1) = 1.0d0; A(4,2) = 1.0d0; 

    b(1) = 1.0d0
    b(2) = 1.0d0
    b(3) = 1.0d0
    b(4) = 1.0d0

    call monolis_optimize_nnls(A, b, x, 4, 2, max_iter, tol, residual)

    call monolis_test_check_eq_R1("monolis_optimize_nnls_test 1a", x(1), 1.0d0)
    call monolis_test_check_eq_R1("monolis_optimize_nnls_test 1b", x(2), 0.0d0)

    A = 0.0d0
    A(1,1) = 1.0d0; A(1,2) = 0.0d0; 
    A(2,1) = 1.0d0; A(2,2) = 0.0d0; 
    A(3,1) = 0.0d0; A(3,2) = 1.0d0; 

    b(1) = 2.0d0
    b(2) = 1.0d0
    b(3) = 1.0d0

    call monolis_optimize_nnls(A(1:3,1:2), b(1:3), x(1:2), 3, 2, max_iter, tol, residual)

    call monolis_test_check_eq_R1("monolis_optimize_nnls_test 2a", x(1), 1.5d0)
    call monolis_test_check_eq_R1("monolis_optimize_nnls_test 2b", x(2), 1.0d0)

    A(1,1) = 1.0d0; A(1,2) = 0.0d0; 
    A(2,1) = 1.0d0; A(2,2) = 0.0d0; 
    A(3,1) = 0.0d0; A(3,2) = 1.0d0; 

    b(1) =-1.0d0
    b(2) =-1.0d0
    b(3) =-1.0d0

    call monolis_optimize_nnls(A(1:3,1:2), b(1:3), x(1:2), 3, 2, max_iter, tol, residual)

    call monolis_test_check_eq_R1("monolis_optimize_nnls_test 3a", x(1), 0.0d0)
    call monolis_test_check_eq_R1("monolis_optimize_nnls_test 3b", x(2), 0.0d0)
  end subroutine monolis_optimize_nnls_test

end module mod_monolis_opt_nnls_test