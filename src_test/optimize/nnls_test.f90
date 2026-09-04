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
    real(kdouble) :: tol, tol_inner, residual
    real(kdouble) :: A(4,2)
    real(kdouble) :: b(4)
    real(kdouble) :: x(2)
    integer(kint) :: self_comm

    call monolis_std_global_log_string("monolis_optimize_nnls_R")
    call monolis_std_global_log_string("monolis_optimize_nnls_R_with_sparse_solution")
    call monolis_std_global_log_string("monolis_optimize_nnls_LH_main")

    self_comm = monolis_mpi_get_self_comm()
    max_iter = 100
    tol = 1.0d-6
    tol_inner = 1.0d-12

    A(1,1) = 1.0d0; A(1,2) = 1.0d0; 
    A(2,1) = 1.0d0; A(2,2) = 1.0d0; 
    A(3,1) = 1.0d0; A(3,2) = 1.0d0; 
    A(4,1) = 1.0d0; A(4,2) = 1.0d0; 

    b(1) = 1.0d0
    b(2) = 1.0d0
    b(3) = 1.0d0
    b(4) = 1.0d0

    call monolis_optimize_nnls_R(A, b, x, 4, 2, max_iter, tol, residual, self_comm)

    call monolis_test_check_eq_R1("monolis_optimize_nnls_test 1a", x(1), 1.0d0)
    call monolis_test_check_eq_R1("monolis_optimize_nnls_test 1b", x(2), 0.0d0)

    call monolis_optimize_nnls_R_with_sparse_solution(A, b, x, 4, 2, max_iter, tol, tol_inner, residual, self_comm)

    call monolis_test_check_eq_R1("monolis_optimize_nnls_test 4a", x(1), 1.0d0)
    call monolis_test_check_eq_R1("monolis_optimize_nnls_test 4b", x(2), 0.0d0)

    A = 0.0d0
    A(1,1) = 1.0d0; A(1,2) = 0.0d0; 
    A(2,1) = 1.0d0; A(2,2) = 0.0d0; 
    A(3,1) = 0.0d0; A(3,2) = 1.0d0; 

    b(1) = 2.0d0
    b(2) = 1.0d0
    b(3) = 1.0d0

    call monolis_optimize_nnls_R(A(1:3,1:2), b(1:3), x(1:2), 3, 2, max_iter, tol, residual, self_comm)

    call monolis_test_check_eq_R1("monolis_optimize_nnls_test 2a", x(1), 1.5d0)
    call monolis_test_check_eq_R1("monolis_optimize_nnls_test 2b", x(2), 1.0d0)

    call monolis_optimize_nnls_R_with_sparse_solution(A(1:3,1:2), b(1:3), x(1:2), 3, 2, max_iter, tol, tol_inner, &
      & residual, self_comm)

    call monolis_test_check_eq_R1("monolis_optimize_nnls_test 5a", x(1), 1.5d0)
    call monolis_test_check_eq_R1("monolis_optimize_nnls_test 5b", x(2), 1.0d0)

    A(1,1) = 1.0d0; A(1,2) = 0.0d0; 
    A(2,1) = 1.0d0; A(2,2) = 0.0d0; 
    A(3,1) = 0.0d0; A(3,2) = 1.0d0; 

    b(1) =-1.0d0
    b(2) =-1.0d0
    b(3) =-1.0d0

    call monolis_optimize_nnls_R(A(1:3,1:2), b(1:3), x(1:2), 3, 2, max_iter, tol, residual, self_comm)

    call monolis_test_check_eq_R1("monolis_optimize_nnls_test 3a", x(1), 0.0d0)
    call monolis_test_check_eq_R1("monolis_optimize_nnls_test 3b", x(2), 0.0d0)

    call monolis_optimize_nnls_R_with_sparse_solution(A(1:3,1:2), b(1:3), x(1:2), 3, 2, max_iter, tol, tol_inner, &
      & residual, self_comm)

    call monolis_test_check_eq_R1("monolis_optimize_nnls_test 6a", x(1), 0.0d0)
    call monolis_test_check_eq_R1("monolis_optimize_nnls_test 6b", x(2), 0.0d0)

    call monolis_optimize_parallel_nnls_test()
  end subroutine monolis_optimize_nnls_test

  !> @ingroup linalg
  !> Non-Negative Least Squares（列分散の並列実行とセルフコミュニケータ実行の比較）
  subroutine monolis_optimize_parallel_nnls_test()
    implicit none
    integer(kint) :: max_iter, n_loc, offset, i, j
    real(kdouble) :: tol, tol_inner, residual, residual_ref
    real(kdouble) :: A_g(4,4)
    real(kdouble) :: b(4)
    real(kdouble) :: x_ref(4)
    real(kdouble), allocatable :: A_loc(:,:)
    real(kdouble), allocatable :: x(:)
    integer(kint) :: comm, self_comm

    if(monolis_mpi_get_global_comm_size() > 2) return

    comm = monolis_mpi_get_global_comm()
    self_comm = monolis_mpi_get_self_comm()
    max_iter = 100
    tol = 1.0d-6
    tol_inner = 1.0d-12

    !> 対角優位な 4x4 行列（全成分が active になるケース）
    A_g = 0.1d0
    do i = 1, 4
      A_g(i,i) = 2.0d0
    enddo

    b(1) = 1.0d0
    b(2) = 2.0d0
    b(3) = 3.0d0
    b(4) = 4.0d0

    !> セルフコミュニケータによる参照解
    call monolis_optimize_nnls_R_with_sparse_solution(A_g, b, x_ref, 4, 4, max_iter, tol, tol_inner, &
      & residual_ref, self_comm)

    !> 列ブロック分散（1 ランク: 4 列、2 ランク: 各 2 列）
    if(monolis_mpi_get_global_comm_size() == 1)then
      n_loc = 4
      offset = 0
    else
      n_loc = 2
      offset = 2*monolis_mpi_get_global_my_rank()
    endif

    call monolis_alloc_R_2d(A_loc, 4, n_loc)
    call monolis_alloc_R_1d(x, n_loc)

    do i = 1, n_loc
      do j = 1, 4
        A_loc(j,i) = A_g(j,offset+i)
      enddo
    enddo

    call monolis_optimize_nnls_R_with_sparse_solution(A_loc, b, x, 4, n_loc, max_iter, tol, tol_inner, residual, comm)

    do i = 1, n_loc
      call monolis_test_check_eq_R1("monolis_optimize_parallel_nnls_test 1a", x(i), x_ref(offset+i))
    enddo
    call monolis_test_check_eq_R1("monolis_optimize_parallel_nnls_test 1b", residual, residual_ref)

    !> 疎な解になるケース（b が第 2 列のスカラー倍）
    do i = 1, 4
      b(i) = 2.0d0*A_g(i,2)
    enddo

    call monolis_optimize_nnls_R_with_sparse_solution(A_g, b, x_ref, 4, 4, max_iter, tol, tol_inner, &
      & residual_ref, self_comm)

    x = 0.0d0
    call monolis_optimize_nnls_R_with_sparse_solution(A_loc, b, x, 4, n_loc, max_iter, tol, tol_inner, residual, comm)

    do i = 1, n_loc
      call monolis_test_check_eq_R1("monolis_optimize_parallel_nnls_test 2a", x(i), x_ref(offset+i))
    enddo
    call monolis_test_check_eq_R1("monolis_optimize_parallel_nnls_test 2b", residual, residual_ref)

    call monolis_dealloc_R_2d(A_loc)
    call monolis_dealloc_R_1d(x)
  end subroutine monolis_optimize_parallel_nnls_test

end module mod_monolis_opt_nnls_test