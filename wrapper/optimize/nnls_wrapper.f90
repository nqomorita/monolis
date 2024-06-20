!> NNLS モジュール
module mod_monolis_nnls_wrapper
  use mod_monolis_utils
  use mod_monolis_opt_nnls
  use iso_c_binding

  implicit none

contains

  subroutine monolis_optimize_nnls_R_c(A, b, x, m, n, max_iter, tol, residual, comm) &
    & bind(c, name = "monolis_optimize_nnls_R_c_main")
    implicit none
    !> 入力行列（m x n）
    real(c_double), target :: A(m*n)
    !> 特異値（m）
    real(c_double), target :: b(m)
    !> 特異値（n）
    real(c_double), target :: x(n)
    !> 行列の大きさ（行数 m）
    integer(c_int), value :: m
    !> 行列の大きさ（列数 n）
    integer(c_int), value :: n
    !> 最大反復回数
    integer(c_int), value :: max_iter
    !> 特異値（n）
    real(c_double), value :: tol
    !> 特異値（n）
    real(c_double), target :: residual
    !> コミュニケータ
    integer(c_int), value :: comm
    integer(kint) :: i, j
    real(kdouble), allocatable :: A_temp(:,:)

    call monolis_alloc_R_2d(A_temp, m, n)

    do i = 1, n
      do j = 1, m
        A_temp(j,i) = A(n*(i-1) + j)
      enddo
    enddo

    call monolis_optimize_nnls_R(A_temp, b, x, m, n, max_iter, tol, residual, comm)
  end subroutine monolis_optimize_nnls_R_c

  subroutine monolis_optimize_nnls_R_with_sparse_solution_c(A, b, x, m, n, max_iter, tol, residual, comm) &
    & bind(c, name = "monolis_optimize_nnls_R_with_sparse_solution_c_main")
    implicit none
    !> 入力行列（m x n）
    real(c_double), target :: A(m*n)
    !> 特異値（m）
    real(c_double), target :: b(m)
    !> 特異値（n）
    real(c_double), target :: x(n)
    !> 行列の大きさ（行数 m）
    integer(c_int), value :: m
    !> 行列の大きさ（列数 n）
    integer(c_int), value :: n
    !> 最大反復回数
    integer(c_int), value :: max_iter
    !> 特異値（n）
    real(c_double), value :: tol
    !> 特異値（n）
    real(c_double), target :: residual
    !> コミュニケータ
    integer(c_int), value :: comm
    integer(kint) :: i, j
    real(kdouble), allocatable :: A_temp(:,:)

    call monolis_alloc_R_2d(A_temp, m, n)

    do i = 1, n
      do j = 1, m
        A_temp(j,i) = A(n*(i-1) + j)
      enddo
    enddo

    call monolis_optimize_nnls_R_with_sparse_solution(A_temp, b, x, m, n, max_iter, tol, residual, comm)
  end subroutine monolis_optimize_nnls_R_with_sparse_solution_c
end module mod_monolis_nnls_wrapper
