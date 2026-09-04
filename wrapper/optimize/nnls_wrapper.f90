!> NNLS モジュール
module mod_monolis_nnls_wrapper
  use mod_monolis_utils
  use mod_monolis_opt_nnls
  use iso_c_binding

  implicit none

contains

  subroutine monolis_optimize_nnls_R_c(A, b, x, m, n_loc, max_iter, tol, residual, comm) &
    & bind(c, name = "monolis_optimize_nnls_R_c_main")
    implicit none
    !> 入力行列（m x n_loc、列ブロック分散）
    real(c_double), target :: A(m*n_loc)
    !> 右辺ベクトル（m、全ランク複製）
    real(c_double), target :: b(m)
    !> 解ベクトル（n_loc）
    real(c_double), target :: x(n_loc)
    !> 行列の大きさ（行数 m）
    integer(kint_c), value :: m
    !> 行列の大きさ（自ランク保持列数 n_loc）
    integer(kint_c), value :: n_loc
    !> 最大反復回数
    integer(kint_c), value :: max_iter
    !> 収束判定閾値
    real(c_double), value :: tol
    !> 残差
    real(c_double), target :: residual
    !> MPI コミュニケータ
    integer(kint_c), value :: comm
    integer(kint) :: i, j
    real(kdouble), allocatable :: A_temp(:,:)

    !> for matrix
    call monolis_alloc_R_2d(A_temp, m, n_loc)

    do i = 1, n_loc
      do j = 1, m
        A_temp(j,i) = A(m*(i-1) + j)
      enddo
    enddo

    call monolis_optimize_nnls_R(A_temp, b, x, m, n_loc, max_iter, tol, residual, comm)
  end subroutine monolis_optimize_nnls_R_c

  subroutine monolis_optimize_nnls_R_with_sparse_solution_c(A, b, x, m, n_loc, max_iter, tol_outer, tol_inner, residual, comm) &
    & bind(c, name = "monolis_optimize_nnls_R_with_sparse_solution_c_main")
    implicit none
    !> 入力行列（m x n_loc、列ブロック分散）
    real(c_double), target :: A(m*n_loc)
    !> 右辺ベクトル（m、全ランク複製）
    real(c_double), target :: b(m)
    !> 解ベクトル（n_loc）
    real(c_double), target :: x(n_loc)
    !> 行列の大きさ（行数 m）
    integer(kint_c), value :: m
    !> 行列の大きさ（自ランク保持列数 n_loc）
    integer(kint_c), value :: n_loc
    !> 最大反復回数
    integer(kint_c), value :: max_iter
    !> 外側反復の収束判定閾値
    real(c_double), value :: tol_outer
    !> 内部 NNLS の閾値
    real(c_double), value :: tol_inner
    !> 残差
    real(c_double), target :: residual
    !> MPI コミュニケータ
    integer(kint_c), value :: comm
    integer(kint) :: i, j
    real(kdouble), allocatable :: A_temp(:,:)

    !> for matrix
    call monolis_alloc_R_2d(A_temp, m, n_loc)

    do i = 1, n_loc
      do j = 1, m
        A_temp(j,i) = A(m*(i-1) + j)
      enddo
    enddo

    call monolis_optimize_nnls_R_with_sparse_solution(A_temp, b, x, m, n_loc, max_iter, tol_outer, tol_inner, residual, comm)
  end subroutine monolis_optimize_nnls_R_with_sparse_solution_c
end module mod_monolis_nnls_wrapper
