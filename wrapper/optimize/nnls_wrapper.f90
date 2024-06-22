!> NNLS モジュール
module mod_monolis_nnls_wrapper
  use mod_monolis_utils
  use mod_monolis_opt_nnls
  use iso_c_binding

  implicit none

contains

  subroutine monolis_optimize_nnls_R_c(A, b, x, m, n, max_iter, tol, residual, &
    my_rank, comm, comm_size, &
    recv_n_neib, recv_nitem, recv_neib_pe, recv_index, recv_item, &
    send_n_neib, send_nitem, send_neib_pe, send_index, send_item) &
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
    integer(c_int), intent(in), value :: my_rank, comm, comm_size
    integer(c_int), intent(in), value :: recv_n_neib, send_n_neib, recv_nitem, send_nitem
    integer(c_int), intent(in), target :: recv_neib_pe(recv_n_neib)
    integer(c_int), intent(in), target :: recv_index(recv_n_neib + 1), recv_item(recv_nitem)
    integer(c_int), intent(in), target :: send_neib_pe(send_n_neib)
    integer(c_int), intent(in), target :: send_index(send_n_neib + 1), send_item(send_nitem)
    type(monolis_com) :: monoCOM
    integer(kint) :: i, j
    real(kdouble), allocatable :: A_temp(:,:)

    !> for monoCOM
    monoCOM%my_rank = my_rank
    monoCOM%comm = comm
    monoCOM%comm_size = comm_size
    monoCOM%recv_n_neib = recv_n_neib
    monoCOM%recv_neib_pe => recv_neib_pe
    monoCOM%recv_index => recv_index
    monoCOM%recv_item => recv_item
    monoCOM%send_n_neib = send_n_neib
    monoCOM%send_neib_pe => send_neib_pe
    monoCOM%send_index => send_index
    monoCOM%send_item => send_item

    !> for matrix
    call monolis_alloc_R_2d(A_temp, m, n)

    do i = 1, n
      do j = 1, m
        A_temp(j,i) = A(m*(i-1) + j)
      enddo
    enddo

    call monolis_optimize_nnls_R(A_temp, b, x, m, n, max_iter, tol, residual, monoCOM)
  end subroutine monolis_optimize_nnls_R_c

  subroutine monolis_optimize_nnls_R_with_sparse_solution_c(A, b, x, m, n, max_iter, tol, residual, &
    my_rank, comm, comm_size, &
    recv_n_neib, recv_nitem, recv_neib_pe, recv_index, recv_item, &
    send_n_neib, send_nitem, send_neib_pe, send_index, send_item) &
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
    integer(c_int), intent(in), value :: my_rank, comm, comm_size
    integer(c_int), intent(in), value :: recv_n_neib, send_n_neib, recv_nitem, send_nitem
    integer(c_int), intent(in), target :: recv_neib_pe(recv_n_neib)
    integer(c_int), intent(in), target :: recv_index(recv_n_neib + 1), recv_item(recv_nitem)
    integer(c_int), intent(in), target :: send_neib_pe(send_n_neib)
    integer(c_int), intent(in), target :: send_index(send_n_neib + 1), send_item(send_nitem)
    type(monolis_com) :: monoCOM
    integer(kint) :: i, j
    real(kdouble), allocatable :: A_temp(:,:)

    !> for monoCOM
    monoCOM%my_rank = my_rank
    monoCOM%comm = comm
    monoCOM%comm_size = comm_size
    monoCOM%recv_n_neib = recv_n_neib
    monoCOM%recv_neib_pe => recv_neib_pe
    monoCOM%recv_index => recv_index
    monoCOM%recv_item => recv_item
    monoCOM%send_n_neib = send_n_neib
    monoCOM%send_neib_pe => send_neib_pe
    monoCOM%send_index => send_index
    monoCOM%send_item => send_item

    !> for matrix
    call monolis_alloc_R_2d(A_temp, m, n)

    do i = 1, n
      do j = 1, m
        A_temp(j,i) = A(m*(i-1) + j)
      enddo
    enddo

    call monolis_optimize_nnls_R_with_sparse_solution(A_temp, b, x, m, n, max_iter, tol, residual, monoCOM)
  end subroutine monolis_optimize_nnls_R_with_sparse_solution_c
end module mod_monolis_nnls_wrapper
