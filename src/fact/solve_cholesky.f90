!> 多重フロント法 Cholesky の前進・後退代入
module mod_monolis_fact_solve_cholesky
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_lapack

  implicit none

  private
  public :: monolis_fact_solve_cholesky

contains

  !> @ingroup fact
  !> 前進・後退代入により L L^T x = rhs を解く
  subroutine monolis_fact_solve_cholesky(lu, rhs)
    implicit none
    !> [in] 分解構造体（Cholesky 因子を保持）
    type(monolis_mat_lu), intent(in)    :: lu
    !> [in,out] 入力時に rhs、出力時に解
    real(kdouble),        intent(inout) :: rhs(:)

    integer(kint) :: n, nfronts, order_pos, front
    integer(kint) :: fs, npiv, nupd, first_col, i, idx, ldp
    real(kdouble), allocatable :: work(:), update_work(:), pivot_work(:,:)

    call monolis_std_debug_log_header("monolis_fact_solve_cholesky")

    if (.not. lu%factorized) then
      call monolis_std_error_string("monolis_fact_solve_cholesky: factorization not done")
      call monolis_std_error_stop()
    end if

    n = lu%N
    nfronts = lu%nfronts
    if (n <= 0) return

    ldp = max(1, lu%max_front_size)
    call monolis_alloc_R_1d(work, n)
    call monolis_alloc_R_1d(update_work, ldp)
    call monolis_alloc_R_2d(pivot_work, ldp, 1)

    do i = 1, n
      work(i) = rhs(lu%iperm(i))
    end do

    !> 前進代入 L y = b（後順走査）
    do order_pos = 1, nfronts
      front = lu%front_postorder(order_pos)
      fs   = lu%factors(front)%front_size
      npiv = lu%factors(front)%pivot_size
      nupd = lu%factors(front)%update_size
      first_col = lu%super_start(front)

      if (npiv > 0) then
        pivot_work(1:npiv, 1) = work(first_col:first_col + npiv - 1)
        call dtrsm('L', 'L', 'N', 'N', npiv, 1, 1.0d0, &
            lu%factors(front)%factor, max(1, fs), pivot_work, ldp)
        work(first_col:first_col + npiv - 1) = pivot_work(1:npiv, 1)
      end if

      if (nupd > 0) then
        do i = 1, nupd
          idx = lu%front_ind(lu%front_ptr(front) + npiv + i - 1)
          update_work(i) = work(idx)
        end do
        call dgemv('N', nupd, npiv, -1.0d0, &
            lu%factors(front)%factor(npiv + 1, 1), max(1, fs), &
            pivot_work(1:npiv, 1), 1, 1.0d0, update_work, 1)
        do i = 1, nupd
          idx = lu%front_ind(lu%front_ptr(front) + npiv + i - 1)
          work(idx) = update_work(i)
        end do
      end if
    end do

    !> 後退代入 L^T x = y（逆順、上三角因子は L の転置で代用）
    do order_pos = nfronts, 1, -1
      front = lu%front_postorder(order_pos)
      fs   = lu%factors(front)%front_size
      npiv = lu%factors(front)%pivot_size
      nupd = lu%factors(front)%update_size
      first_col = lu%super_start(front)

      if (nupd > 0) then
        pivot_work(1:npiv, 1) = work(first_col:first_col + npiv - 1)
        do i = 1, nupd
          idx = lu%front_ind(lu%front_ptr(front) + npiv + i - 1)
          update_work(i) = work(idx)
        end do
        call dgemv('T', nupd, npiv, -1.0d0, &
            lu%factors(front)%factor(npiv + 1, 1), max(1, fs), &
            update_work, 1, 1.0d0, pivot_work(1:npiv, 1), 1)
        work(first_col:first_col + npiv - 1) = pivot_work(1:npiv, 1)
      end if

      if (npiv > 0) then
        pivot_work(1:npiv, 1) = work(first_col:first_col + npiv - 1)
        call dtrsm('L', 'L', 'T', 'N', npiv, 1, 1.0d0, &
            lu%factors(front)%factor, max(1, fs), pivot_work, ldp)
        work(first_col:first_col + npiv - 1) = pivot_work(1:npiv, 1)
      end if
    end do

    do i = 1, n
      rhs(lu%iperm(i)) = work(i)
    end do

    call monolis_dealloc_R_1d(work)
    call monolis_dealloc_R_1d(update_work)
    call monolis_dealloc_R_2d(pivot_work)
  end subroutine monolis_fact_solve_cholesky

end module mod_monolis_fact_solve_cholesky
