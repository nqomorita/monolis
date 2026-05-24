!> 多重フロント法 LU の前進・後退代入
module mod_monolis_fact_solve
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_fact_blas

  implicit none

  private
  public :: monolis_fact_solve

contains

  !> @ingroup fact
  !> 前進・後退代入により LU x = rhs を解く
  subroutine monolis_fact_solve(lu, rhs)
    implicit none
    !> [in] LU 分解構造体
    type(monolis_mat_lu), intent(in)    :: lu
    !> [in,out] 入力時に rhs、出力時に解
    real(kdouble),        intent(inout) :: rhs(:)

    integer(kint) :: n, nfronts, order_pos, front
    integer(kint) :: fs, npiv, nupd, first_col, i, idx, ldp
    real(kdouble), allocatable :: work(:), update_work(:), pivot_work(:,:)

    call monolis_std_debug_log_header("monolis_fact_solve")

    if (.not. lu%factorized) then
      call monolis_std_error_string("monolis_fact_solve: factorization not done")
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

    !> 前進代入（後順走査）
    do order_pos = 1, nfronts
      front = lu%front_postorder(order_pos)
      fs   = lu%factors(front)%front_size
      npiv = lu%factors(front)%pivot_size
      nupd = lu%factors(front)%update_size
      first_col = lu%super_start(front)

      if (npiv > 0) then
        pivot_work(1:npiv, 1) = work(first_col:first_col + npiv - 1)
        call dtrsm('L', 'L', 'N', 'U', npiv, 1, 1.0d0, &
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

    !> 後退代入（逆順）
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
        call dgemv('N', npiv, nupd, -1.0d0, &
            lu%factors(front)%upper_update, max(1, npiv), update_work, 1, &
            1.0d0, pivot_work(1:npiv, 1), 1)
        work(first_col:first_col + npiv - 1) = pivot_work(1:npiv, 1)
      end if

      if (npiv > 0) then
        pivot_work(1:npiv, 1) = work(first_col:first_col + npiv - 1)
        call dtrsm('L', 'U', 'N', 'N', npiv, 1, 1.0d0, &
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
  end subroutine monolis_fact_solve

end module mod_monolis_fact_solve
