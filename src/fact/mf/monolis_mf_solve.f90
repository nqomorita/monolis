module monolis_mf_solve_module
  use iso_fortran_env, only: real64
  use monolis_mf_types_module, only: monolis_mf_handle, &
      d_lu_only_symbolic_analysis
  use monolis_mf_lapack_wrapper_module, only: dgemv, dtrsm
  implicit none

  private
  public :: solve_dense_front_numeric

contains

  subroutine solve_dense_front_numeric(handle, rhs, ierr)
    type(monolis_mf_handle), intent(in), target :: handle
    real(real64), intent(inout) :: rhs(:)
    integer, intent(out) :: ierr

    type(d_lu_only_symbolic_analysis), pointer :: symbolic
    integer :: n
    integer :: nfronts
    integer :: order_pos
    integer :: front
    integer :: fs
    integer :: npiv
    integer :: nupd
    integer :: first_col
    integer :: i
    integer :: idx
    integer :: alloc_stat
    real(real64), allocatable :: work(:)
    real(real64), allocatable :: update_work(:)
    real(real64), allocatable :: pivot_work(:,:)

    ierr = 0
    symbolic => handle%symbolic
    n = handle%n
    nfronts = symbolic%nfronts

    allocate(work(n), update_work(max(1, symbolic%max_front_size)), &
        pivot_work(max(1, symbolic%max_front_size), 1), stat=alloc_stat)
    if (alloc_stat /= 0) then
      ierr = -3231
      return
    end if

    do i = 1, n
      work(i) = rhs(symbolic%invp(i))
    end do

    do order_pos = 1, nfronts
      front = symbolic%front_postorder(order_pos)
      fs = handle%fronts(front)%front_size
      npiv = handle%fronts(front)%pivot_size
      nupd = handle%fronts(front)%update_size
      first_col = symbolic%super_start(front)

      if (npiv > 0) then
        pivot_work(1:npiv, 1) = work(first_col:first_col + npiv - 1)
        call dtrsm('L', 'L', 'N', 'U', npiv, 1, 1.0_real64, &
            handle%fronts(front)%factor, max(1, fs), pivot_work, &
            max(1, symbolic%max_front_size))
        work(first_col:first_col + npiv - 1) = pivot_work(1:npiv, 1)
      end if
      if (nupd > 0) then
        do i = 1, nupd
          idx = symbolic%front_ind(int(symbolic%front_ptr(front)) + npiv + i - 1)
          update_work(i) = work(idx)
        end do
        call dgemv('N', nupd, npiv, -1.0_real64, &
            handle%fronts(front)%factor(npiv + 1, 1), max(1, fs), &
            pivot_work(1:npiv, 1), 1, &
            1.0_real64, update_work, 1)
        do i = 1, nupd
          idx = symbolic%front_ind(int(symbolic%front_ptr(front)) + npiv + i - 1)
          work(idx) = update_work(i)
        end do
      end if
    end do

    do order_pos = nfronts, 1, -1
      front = symbolic%front_postorder(order_pos)
      fs = handle%fronts(front)%front_size
      npiv = handle%fronts(front)%pivot_size
      nupd = handle%fronts(front)%update_size
      first_col = symbolic%super_start(front)

      if (nupd > 0) then
        pivot_work(1:npiv, 1) = work(first_col:first_col + npiv - 1)
        do i = 1, nupd
          idx = symbolic%front_ind(int(symbolic%front_ptr(front)) + npiv + i - 1)
          update_work(i) = work(idx)
        end do
        call dgemv('N', npiv, nupd, -1.0_real64, &
            handle%fronts(front)%upper_update, max(1, npiv), update_work, 1, &
            1.0_real64, pivot_work(1:npiv, 1), 1)
        work(first_col:first_col + npiv - 1) = pivot_work(1:npiv, 1)
      end if
      if (npiv > 0) then
        pivot_work(1:npiv, 1) = work(first_col:first_col + npiv - 1)
        call dtrsm('L', 'U', 'N', 'N', npiv, 1, 1.0_real64, &
            handle%fronts(front)%factor, max(1, fs), pivot_work, &
            max(1, symbolic%max_front_size))
        work(first_col:first_col + npiv - 1) = pivot_work(1:npiv, 1)
      end if
    end do

    do i = 1, n
      rhs(symbolic%invp(i)) = work(i)
    end do

    deallocate(work, update_work, pivot_work)
  end subroutine solve_dense_front_numeric

end module monolis_mf_solve_module
