module monolis_mf_driver_module
  use iso_fortran_env, only: int64, real64
  use monolis_mf_types_module, only: monolis_mf_handle, &
      release_all, release_numeric
  use monolis_mf_analysis_module, only: build_symbolic_analysis
  use monolis_mf_factor_module, only: factor_dense_front_numeric
  use monolis_mf_solve_module, only: solve_dense_front_numeric
  implicit none

  private
  public :: monolis_mf_analyze_csr
  public :: monolis_mf_factor
  public :: monolis_mf_finalize
  public :: monolis_mf_solve
  public :: monolis_mf_symbolic_stats

contains

  subroutine monolis_mf_analyze_csr(handle, n, row_ptr, col_ind, ierr)
    type(monolis_mf_handle), intent(inout) :: handle
    integer, intent(in) :: n
    integer, intent(in) :: row_ptr(:)
    integer, intent(in) :: col_ind(:)
    integer, intent(out) :: ierr

    integer :: row
    integer :: entry
    integer :: col
    integer :: alloc_stat
    integer :: kl
    integer :: ku
    integer(int64) :: nnz

    ierr = 0
    call release_all(handle)

    if (n <= 0) then
      ierr = -2001
      return
    end if
    if (size(row_ptr) /= n + 1) then
      ierr = -2002
      return
    end if
    if (row_ptr(1) /= 1) then
      ierr = -2003
      return
    end if

    do row = 1, n
      if (row_ptr(row + 1) < row_ptr(row)) then
        ierr = -2004
        return
      end if
    end do

    nnz = int(row_ptr(n + 1), int64) - 1_int64
    if (nnz <= 0_int64) then
      ierr = -2005
      return
    end if
    if (nnz > int(huge(1), int64)) then
      ierr = -2006
      return
    end if
    if (size(col_ind, kind=int64) /= nnz) then
      ierr = -2007
      return
    end if

    kl = 0
    ku = 0
    do row = 1, n
      do entry = row_ptr(row), row_ptr(row + 1) - 1
        col = col_ind(entry)
        if (col < 1 .or. col > n) then
          ierr = -2008
          return
        end if
        kl = max(kl, row - col)
        ku = max(ku, col - row)
      end do
    end do

    allocate(handle%row_ptr(size(row_ptr)), handle%col_ind(size(col_ind)), &
        stat=alloc_stat)
    if (alloc_stat /= 0) then
      ierr = -2009
      call release_all(handle)
      return
    end if

    handle%n = n
    handle%nnz = nnz
    handle%kl = kl
    handle%ku = ku
    handle%ldab = 2 * kl + ku + 1
    handle%ldlu = kl + ku + 1
    handle%row_ptr = row_ptr
    handle%col_ind = col_ind

    call build_symbolic_analysis(handle%symbolic, n, row_ptr, col_ind, &
        handle%ordering_method, ierr)
    if (ierr /= 0) then
      call release_all(handle)
      return
    end if

    handle%analyzed = .true.
    handle%factorized = .false.
  end subroutine monolis_mf_analyze_csr

  subroutine monolis_mf_factor(handle, values, ierr)
    type(monolis_mf_handle), intent(inout) :: handle
    real(real64), intent(in) :: values(:)
    integer, intent(out) :: ierr

    ierr = 0
    if (.not. handle%analyzed) then
      ierr = -2010
      return
    end if
    if (size(values, kind=int64) /= handle%nnz) then
      ierr = -2011
      return
    end if

    call release_numeric(handle)
    call factor_dense_front_numeric(handle, values, ierr)
    if (ierr /= 0) return

    handle%factorized = .true.
  end subroutine monolis_mf_factor

  subroutine monolis_mf_solve(handle, rhs, ierr)
    type(monolis_mf_handle), intent(inout) :: handle
    real(real64), intent(inout) :: rhs(:)
    integer, intent(out) :: ierr

    ierr = 0
    if (.not. handle%factorized) then
      ierr = -2012
      return
    end if
    if (size(rhs) /= handle%n) then
      ierr = -2013
      return
    end if

    call solve_dense_front_numeric(handle, rhs, ierr)
  end subroutine monolis_mf_solve

  subroutine monolis_mf_symbolic_stats(handle, sym_nnz, l_pattern_nnz, &
      nsuper, nfronts, max_front_size, ierr)
    type(monolis_mf_handle), intent(in) :: handle
    integer, intent(out) :: sym_nnz
    integer(int64), intent(out) :: l_pattern_nnz
    integer, intent(out) :: nsuper
    integer, intent(out) :: nfronts
    integer, intent(out) :: max_front_size
    integer, intent(out) :: ierr

    ierr = 0
    if (.not. handle%analyzed) then
      ierr = -2301
      sym_nnz = 0
      l_pattern_nnz = 0_int64
      nsuper = 0
      nfronts = 0
      max_front_size = 0
      return
    end if

    sym_nnz = handle%symbolic%sym_nnz
    l_pattern_nnz = handle%symbolic%l_pattern_nnz
    nsuper = handle%symbolic%nsuper
    nfronts = handle%symbolic%nfronts
    max_front_size = handle%symbolic%max_front_size
  end subroutine monolis_mf_symbolic_stats

  subroutine monolis_mf_finalize(handle, ierr)
    type(monolis_mf_handle), intent(inout) :: handle
    integer, intent(out) :: ierr

    ierr = 0
    call release_all(handle)
  end subroutine monolis_mf_finalize

end module monolis_mf_driver_module
