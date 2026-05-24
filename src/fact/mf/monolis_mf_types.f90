module monolis_mf_types_module
  use iso_fortran_env, only: int64, real64
  implicit none

  private
  public :: front_block_size
  public :: dlu_ordering_pord
  public :: d_lu_only_symbolic_analysis
  public :: d_lu_only_numeric_front
  public :: monolis_mf_handle
  public :: release_all
  public :: release_numeric
  public :: release_symbolic

  integer, parameter :: front_block_size = 128
  integer, parameter :: dlu_ordering_pord = 1

  type :: d_lu_only_symbolic_analysis
    integer :: n = 0
    integer :: sym_nnz = 0
    integer :: nsuper = 0
    integer :: nfronts = 0
    integer :: max_front_size = 0
    integer(int64) :: l_pattern_nnz = 0_int64
    integer, allocatable :: sym_row_ptr(:)
    integer, allocatable :: sym_col_ind(:)
    integer, allocatable :: perm(:)
    integer, allocatable :: invp(:)
    integer, allocatable :: vertex_front(:)
    integer, allocatable :: super_start(:)
    integer, allocatable :: super_end(:)
    integer, allocatable :: column_to_super(:)
    integer, allocatable :: super_parent(:)
    integer(int64), allocatable :: front_ptr(:)
    integer, allocatable :: front_ind(:)
    integer, allocatable :: front_parent(:)
    integer, allocatable :: front_first_child(:)
    integer, allocatable :: front_next_sibling(:)
    integer, allocatable :: front_postorder(:)
    integer, allocatable :: front_size(:)
    integer, allocatable :: front_pivot_size(:)
    integer, allocatable :: front_update_size(:)
    integer(int64), allocatable :: orig_ptr(:)
    integer, allocatable :: orig_row_pos(:)
    integer, allocatable :: orig_col_pos(:)
    integer, allocatable :: orig_entry(:)
    integer(int64), allocatable :: contrib_pos_ptr(:)
    integer, allocatable :: contrib_parent_pos(:)
    integer(int64), allocatable :: contrib_run_ptr(:)
    integer, allocatable :: contrib_run_first(:)
    integer, allocatable :: contrib_run_len(:)
    integer, allocatable :: contrib_run_parent_first(:)
  end type d_lu_only_symbolic_analysis

  type :: d_lu_only_numeric_front
    integer :: front_size = 0
    integer :: pivot_size = 0
    integer :: update_size = 0
    real(real64), allocatable :: factor(:,:)
    real(real64), allocatable :: upper_update(:,:)
    real(real64), allocatable :: contribution(:,:)
  end type d_lu_only_numeric_front

  type :: monolis_mf_handle
    integer :: n = 0
    integer(int64) :: nnz = 0_int64
    integer :: kl = 0
    integer :: ku = 0
    integer :: ldab = 0
    integer :: ldlu = 0
    integer, allocatable :: row_ptr(:)
    integer, allocatable :: col_ind(:)
    real(real64), allocatable :: lu_band(:,:)
    type(d_lu_only_numeric_front), allocatable :: fronts(:)
    type(d_lu_only_symbolic_analysis) :: symbolic
    integer :: ordering_method = dlu_ordering_pord
    logical :: analyzed = .false.
    logical :: factorized = .false.
  end type monolis_mf_handle

contains

  subroutine release_all(handle)
    type(monolis_mf_handle), intent(inout) :: handle

    if (allocated(handle%row_ptr)) deallocate(handle%row_ptr)
    if (allocated(handle%col_ind)) deallocate(handle%col_ind)
    call release_numeric(handle)
    call release_symbolic(handle%symbolic)
    handle%n = 0
    handle%nnz = 0_int64
    handle%kl = 0
    handle%ku = 0
    handle%ldab = 0
    handle%ldlu = 0
    handle%analyzed = .false.
    handle%factorized = .false.
  end subroutine release_all

  subroutine release_numeric(handle)
    type(monolis_mf_handle), intent(inout) :: handle

    integer :: front

    if (allocated(handle%lu_band)) deallocate(handle%lu_band)
    if (allocated(handle%fronts)) then
      do front = 1, size(handle%fronts)
        if (allocated(handle%fronts(front)%factor)) deallocate(handle%fronts(front)%factor)
        if (allocated(handle%fronts(front)%upper_update)) then
          deallocate(handle%fronts(front)%upper_update)
        end if
        if (allocated(handle%fronts(front)%contribution)) then
          deallocate(handle%fronts(front)%contribution)
        end if
        handle%fronts(front)%front_size = 0
        handle%fronts(front)%pivot_size = 0
        handle%fronts(front)%update_size = 0
      end do
      deallocate(handle%fronts)
    end if
    handle%factorized = .false.
  end subroutine release_numeric

  subroutine release_symbolic(symbolic)
    type(d_lu_only_symbolic_analysis), intent(inout) :: symbolic

    if (allocated(symbolic%sym_row_ptr)) deallocate(symbolic%sym_row_ptr)
    if (allocated(symbolic%sym_col_ind)) deallocate(symbolic%sym_col_ind)
    if (allocated(symbolic%super_start)) deallocate(symbolic%super_start)
    if (allocated(symbolic%super_end)) deallocate(symbolic%super_end)
    if (allocated(symbolic%column_to_super)) deallocate(symbolic%column_to_super)
    if (allocated(symbolic%super_parent)) deallocate(symbolic%super_parent)
    if (allocated(symbolic%front_ptr)) deallocate(symbolic%front_ptr)
    if (allocated(symbolic%front_ind)) deallocate(symbolic%front_ind)
    if (allocated(symbolic%front_parent)) deallocate(symbolic%front_parent)
    if (allocated(symbolic%front_first_child)) deallocate(symbolic%front_first_child)
    if (allocated(symbolic%front_next_sibling)) deallocate(symbolic%front_next_sibling)
    if (allocated(symbolic%front_postorder)) deallocate(symbolic%front_postorder)
    if (allocated(symbolic%front_size)) deallocate(symbolic%front_size)
    if (allocated(symbolic%front_pivot_size)) deallocate(symbolic%front_pivot_size)
    if (allocated(symbolic%front_update_size)) deallocate(symbolic%front_update_size)
    if (allocated(symbolic%orig_ptr)) deallocate(symbolic%orig_ptr)
    if (allocated(symbolic%orig_row_pos)) deallocate(symbolic%orig_row_pos)
    if (allocated(symbolic%orig_col_pos)) deallocate(symbolic%orig_col_pos)
    if (allocated(symbolic%orig_entry)) deallocate(symbolic%orig_entry)
    if (allocated(symbolic%contrib_pos_ptr)) deallocate(symbolic%contrib_pos_ptr)
    if (allocated(symbolic%contrib_parent_pos)) deallocate(symbolic%contrib_parent_pos)
    if (allocated(symbolic%contrib_run_ptr)) deallocate(symbolic%contrib_run_ptr)
    if (allocated(symbolic%contrib_run_first)) deallocate(symbolic%contrib_run_first)
    if (allocated(symbolic%contrib_run_len)) deallocate(symbolic%contrib_run_len)
    if (allocated(symbolic%contrib_run_parent_first)) then
      deallocate(symbolic%contrib_run_parent_first)
    end if
    symbolic%n = 0
    symbolic%sym_nnz = 0
    symbolic%nsuper = 0
    symbolic%nfronts = 0
    symbolic%max_front_size = 0
    symbolic%l_pattern_nnz = 0_int64
  end subroutine release_symbolic

end module monolis_mf_types_module
