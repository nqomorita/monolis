module monolis_mf_analysis_module
  use iso_fortran_env, only: int64
  use monolis_mf_types_module, only: d_lu_only_symbolic_analysis, &
      release_symbolic
  use monolis_mf_reorder_module, only: build_front_postorder, &
      build_pord_front_structure, front_local_position, &
      sort_int_range
  implicit none

  private
  public :: build_symbolic_analysis

contains

  subroutine build_original_front_entries(symbolic, n, row_ptr, col_ind, ierr)
    type(d_lu_only_symbolic_analysis), intent(inout), target :: symbolic
    integer, intent(in) :: n
    integer, intent(in) :: row_ptr(:)
    integer, intent(in) :: col_ind(:)
    integer, intent(out) :: ierr

    integer :: nfronts
    integer :: row
    integer :: entry
    integer :: col
    integer :: new_row
    integer :: new_col
    integer :: first_col
    integer :: front
    integer :: pos
    integer :: first_pos
    integer :: last_pos
    integer :: row_local
    integer :: col_local
    integer :: alloc_stat
    integer, allocatable :: counts(:)
    integer, allocatable :: next_pos(:)
    integer, allocatable :: local_pos(:)

    ierr = 0
    nfronts = symbolic%nfronts

    allocate(counts(nfronts), symbolic%orig_ptr(nfronts + 1), stat=alloc_stat)
    if (alloc_stat /= 0) then
      ierr = -3211
      return
    end if
    counts = 0

    do row = 1, n
      new_row = symbolic%perm(row)
      do entry = row_ptr(row), row_ptr(row + 1) - 1
        new_col = symbolic%perm(col_ind(entry))
        first_col = min(new_row, new_col)
        front = symbolic%column_to_super(first_col)
        if (front <= 0) then
          ierr = -3212
          deallocate(counts, symbolic%orig_ptr)
          return
        end if
        counts(front) = counts(front) + 1
      end do
    end do

    symbolic%orig_ptr(1) = 1_int64
    do front = 1, nfronts
      symbolic%orig_ptr(front + 1) = symbolic%orig_ptr(front) + int(counts(front), int64)
    end do

    allocate(symbolic%orig_row_pos(max(1, row_ptr(n + 1) - 1)), &
        symbolic%orig_col_pos(max(1, row_ptr(n + 1) - 1)), &
        symbolic%orig_entry(max(1, row_ptr(n + 1) - 1)), next_pos(nfronts), &
        stat=alloc_stat)
    if (alloc_stat /= 0) then
      ierr = -3213
      deallocate(counts, symbolic%orig_ptr)
      return
    end if

    next_pos = int(symbolic%orig_ptr(1:nfronts))
    do row = 1, n
      new_row = symbolic%perm(row)
      do entry = row_ptr(row), row_ptr(row + 1) - 1
        col = col_ind(entry)
        new_col = symbolic%perm(col)
        first_col = min(new_row, new_col)
        front = symbolic%column_to_super(first_col)
        pos = next_pos(front)

        symbolic%orig_row_pos(pos) = new_row
        symbolic%orig_col_pos(pos) = new_col
        symbolic%orig_entry(pos) = entry
        next_pos(front) = pos + 1
      end do
    end do

    allocate(local_pos(n), stat=alloc_stat)
    if (alloc_stat /= 0) then
      ierr = -3214
      deallocate(counts, next_pos, symbolic%orig_ptr, symbolic%orig_row_pos, &
          symbolic%orig_col_pos, symbolic%orig_entry)
      return
    end if
    local_pos = 0

    do front = 1, nfronts
      first_pos = int(symbolic%front_ptr(front))
      last_pos = int(symbolic%front_ptr(front + 1)) - 1
      do pos = first_pos, last_pos
        local_pos(symbolic%front_ind(pos)) = pos - first_pos + 1
      end do

      do pos = int(symbolic%orig_ptr(front)), int(symbolic%orig_ptr(front + 1) - 1_int64)
        row_local = local_pos(symbolic%orig_row_pos(pos))
        col_local = local_pos(symbolic%orig_col_pos(pos))
        if (row_local <= 0 .or. col_local <= 0) then
          ierr = -3215
          deallocate(counts, next_pos, local_pos, symbolic%orig_ptr, &
              symbolic%orig_row_pos, symbolic%orig_col_pos, symbolic%orig_entry)
          return
        end if
        symbolic%orig_row_pos(pos) = row_local
        symbolic%orig_col_pos(pos) = col_local
      end do

      do pos = first_pos, last_pos
        local_pos(symbolic%front_ind(pos)) = 0
      end do
    end do

    deallocate(counts, next_pos, local_pos)
  end subroutine build_original_front_entries

  subroutine build_child_contribution_maps(symbolic, ierr)
    type(d_lu_only_symbolic_analysis), intent(inout) :: symbolic
    integer, intent(out) :: ierr

    integer :: nfronts
    integer :: front
    integer :: parent
    integer :: npiv
    integer :: nupd
    integer :: i
    integer :: run
    integer :: pos
    integer :: child_pos
    integer :: run_count
    integer :: alloc_stat
    integer :: max_update_size
    integer(int64) :: total_update
    integer(int64) :: total_runs
    integer, allocatable :: parent_pos(:)

    ierr = 0
    nfronts = symbolic%nfronts
    total_update = 0_int64
    total_runs = 0_int64
    max_update_size = 0

    do front = 1, nfronts
      nupd = symbolic%front_update_size(front)
      max_update_size = max(max_update_size, nupd)
      total_update = total_update + int(nupd, int64)
    end do
    if (total_update > int(huge(1), int64)) then
      ierr = -3251
      return
    end if

    allocate(symbolic%contrib_pos_ptr(nfronts + 1), &
        symbolic%contrib_run_ptr(nfronts + 1), parent_pos(max(1, max_update_size)), &
        stat=alloc_stat)
    if (alloc_stat /= 0) then
      ierr = -3252
      return
    end if

    symbolic%contrib_pos_ptr(1) = 1_int64
    symbolic%contrib_run_ptr(1) = 1_int64
    do front = 1, nfronts
      parent = symbolic%front_parent(front)
      npiv = symbolic%front_pivot_size(front)
      nupd = symbolic%front_update_size(front)
      symbolic%contrib_pos_ptr(front + 1) = symbolic%contrib_pos_ptr(front) + &
          int(nupd, int64)

      run_count = 0
      if (parent > 0 .and. nupd > 0) then
        do i = 1, nupd
          child_pos = int(symbolic%front_ptr(front)) + npiv + i - 1
          parent_pos(i) = front_local_position(symbolic, parent, &
              symbolic%front_ind(child_pos))
          if (parent_pos(i) <= 0) then
            ierr = -3253
            deallocate(parent_pos)
            return
          end if
        end do

        i = 1
        do while (i <= nupd)
          run_count = run_count + 1
          do while (i < nupd)
            if (parent_pos(i + 1) /= parent_pos(i) + 1) exit
            i = i + 1
          end do
          i = i + 1
        end do
      end if

      total_runs = total_runs + int(run_count, int64)
      symbolic%contrib_run_ptr(front + 1) = symbolic%contrib_run_ptr(front) + &
          int(run_count, int64)
    end do

    if (total_runs > int(huge(1), int64)) then
      ierr = -3254
      deallocate(parent_pos)
      return
    end if

    allocate(symbolic%contrib_parent_pos(max(1, int(total_update))), &
        symbolic%contrib_run_first(max(1, int(total_runs))), &
        symbolic%contrib_run_len(max(1, int(total_runs))), &
        symbolic%contrib_run_parent_first(max(1, int(total_runs))), &
        stat=alloc_stat)
    if (alloc_stat /= 0) then
      ierr = -3255
      deallocate(parent_pos)
      return
    end if

    do front = 1, nfronts
      parent = symbolic%front_parent(front)
      npiv = symbolic%front_pivot_size(front)
      nupd = symbolic%front_update_size(front)
      pos = int(symbolic%contrib_pos_ptr(front))
      run = int(symbolic%contrib_run_ptr(front))

      if (parent > 0 .and. nupd > 0) then
        do i = 1, nupd
          child_pos = int(symbolic%front_ptr(front)) + npiv + i - 1
          parent_pos(i) = front_local_position(symbolic, parent, &
              symbolic%front_ind(child_pos))
          symbolic%contrib_parent_pos(pos + i - 1) = parent_pos(i)
        end do

        i = 1
        do while (i <= nupd)
          symbolic%contrib_run_first(run) = i
          symbolic%contrib_run_parent_first(run) = parent_pos(i)
          do while (i < nupd)
            if (parent_pos(i + 1) /= parent_pos(i) + 1) exit
            i = i + 1
          end do
          symbolic%contrib_run_len(run) = i - symbolic%contrib_run_first(run) + 1
          run = run + 1
          i = i + 1
        end do
      else if (nupd > 0) then
        symbolic%contrib_parent_pos(pos:pos + nupd - 1) = 0
      end if
    end do

    deallocate(parent_pos)
  end subroutine build_child_contribution_maps
  subroutine build_symbolic_analysis(symbolic, n, row_ptr, col_ind, ordering_method, ierr)
    type(d_lu_only_symbolic_analysis), intent(inout) :: symbolic
    integer, intent(in) :: n
    integer, intent(in) :: row_ptr(:)
    integer, intent(in) :: col_ind(:)
    integer, intent(in) :: ordering_method
    integer, intent(out) :: ierr

    ierr = 0
    call release_symbolic(symbolic)

    symbolic%n = n
    call build_symmetric_pattern(symbolic, n, row_ptr, col_ind, ierr)
    if (ierr /= 0) return

    call build_pord_front_structure(symbolic, ierr)
    if (ierr /= 0) return

    call build_original_front_entries(symbolic, n, row_ptr, col_ind, ierr)
    if (ierr /= 0) return

    call build_child_contribution_maps(symbolic, ierr)
  end subroutine build_symbolic_analysis

  subroutine build_symmetric_pattern(symbolic, n, row_ptr, col_ind, ierr)
    type(d_lu_only_symbolic_analysis), intent(inout) :: symbolic
    integer, intent(in) :: n
    integer, intent(in) :: row_ptr(:)
    integer, intent(in) :: col_ind(:)
    integer, intent(out) :: ierr

    integer :: row
    integer :: entry
    integer :: col
    integer :: pos
    integer :: write_pos
    integer :: count
    integer :: prev
    integer :: alloc_stat
    integer(int64) :: raw_nnz_64
    integer, allocatable :: counts(:)
    integer, allocatable :: raw_ptr(:)
    integer, allocatable :: next_pos(:)
    integer, allocatable :: raw_cols(:)
    integer, allocatable :: compact_cols(:)

    ierr = 0
    allocate(counts(n), stat=alloc_stat)
    if (alloc_stat /= 0) then
      ierr = -2401
      return
    end if
    counts = 0

    do row = 1, n
      do entry = row_ptr(row), row_ptr(row + 1) - 1
        col = col_ind(entry)
        counts(row) = counts(row) + 1
        if (col /= row) counts(col) = counts(col) + 1
      end do
    end do

    allocate(raw_ptr(n + 1), next_pos(n), stat=alloc_stat)
    if (alloc_stat /= 0) then
      ierr = -2402
      deallocate(counts)
      return
    end if

    raw_ptr(1) = 1
    do row = 1, n
      raw_ptr(row + 1) = raw_ptr(row) + counts(row)
    end do
    raw_nnz_64 = int(raw_ptr(n + 1), int64) - 1_int64
    if (raw_nnz_64 > int(huge(1), int64)) then
      ierr = -2403
      deallocate(counts, raw_ptr, next_pos)
      return
    end if

    allocate(raw_cols(int(raw_nnz_64)), stat=alloc_stat)
    if (alloc_stat /= 0) then
      ierr = -2404
      deallocate(counts, raw_ptr, next_pos)
      return
    end if

    next_pos = raw_ptr(1:n)
    do row = 1, n
      do entry = row_ptr(row), row_ptr(row + 1) - 1
        col = col_ind(entry)
        pos = next_pos(row)
        raw_cols(pos) = col
        next_pos(row) = pos + 1

        if (col /= row) then
          pos = next_pos(col)
          raw_cols(pos) = row
          next_pos(col) = pos + 1
        end if
      end do
    end do

    allocate(symbolic%sym_row_ptr(n + 1), symbolic%sym_col_ind(size(raw_cols)), stat=alloc_stat)
    if (alloc_stat /= 0) then
      ierr = -2405
      deallocate(counts, raw_ptr, next_pos, raw_cols)
      return
    end if

    write_pos = 1
    symbolic%sym_row_ptr(1) = 1
    do row = 1, n
      call sort_int_range(raw_cols, raw_ptr(row), raw_ptr(row + 1) - 1)
      prev = 0
      count = 0
      do entry = raw_ptr(row), raw_ptr(row + 1) - 1
        if (count == 0 .or. raw_cols(entry) /= prev) then
          symbolic%sym_col_ind(write_pos) = raw_cols(entry)
          write_pos = write_pos + 1
          count = count + 1
          prev = raw_cols(entry)
        end if
      end do
      symbolic%sym_row_ptr(row + 1) = write_pos
    end do

    symbolic%sym_nnz = write_pos - 1
    allocate(compact_cols(symbolic%sym_nnz), stat=alloc_stat)
    if (alloc_stat /= 0) then
      ierr = -2406
      deallocate(counts, raw_ptr, next_pos, raw_cols)
      return
    end if
    compact_cols = symbolic%sym_col_ind(1:symbolic%sym_nnz)
    call move_alloc(compact_cols, symbolic%sym_col_ind)

    deallocate(counts, raw_ptr, next_pos, raw_cols)
  end subroutine build_symmetric_pattern


end module monolis_mf_analysis_module
