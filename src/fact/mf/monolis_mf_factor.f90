module monolis_mf_factor_module
  use iso_fortran_env, only: int64, real64
  use monolis_mf_types_module, only: monolis_mf_handle, &
      d_lu_only_numeric_front, d_lu_only_symbolic_analysis, front_block_size, &
      release_numeric
  use monolis_mf_lapack_wrapper_module, only: daxpy, dgemm, dlaset, &
      dscal, dtrsm
  implicit none

  private
  public :: factor_dense_front_numeric

contains

  subroutine factor_dense_front_numeric(handle, values, ierr)
    type(monolis_mf_handle), intent(inout), target :: handle
    real(real64), intent(in) :: values(:)
    integer, intent(out) :: ierr

    type(d_lu_only_symbolic_analysis), pointer :: symbolic
    integer :: nfronts
    integer :: order_pos
    integer :: front
    integer :: child
    integer :: fs
    integer :: npiv
    integer :: nupd
    integer :: alloc_stat
    integer(int64) :: entry_pos

    ierr = 0
    symbolic => handle%symbolic
    nfronts = symbolic%nfronts

    allocate(handle%fronts(nfronts), stat=alloc_stat)
    if (alloc_stat /= 0) then
      ierr = -3201
      call release_numeric(handle)
      return
    end if

    do order_pos = 1, nfronts
      front = symbolic%front_postorder(order_pos)
      fs = symbolic%front_size(front)
      npiv = symbolic%front_pivot_size(front)
      nupd = symbolic%front_update_size(front)

      call allocate_front_storage(handle%fronts(front), fs, npiv, nupd, ierr)
      if (ierr /= 0) then
        call release_numeric(handle)
        return
      end if

      child = symbolic%front_first_child(front)
      do while (child /= 0)
        call assemble_child_contribution(handle, child, front, handle%fronts(front), ierr)
        if (ierr /= 0) then
          call release_numeric(handle)
          return
        end if
        child = symbolic%front_next_sibling(child)
      end do

      do entry_pos = symbolic%orig_ptr(front), symbolic%orig_ptr(front + 1) - 1_int64
        call add_front_entry(handle%fronts(front), &
            symbolic%orig_row_pos(entry_pos), symbolic%orig_col_pos(entry_pos), &
            values(symbolic%orig_entry(entry_pos)))
      end do

      call factor_one_front(handle%fronts(front), front, ierr)
      if (ierr /= 0) then
        call release_numeric(handle)
        return
      end if
    end do

    do front = 1, nfronts
      if (allocated(handle%fronts(front)%contribution)) then
        deallocate(handle%fronts(front)%contribution)
      end if
    end do
  end subroutine factor_dense_front_numeric
  subroutine assemble_child_contribution(handle, child, parent, parent_data, ierr)
    type(monolis_mf_handle), intent(inout), target :: handle
    integer, intent(in) :: child
    integer, intent(in) :: parent
    type(d_lu_only_numeric_front), intent(inout) :: parent_data
    integer, intent(out) :: ierr

    type(d_lu_only_symbolic_analysis), pointer :: symbolic
    integer :: child_nupd
    integer :: j
    integer :: parent_j
    integer :: run
    integer :: pos_base
    integer :: run_first
    integer :: run_last

    ierr = 0
    symbolic => handle%symbolic
    child_nupd = symbolic%front_update_size(child)
    if (child_nupd <= 0) return
    if (symbolic%front_parent(child) /= parent) then
      ierr = -3224
      return
    end if
    if (.not. allocated(handle%fronts(child)%contribution)) then
      ierr = -3223
      return
    end if

    pos_base = int(symbolic%contrib_pos_ptr(child)) - 1
    run_first = int(symbolic%contrib_run_ptr(child))
    run_last = int(symbolic%contrib_run_ptr(child + 1)) - 1
    do j = 1, child_nupd
      parent_j = symbolic%contrib_parent_pos(pos_base + j)
      if (parent_j <= 0) then
        ierr = -3222
        return
      end if
      do run = run_first, run_last
        call add_front_column_run(parent_data, &
            symbolic%contrib_run_parent_first(run), parent_j, &
            symbolic%contrib_run_len(run), &
            handle%fronts(child)%contribution(symbolic%contrib_run_first(run), j))
      end do
    end do

    deallocate(handle%fronts(child)%contribution)
  end subroutine assemble_child_contribution

  subroutine allocate_front_storage(front_data, fs, npiv, nupd, ierr)
    type(d_lu_only_numeric_front), intent(inout) :: front_data
    integer, intent(in) :: fs
    integer, intent(in) :: npiv
    integer, intent(in) :: nupd
    integer, intent(out) :: ierr

    integer :: alloc_stat

    ierr = 0
    front_data%front_size = fs
    front_data%pivot_size = npiv
    front_data%update_size = nupd

    allocate(front_data%factor(max(1, fs), max(1, npiv)), stat=alloc_stat)
    if (alloc_stat /= 0) then
      ierr = -3241
      return
    end if
    call dlaset('A', max(1, fs), max(1, npiv), 0.0_real64, 0.0_real64, &
        front_data%factor, max(1, fs))

    if (nupd > 0) then
      allocate(front_data%upper_update(max(1, npiv), nupd), &
          front_data%contribution(nupd, nupd), stat=alloc_stat)
      if (alloc_stat /= 0) then
        ierr = -3242
        return
      end if
      call dlaset('A', max(1, npiv), nupd, 0.0_real64, 0.0_real64, &
          front_data%upper_update, max(1, npiv))
      call dlaset('A', nupd, nupd, 0.0_real64, 0.0_real64, &
          front_data%contribution, max(1, nupd))
    end if
  end subroutine allocate_front_storage

  subroutine add_front_entry(front_data, row_pos, col_pos, value)
    type(d_lu_only_numeric_front), intent(inout) :: front_data
    integer, intent(in) :: row_pos
    integer, intent(in) :: col_pos
    real(real64), intent(in) :: value

    integer :: npiv

    npiv = front_data%pivot_size
    if (col_pos <= npiv) then
      front_data%factor(row_pos, col_pos) = front_data%factor(row_pos, col_pos) + value
    else if (row_pos <= npiv) then
      front_data%upper_update(row_pos, col_pos - npiv) = &
          front_data%upper_update(row_pos, col_pos - npiv) + value
    else
      front_data%contribution(row_pos - npiv, col_pos - npiv) = &
          front_data%contribution(row_pos - npiv, col_pos - npiv) + value
    end if
  end subroutine add_front_entry

  subroutine add_front_column_run(front_data, row_first, col_pos, len, source)
    type(d_lu_only_numeric_front), intent(inout) :: front_data
    integer, intent(in) :: row_first
    integer, intent(in) :: col_pos
    integer, intent(in) :: len
    real(real64), intent(in) :: source(*)

    integer :: npiv
    integer :: pivot_len
    integer :: update_len
    integer :: update_row
    integer :: update_col

    if (len <= 0) return
    npiv = front_data%pivot_size
    if (col_pos <= npiv) then
      call daxpy(len, 1.0_real64, source, 1, &
          front_data%factor(row_first, col_pos), 1)
      return
    end if

    update_col = col_pos - npiv
    pivot_len = 0
    if (row_first <= npiv) pivot_len = min(len, npiv - row_first + 1)
    if (pivot_len > 0) then
      call daxpy(pivot_len, 1.0_real64, source, 1, &
          front_data%upper_update(row_first, update_col), 1)
    end if

    update_len = len - pivot_len
    if (update_len > 0) then
      update_row = row_first + pivot_len - npiv
      call daxpy(update_len, 1.0_real64, source(1 + pivot_len), 1, &
          front_data%contribution(update_row, update_col), 1)
    end if
  end subroutine add_front_column_run

  subroutine factor_one_front(front_data, front, ierr)
    type(d_lu_only_numeric_front), intent(inout) :: front_data
    integer, intent(in) :: front
    integer, intent(out) :: ierr

    integer :: fs
    integer :: npiv
    integer :: nupd
    integer :: ldf
    integer :: ldu
    integer :: k
    integer :: j
    integer :: panel_end
    integer :: block_cols
    integer :: len
    integer :: trailing_pivots
    integer :: update_pivot_rows
    real(real64) :: pivot

    ierr = 0
    fs = front_data%front_size
    npiv = front_data%pivot_size
    nupd = front_data%update_size
    ldf = max(1, fs)
    ldu = max(1, npiv)
    if (npiv <= 0) return

    k = 1
    do while (k <= npiv)
      panel_end = min(npiv, k + front_block_size - 1)
      block_cols = panel_end - k + 1

      do j = k, panel_end
        pivot = front_data%factor(j, j)
        if (abs(pivot) <= 100.0_real64 * epsilon(1.0_real64)) then
          ierr = front
          return
        end if

        len = fs - j
        if (len > 0) then
          call dscal(len, 1.0_real64 / pivot, front_data%factor(j + 1, j), 1)
          if (j < panel_end) then
            call panel_column_update(front_data, j, j + 1, panel_end, len)
          end if
        end if
      end do

      trailing_pivots = npiv - panel_end
      if (trailing_pivots > 0) then
        call dtrsm('L', 'L', 'N', 'U', block_cols, trailing_pivots, 1.0_real64, &
            front_data%factor(k, k), ldf, front_data%factor(k, panel_end + 1), ldf)
      end if
      if (nupd > 0) then
        call dtrsm('L', 'L', 'N', 'U', block_cols, nupd, 1.0_real64, &
            front_data%factor(k, k), ldf, front_data%upper_update(k, 1), ldu)
      end if

      if (trailing_pivots > 0) then
        call dgemm('N', 'N', fs - panel_end, trailing_pivots, block_cols, &
            -1.0_real64, front_data%factor(panel_end + 1, k), ldf, &
            front_data%factor(k, panel_end + 1), ldf, 1.0_real64, &
            front_data%factor(panel_end + 1, panel_end + 1), ldf)
      end if
      if (nupd > 0) then
        update_pivot_rows = npiv - panel_end
        if (update_pivot_rows > 0) then
          call dgemm('N', 'N', update_pivot_rows, nupd, block_cols, &
              -1.0_real64, front_data%factor(panel_end + 1, k), ldf, &
              front_data%upper_update(k, 1), ldu, 1.0_real64, &
              front_data%upper_update(panel_end + 1, 1), ldu)
        end if
        call dgemm('N', 'N', nupd, nupd, block_cols, -1.0_real64, &
            front_data%factor(npiv + 1, k), ldf, front_data%upper_update(k, 1), ldu, &
            1.0_real64, front_data%contribution, max(1, nupd))
      end if

      k = panel_end + 1
    end do
  end subroutine factor_one_front

  subroutine panel_column_update(front_data, pivot_col, first_col, last_col, len)
    type(d_lu_only_numeric_front), intent(inout) :: front_data
    integer, intent(in) :: pivot_col
    integer, intent(in) :: first_col
    integer, intent(in) :: last_col
    integer, intent(in) :: len

    integer :: col

    do col = first_col, last_col
      if (front_data%factor(pivot_col, col) /= 0.0_real64) then
        call daxpy(len, -front_data%factor(pivot_col, col), &
            front_data%factor(pivot_col + 1, pivot_col), 1, &
            front_data%factor(pivot_col + 1, col), 1)
      end if
    end do
  end subroutine panel_column_update
  subroutine assemble_band_numeric(handle, values, ierr)
    type(monolis_mf_handle), intent(inout) :: handle
    real(real64), intent(in) :: values(:)
    integer, intent(out) :: ierr

    integer :: row
    integer :: entry
    integer :: col
    integer :: band_row

    ierr = 0
    call dlaset('A', handle%ldlu, handle%n, 0.0_real64, 0.0_real64, &
        handle%lu_band, handle%ldlu)
    do row = 1, handle%n
      do entry = handle%row_ptr(row), handle%row_ptr(row + 1) - 1
        col = handle%col_ind(entry)
        band_row = handle%ku + 1 + row - col
        if (band_row < 1 .or. band_row > handle%ldlu) then
          ierr = -3101
          return
        end if
        handle%lu_band(band_row, col) = handle%lu_band(band_row, col) + values(entry)
      end do
    end do
  end subroutine assemble_band_numeric

  subroutine factor_band_front_numeric(handle, ierr)
    type(monolis_mf_handle), intent(inout) :: handle
    integer, intent(out) :: ierr

    integer :: pivot_col
    integer :: row
    integer :: col
    integer :: row_end
    integer :: col_end
    integer :: row_start
    integer :: update_len
    integer :: src_band_row
    integer :: dst_band_row
    real(real64) :: pivot
    real(real64) :: upper_value

    ierr = 0
    do pivot_col = 1, handle%n
      pivot = band_value(handle, pivot_col, pivot_col)
      if (abs(pivot) <= 100.0_real64 * epsilon(1.0_real64)) then
        ierr = pivot_col
        return
      end if

      row_end = min(handle%n, pivot_col + handle%kl)
      col_end = min(handle%n, pivot_col + handle%ku)
      update_len = row_end - pivot_col

      if (update_len > 0) then
        call dscal(update_len, 1.0_real64 / pivot, &
            handle%lu_band(handle%ku + 2, pivot_col), 1)
      end if

      do col = pivot_col + 1, col_end
        upper_value = band_value(handle, pivot_col, col)
        if (upper_value == 0.0_real64) cycle

        row_start = max(pivot_col + 1, col - handle%ku)
        update_len = row_end - row_start + 1
        if (update_len <= 0) cycle

        src_band_row = handle%ku + 1 + row_start - pivot_col
        dst_band_row = handle%ku + 1 + row_start - col
        call daxpy(update_len, -upper_value, handle%lu_band(src_band_row, pivot_col), 1, &
            handle%lu_band(dst_band_row, col), 1)
      end do
    end do
  end subroutine factor_band_front_numeric
  logical function in_band(handle, row, col)
    type(monolis_mf_handle), intent(in) :: handle
    integer, intent(in) :: row
    integer, intent(in) :: col

    integer :: band_row

    band_row = handle%ku + 1 + row - col
    in_band = band_row >= 1 .and. band_row <= handle%ldlu
  end function in_band

  real(real64) function band_value(handle, row, col) result(value)
    type(monolis_mf_handle), intent(in) :: handle
    integer, intent(in) :: row
    integer, intent(in) :: col

    integer :: band_row

    band_row = handle%ku + 1 + row - col
    if (band_row < 1 .or. band_row > handle%ldlu) then
      value = 0.0_real64
    else
      value = handle%lu_band(band_row, col)
    end if
  end function band_value

  subroutine set_band_value(handle, row, col, value)
    type(monolis_mf_handle), intent(inout) :: handle
    integer, intent(in) :: row
    integer, intent(in) :: col
    real(real64), intent(in) :: value

    integer :: band_row

    band_row = handle%ku + 1 + row - col
    if (band_row >= 1 .and. band_row <= handle%ldlu) then
      handle%lu_band(band_row, col) = value
    end if
  end subroutine set_band_value

end module monolis_mf_factor_module
