module monolis_mf_reorder_module
  use iso_c_binding, only: c_int
  use iso_fortran_env, only: int64
  use monolis_mf_types_module, only: d_lu_only_symbolic_analysis
  use monolis_mf_reorder_interface_module, only: dlu_pord_ordering
  implicit none

  private
  public :: build_front_postorder
  public :: front_local_position
  public :: build_pord_front_structure
  public :: sort_int_range

contains

  subroutine build_front_postorder(symbolic, post_count)
    type(d_lu_only_symbolic_analysis), intent(inout) :: symbolic
    integer, intent(out) :: post_count

    integer :: nfronts
    integer :: root
    integer :: top
    integer :: child
    integer :: alloc_stat
    integer, allocatable :: stack_node(:)
    integer, allocatable :: stack_next_child(:)

    nfronts = symbolic%nfronts
    post_count = 0
    if (nfronts <= 0) return

    allocate(stack_node(nfronts), stack_next_child(nfronts), stat=alloc_stat)
    if (alloc_stat /= 0) return

    do root = 1, nfronts
      if (symbolic%front_parent(root) /= 0) cycle

      top = 1
      stack_node(top) = root
      stack_next_child(top) = symbolic%front_first_child(root)
      do while (top > 0)
        child = stack_next_child(top)
        if (child /= 0) then
          stack_next_child(top) = symbolic%front_next_sibling(child)
          top = top + 1
          stack_node(top) = child
          stack_next_child(top) = symbolic%front_first_child(child)
        else
          post_count = post_count + 1
          symbolic%front_postorder(post_count) = stack_node(top)
          top = top - 1
        end if
      end do
    end do

    deallocate(stack_node, stack_next_child)
  end subroutine build_front_postorder

  integer function front_local_position(symbolic, front, col) result(local_pos)
    type(d_lu_only_symbolic_analysis), intent(in) :: symbolic
    integer, intent(in) :: front
    integer, intent(in) :: col

    integer :: left
    integer :: right
    integer :: mid
    integer :: base
    integer :: value

    local_pos = 0
    left = int(symbolic%front_ptr(front))
    right = int(symbolic%front_ptr(front + 1)) - 1
    base = left - 1
    do while (left <= right)
      mid = (left + right) / 2
      value = symbolic%front_ind(mid)
      if (value == col) then
        local_pos = mid - base
        return
      else if (value < col) then
        left = mid + 1
      else
        right = mid - 1
      end if
    end do
  end function front_local_position

  subroutine build_pord_front_structure(symbolic, ierr)
    type(d_lu_only_symbolic_analysis), intent(inout) :: symbolic
    integer, intent(out) :: ierr

    integer :: n
    integer :: row
    integer :: entry
    integer :: col
    integer :: pos
    integer :: front
    integer :: new_col
    integer :: nfronts
    integer :: pivot
    integer :: child
    integer :: count
    integer :: first_col
    integer :: pord_status
    integer :: alloc_stat
    integer :: graph_nnz
    integer :: post_pos
    integer :: total_front_vars_i
    integer(int64) :: front_entries
    integer(int64) :: total_front_vars
    integer, allocatable :: counts(:)
    integer, allocatable :: next_pos(:)
    integer, allocatable :: marker(:)
    integer, allocatable :: work(:)
    integer(c_int), allocatable :: pord_xadj(:)
    integer(c_int), allocatable :: pord_adjncy(:)
    integer(c_int), allocatable :: pord_perm(:)
    integer(c_int), allocatable :: pord_invp(:)
    integer(c_int), allocatable :: pord_parent(:)
    integer(c_int), allocatable :: pord_size(:)
    integer(c_int), allocatable :: pord_update(:)
    integer(c_int), allocatable :: pord_vertex_front(:)
    integer(c_int) :: c_nfronts

    ierr = 0
    n = symbolic%n

    allocate(counts(n), stat=alloc_stat)
    if (alloc_stat /= 0) then
      ierr = -2851
      return
    end if
    counts = 0

    do row = 1, n
      do entry = symbolic%sym_row_ptr(row), symbolic%sym_row_ptr(row + 1) - 1
        col = symbolic%sym_col_ind(entry)
        if (col /= row) counts(row) = counts(row) + 1
      end do
    end do

    allocate(pord_xadj(n + 1), next_pos(n), stat=alloc_stat)
    if (alloc_stat /= 0) then
      ierr = -2852
      deallocate(counts)
      return
    end if

    pord_xadj(1) = 1
    do row = 1, n
      pord_xadj(row + 1) = pord_xadj(row) + counts(row)
    end do
    graph_nnz = pord_xadj(n + 1) - 1

    allocate(pord_adjncy(max(1, graph_nnz)), pord_perm(n), pord_invp(n), &
        pord_parent(n), pord_size(n), pord_update(n), pord_vertex_front(n), &
        stat=alloc_stat)
    if (alloc_stat /= 0) then
      ierr = -2853
      deallocate(counts, next_pos, pord_xadj)
      return
    end if

    next_pos = int(pord_xadj(1:n))
    do row = 1, n
      do entry = symbolic%sym_row_ptr(row), symbolic%sym_row_ptr(row + 1) - 1
        col = symbolic%sym_col_ind(entry)
        if (col == row) cycle
        pos = next_pos(row)
        pord_adjncy(pos) = col
        next_pos(row) = pos + 1
      end do
    end do

    pord_status = int(dlu_pord_ordering(int(n, c_int), int(graph_nnz, c_int), &
        pord_xadj, pord_adjncy, pord_perm, pord_invp, c_nfronts, &
        pord_parent, pord_size, pord_update, pord_vertex_front))
    if (pord_status /= 0) then
      ierr = -2860 + pord_status
      deallocate(counts, next_pos, pord_xadj, pord_adjncy, pord_perm, pord_invp, &
          pord_parent, pord_size, pord_update, pord_vertex_front)
      return
    end if

    nfronts = int(c_nfronts)
    symbolic%nsuper = nfronts
    symbolic%nfronts = nfronts

    allocate(symbolic%perm(n), symbolic%invp(n), symbolic%vertex_front(n), &
        symbolic%super_start(nfronts), symbolic%super_end(nfronts), &
        symbolic%column_to_super(n), symbolic%super_parent(nfronts), &
        symbolic%front_ptr(nfronts + 1), symbolic%front_parent(nfronts), &
        symbolic%front_first_child(nfronts), symbolic%front_next_sibling(nfronts), &
        symbolic%front_postorder(nfronts), &
        symbolic%front_size(nfronts), symbolic%front_pivot_size(nfronts), &
        symbolic%front_update_size(nfronts), stat=alloc_stat)
    if (alloc_stat /= 0) then
      ierr = -2854
      deallocate(counts, next_pos, pord_xadj, pord_adjncy, pord_perm, pord_invp, &
          pord_parent, pord_size, pord_update, pord_vertex_front)
      return
    end if

    symbolic%perm = int(pord_perm)
    symbolic%invp = int(pord_invp)
    symbolic%vertex_front = int(pord_vertex_front)
    symbolic%super_start = n + 1
    symbolic%super_end = 0
    symbolic%column_to_super = 0
    symbolic%front_first_child = 0
    symbolic%front_next_sibling = 0
    symbolic%front_postorder = 0
    symbolic%l_pattern_nnz = 0_int64
    symbolic%max_front_size = 0

    symbolic%front_ptr(1) = 1_int64
    total_front_vars = 0_int64
    do front = 1, nfronts
      symbolic%front_parent(front) = int(pord_parent(front))
      symbolic%super_parent(front) = symbolic%front_parent(front)
      symbolic%front_pivot_size(front) = int(pord_size(front))
      symbolic%front_update_size(front) = int(pord_update(front))
      symbolic%front_size(front) = symbolic%front_pivot_size(front) + symbolic%front_update_size(front)
      symbolic%max_front_size = max(symbolic%max_front_size, symbolic%front_size(front))
      symbolic%l_pattern_nnz = symbolic%l_pattern_nnz + &
          int(symbolic%front_pivot_size(front), int64) * int(symbolic%front_size(front), int64)
      total_front_vars = total_front_vars + int(symbolic%front_size(front), int64)
      symbolic%front_ptr(front + 1) = symbolic%front_ptr(front) + &
          int(symbolic%front_size(front), int64)
    end do

    if (total_front_vars > int(huge(1), int64)) then
      ierr = -2855
      deallocate(counts, next_pos, pord_xadj, pord_adjncy, pord_perm, pord_invp, &
          pord_parent, pord_size, pord_update, pord_vertex_front)
      return
    end if
    total_front_vars_i = int(total_front_vars)

    allocate(symbolic%front_ind(max(1, total_front_vars_i)), marker(n), work(n), &
        stat=alloc_stat)
    if (alloc_stat /= 0) then
      ierr = -2856
      deallocate(counts, next_pos, pord_xadj, pord_adjncy, pord_perm, pord_invp, &
          pord_parent, pord_size, pord_update, pord_vertex_front)
      return
    end if

    do row = 1, n
      new_col = symbolic%perm(row)
      front = symbolic%vertex_front(row)
      symbolic%column_to_super(new_col) = front
      symbolic%super_start(front) = min(symbolic%super_start(front), new_col)
      symbolic%super_end(front) = max(symbolic%super_end(front), new_col)
    end do

    do front = nfronts, 1, -1
      child = symbolic%front_parent(front)
      if (child > 0) then
        symbolic%front_next_sibling(front) = symbolic%front_first_child(child)
        symbolic%front_first_child(child) = front
      end if
    end do

    call build_front_postorder(symbolic, post_pos)
    if (post_pos /= nfronts) then
      ierr = -2857
      deallocate(counts, next_pos, pord_xadj, pord_adjncy, pord_perm, pord_invp, &
          pord_parent, pord_size, pord_update, pord_vertex_front, marker, work)
      return
    end if

    marker = 0
    do post_pos = 1, nfronts
      front = symbolic%front_postorder(post_pos)
      count = 0
      first_col = symbolic%super_start(front)

      do pivot = 0, symbolic%front_pivot_size(front) - 1
        count = count + 1
        work(count) = first_col + pivot
        marker(first_col + pivot) = front
      end do

      child = symbolic%front_first_child(front)
      do while (child /= 0)
        do pos = int(symbolic%front_ptr(child)), int(symbolic%front_ptr(child + 1)) - 1
          col = symbolic%front_ind(pos)
          if (col > first_col .and. marker(col) /= front) then
            count = count + 1
            work(count) = col
            marker(col) = front
          end if
        end do
        child = symbolic%front_next_sibling(child)
      end do

      do pivot = 0, symbolic%front_pivot_size(front) - 1
        row = symbolic%invp(first_col + pivot)
        do entry = symbolic%sym_row_ptr(row), symbolic%sym_row_ptr(row + 1) - 1
          col = symbolic%perm(symbolic%sym_col_ind(entry))
          if (col > first_col .and. marker(col) /= front) then
            count = count + 1
            work(count) = col
            marker(col) = front
          end if
        end do
      end do

      if (count /= symbolic%front_size(front)) then
        ierr = -2858
        deallocate(counts, next_pos, pord_xadj, pord_adjncy, pord_perm, pord_invp, &
            pord_parent, pord_size, pord_update, pord_vertex_front, marker, work)
        return
      end if

      call sort_int_range(work, 1, count)
      pos = int(symbolic%front_ptr(front))
      symbolic%front_ind(pos:pos + count - 1) = work(1:count)
    end do

    deallocate(counts, next_pos, pord_xadj, pord_adjncy, pord_perm, pord_invp, &
        pord_parent, pord_size, pord_update, pord_vertex_front, marker, work)
  end subroutine build_pord_front_structure
  subroutine sort_int_range(values, first, last)
    integer, intent(inout) :: values(:)
    integer, intent(in) :: first
    integer, intent(in) :: last

    integer :: gap
    integer :: i
    integer :: j
    integer :: temp
    integer :: len

    len = last - first + 1
    if (len <= 1) return

    gap = len / 2
    do while (gap > 0)
      do i = first + gap, last
        temp = values(i)
        j = i
        do while (j - gap >= first)
          if (values(j - gap) <= temp) exit
          values(j) = values(j - gap)
          j = j - gap
        end do
        values(j) = temp
      end do
      gap = gap / 2
    end do
  end subroutine sort_int_range

end module monolis_mf_reorder_module
