!> LU 解析フェーズ
!>
!> 元 CSR 行列から、フィル削減順序付け（PORD）→ フロント木構築
!> → 元非ゼロ要素のフロント割当 → 子フロント寄与の親への対応表
!> までを行い、すべての結果を [[monolis_mat_lu]] に書き込む。
module mod_monolis_fact_analysis
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_pord_const
  use mod_monolis_pord_types
  use mod_monolis_pord_graph
  use mod_monolis_pord_tree
  use mod_monolis_pord_ordering

  implicit none

  private
  public :: monolis_fact_analysis

contains

  !> @ingroup fact
  !> 解析フェーズのエントリ
  subroutine monolis_fact_analysis(monoMAT, lu)
    implicit none
    !> [in] 行列構造体
    type(monolis_mat),    intent(in)    :: monoMAT
    !> [in,out] LU 分解構造体
    type(monolis_mat_lu), intent(inout) :: lu

    integer(kint) :: n
    integer(kint), allocatable :: row_ptr(:)
    integer(kint), allocatable :: col_ind(:)

    call monolis_std_debug_log_header("monolis_fact_analysis")

    call monolis_mat_finalize_LU(lu)

    n = monoMAT%N
    lu%N = n
    lu%NDOF = monoMAT%NDOF
    if (n <= 0) return

    call build_local_csr(monoMAT, n, row_ptr, col_ind)

    call build_symmetric_pattern(lu, n, row_ptr, col_ind)

    call build_pord_front_structure(lu)

    call build_original_front_entries(lu, n, row_ptr, col_ind)

    call build_child_contribution_maps(lu)

    call monolis_dealloc_I_1d(row_ptr)
    call monolis_dealloc_I_1d(col_ind)

    lu%analyzed = .true.
    lu%factorized = .false.
  end subroutine monolis_fact_analysis

  !> monoMAT の CSR から 1-based の row_ptr/col_ind を作る
  subroutine build_local_csr(monoMAT, n, row_ptr, col_ind)
    implicit none
    type(monolis_mat), intent(in) :: monoMAT
    integer(kint),     intent(in) :: n
    integer(kint), allocatable, intent(out) :: row_ptr(:)
    integer(kint), allocatable, intent(out) :: col_ind(:)
    integer(kint) :: i, nz

    nz = monoMAT%CSR%index(n + 1)
    call monolis_alloc_I_1d(row_ptr, n + 1)
    call monolis_alloc_I_1d(col_ind, max(1, nz))

    do i = 1, n + 1
      row_ptr(i) = monoMAT%CSR%index(i) + 1
    end do
    do i = 1, nz
      col_ind(i) = monoMAT%CSR%item(i)
    end do
  end subroutine build_local_csr

  !> 対称化パターンを CSR で構築（1-based）
  subroutine build_symmetric_pattern(lu, n, row_ptr, col_ind)
    implicit none
    type(monolis_mat_lu), intent(inout) :: lu
    integer(kint), intent(in) :: n
    integer(kint), intent(in) :: row_ptr(:)
    integer(kint), intent(in) :: col_ind(:)

    integer(kint) :: row, entry, col, pos, write_pos, cnt, prev, raw_nnz
    integer(kint), allocatable :: counts(:), raw_ptr(:), next_pos(:), raw_cols(:)
    integer(kint), allocatable :: compact_cols(:)

    call monolis_alloc_I_1d(counts, n)

    do row = 1, n
      do entry = row_ptr(row), row_ptr(row + 1) - 1
        col = col_ind(entry)
        counts(row) = counts(row) + 1
        if (col /= row) counts(col) = counts(col) + 1
      end do
    end do

    call monolis_alloc_I_1d(raw_ptr, n + 1)
    call monolis_alloc_I_1d(next_pos, n)

    raw_ptr(1) = 1
    do row = 1, n
      raw_ptr(row + 1) = raw_ptr(row) + counts(row)
    end do
    raw_nnz = raw_ptr(n + 1) - 1

    call monolis_alloc_I_1d(raw_cols, max(1, raw_nnz))

    next_pos(1:n) = raw_ptr(1:n)
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

    call monolis_alloc_I_1d(lu%sym_row_ptr, n + 1)
    call monolis_alloc_I_1d(lu%sym_col_ind, max(1, raw_nnz))

    write_pos = 1
    lu%sym_row_ptr(1) = 1
    do row = 1, n
      call sort_int_range(raw_cols, raw_ptr(row), raw_ptr(row + 1) - 1)
      prev = 0
      cnt = 0
      do entry = raw_ptr(row), raw_ptr(row + 1) - 1
        if (cnt == 0 .or. raw_cols(entry) /= prev) then
          lu%sym_col_ind(write_pos) = raw_cols(entry)
          write_pos = write_pos + 1
          cnt = cnt + 1
          prev = raw_cols(entry)
        end if
      end do
      lu%sym_row_ptr(row + 1) = write_pos
    end do

    lu%sym_nnz = write_pos - 1
    call monolis_alloc_I_1d(compact_cols, lu%sym_nnz)
    compact_cols(1:lu%sym_nnz) = lu%sym_col_ind(1:lu%sym_nnz)
    call move_alloc(compact_cols, lu%sym_col_ind)

    call monolis_dealloc_I_1d(counts)
    call monolis_dealloc_I_1d(raw_ptr)
    call monolis_dealloc_I_1d(next_pos)
    call monolis_dealloc_I_1d(raw_cols)
  end subroutine build_symmetric_pattern

  !> PORD を呼び出し、perm/iperm/front_*/super_* を構築
  subroutine build_pord_front_structure(lu)
    implicit none
    type(monolis_mat_lu), intent(inout) :: lu

    type(graph_t)    :: G
    type(elimtree_t) :: T
    integer(kint) :: opts_dummy(1)
    integer(kint) :: n, nfronts, nz_off, row, entry, col, cnt, pos
    integer(kint) :: front, child, total_front_vars, pivot, first_col
    integer(kint) :: post_pos, new_col, i
    integer(kint), allocatable :: marker(:), work(:)

    n = lu%N

    nz_off = 0
    do row = 1, n
      do entry = lu%sym_row_ptr(row), lu%sym_row_ptr(row + 1) - 1
        if (lu%sym_col_ind(entry) /= row) nz_off = nz_off + 1
      end do
    end do

    call newGraph(G, n, nz_off)
    G%xadj(0) = 0
    do row = 1, n
      cnt = 0
      do entry = lu%sym_row_ptr(row), lu%sym_row_ptr(row + 1) - 1
        if (lu%sym_col_ind(entry) /= row) cnt = cnt + 1
      end do
      G%xadj(row) = G%xadj(row - 1) + cnt
    end do
    do row = 1, n
      pos = G%xadj(row - 1)
      do entry = lu%sym_row_ptr(row), lu%sym_row_ptr(row + 1) - 1
        col = lu%sym_col_ind(entry)
        if (col /= row) then
          G%adjncy(pos) = col - 1
          pos = pos + 1
        end if
      end do
    end do

    opts_dummy = 0
    call monolis_pord_ordering(G, opts_dummy, .true., T)

    nfronts = T%nfronts
    lu%nsuper = nfronts
    lu%nfronts = nfronts

    call monolis_alloc_I_1d(lu%perm,  n)
    call monolis_alloc_I_1d(lu%iperm, n)
    !> permFromElimTree は old -> new を返す（1-based）。
    !> perm として保存し、iperm はその逆。
    call monolis_pord_perm_from_elimtree(T, n, lu%perm)
    do i = 1, n
      lu%iperm(lu%perm(i)) = i
    end do

    call monolis_alloc_I_1d(lu%vertex_front,        n)
    call monolis_alloc_I_1d(lu%column_to_super,     n)
    call monolis_alloc_I_1d(lu%super_start,         nfronts)
    call monolis_alloc_I_1d(lu%super_end,           nfronts)
    call monolis_alloc_I_1d(lu%super_parent,        nfronts)
    call monolis_alloc_I_1d(lu%front_parent,        nfronts)
    call monolis_alloc_I_1d(lu%front_first_child,   nfronts)
    call monolis_alloc_I_1d(lu%front_next_sibling,  nfronts)
    call monolis_alloc_I_1d(lu%front_postorder,     nfronts)
    call monolis_alloc_I_1d(lu%front_size,          nfronts)
    call monolis_alloc_I_1d(lu%front_pivot_size,    nfronts)
    call monolis_alloc_I_1d(lu%front_update_size,   nfronts)
    call monolis_alloc_I_1d(lu%front_ptr,           nfronts + 1)

    do row = 1, n
      lu%vertex_front(row) = T%vtx2front(row - 1) + 1
    end do

    lu%super_start(:) = n + 1
    lu%super_end(:) = 0
    lu%column_to_super(:) = 0
    lu%front_first_child(:) = 0
    lu%front_next_sibling(:) = 0
    lu%front_postorder(:) = 0
    lu%l_pattern_nnz = 0
    lu%max_front_size = 0

    lu%front_ptr(1) = 1
    total_front_vars = 0
    do front = 1, nfronts
      lu%front_pivot_size(front)  = T%ncolfactor(front - 1)
      lu%front_update_size(front) = T%ncolupdate(front - 1)
      lu%front_size(front) = lu%front_pivot_size(front) + lu%front_update_size(front)
      lu%max_front_size = max(lu%max_front_size, lu%front_size(front))
      lu%l_pattern_nnz = lu%l_pattern_nnz + &
          lu%front_pivot_size(front) * lu%front_size(front)
      total_front_vars = total_front_vars + lu%front_size(front)
      lu%front_ptr(front + 1) = lu%front_ptr(front) + lu%front_size(front)

      if (T%parent(front - 1) == -1) then
        lu%front_parent(front) = 0
      else
        lu%front_parent(front) = T%parent(front - 1) + 1
      end if
      lu%super_parent(front) = lu%front_parent(front)
    end do

    call monolis_alloc_I_1d(lu%front_ind, max(1, total_front_vars))

    do row = 1, n
      new_col = lu%perm(row)
      front = lu%vertex_front(row)
      lu%column_to_super(new_col) = front
      lu%super_start(front) = min(lu%super_start(front), new_col)
      lu%super_end(front)   = max(lu%super_end(front),   new_col)
    end do

    !> firstchild / next_sibling
    do front = nfronts, 1, -1
      child = lu%front_parent(front)
      if (child > 0) then
        lu%front_next_sibling(front) = lu%front_first_child(child)
        lu%front_first_child(child) = front
      end if
    end do

    !> postorder
    call build_front_postorder(lu, post_pos)
    if (post_pos /= nfronts) then
      call monolis_std_error_string("monolis_fact_analysis: postorder count mismatch")
      call monolis_std_error_stop()
    end if

    !> front_ind を埋める
    call monolis_alloc_I_1d(marker, n)
    call monolis_alloc_I_1d(work, n)
    do post_pos = 1, nfronts
      front = lu%front_postorder(post_pos)
      cnt = 0
      first_col = lu%super_start(front)

      do pivot = 0, lu%front_pivot_size(front) - 1
        cnt = cnt + 1
        work(cnt) = first_col + pivot
        marker(first_col + pivot) = front
      end do

      child = lu%front_first_child(front)
      do while (child /= 0)
        do pos = lu%front_ptr(child), lu%front_ptr(child + 1) - 1
          col = lu%front_ind(pos)
          if (col > first_col .and. marker(col) /= front) then
            cnt = cnt + 1
            work(cnt) = col
            marker(col) = front
          end if
        end do
        child = lu%front_next_sibling(child)
      end do

      do pivot = 0, lu%front_pivot_size(front) - 1
        row = lu%iperm(first_col + pivot)
        do entry = lu%sym_row_ptr(row), lu%sym_row_ptr(row + 1) - 1
          col = lu%perm(lu%sym_col_ind(entry))
          if (col > first_col .and. marker(col) /= front) then
            cnt = cnt + 1
            work(cnt) = col
            marker(col) = front
          end if
        end do
      end do

      if (cnt /= lu%front_size(front)) then
        call monolis_std_error_string("monolis_fact_analysis: front size mismatch")
        call monolis_std_error_stop()
      end if

      call sort_int_range(work, 1, cnt)
      pos = lu%front_ptr(front)
      lu%front_ind(pos:pos + cnt - 1) = work(1:cnt)
    end do

    call monolis_dealloc_I_1d(marker)
    call monolis_dealloc_I_1d(work)
    call freeElimTree(T)
    call freeGraph(G)
  end subroutine build_pord_front_structure

  !> 各非ゼロ要素を所属フロントに対応付け、orig_* を構築
  subroutine build_original_front_entries(lu, n, row_ptr, col_ind)
    implicit none
    type(monolis_mat_lu), intent(inout) :: lu
    integer(kint), intent(in) :: n
    integer(kint), intent(in) :: row_ptr(:)
    integer(kint), intent(in) :: col_ind(:)

    integer(kint) :: nfronts, row, entry, col, new_row, new_col, first_col
    integer(kint) :: front, pos, first_pos, last_pos, row_local, col_local
    integer(kint) :: total_nnz
    integer(kint), allocatable :: counts(:), next_pos(:), local_pos(:)

    nfronts = lu%nfronts
    total_nnz = row_ptr(n + 1) - 1

    call monolis_alloc_I_1d(counts, nfronts)
    call monolis_alloc_I_1d(lu%orig_ptr, nfronts + 1)

    do row = 1, n
      new_row = lu%perm(row)
      do entry = row_ptr(row), row_ptr(row + 1) - 1
        new_col = lu%perm(col_ind(entry))
        first_col = min(new_row, new_col)
        front = lu%column_to_super(first_col)
        counts(front) = counts(front) + 1
      end do
    end do

    lu%orig_ptr(1) = 1
    do front = 1, nfronts
      lu%orig_ptr(front + 1) = lu%orig_ptr(front) + counts(front)
    end do

    call monolis_alloc_I_1d(lu%orig_row_pos, max(1, total_nnz))
    call monolis_alloc_I_1d(lu%orig_col_pos, max(1, total_nnz))
    call monolis_alloc_I_1d(lu%orig_entry,   max(1, total_nnz))
    call monolis_alloc_I_1d(next_pos, nfronts)

    next_pos(1:nfronts) = lu%orig_ptr(1:nfronts)
    do row = 1, n
      new_row = lu%perm(row)
      do entry = row_ptr(row), row_ptr(row + 1) - 1
        col = col_ind(entry)
        new_col = lu%perm(col)
        first_col = min(new_row, new_col)
        front = lu%column_to_super(first_col)
        pos = next_pos(front)

        lu%orig_row_pos(pos) = new_row
        lu%orig_col_pos(pos) = new_col
        lu%orig_entry(pos)   = entry
        next_pos(front) = pos + 1
      end do
    end do

    !> グローバル列番号 → フロント内ローカル位置に変換
    call monolis_alloc_I_1d(local_pos, n)
    do front = 1, nfronts
      first_pos = lu%front_ptr(front)
      last_pos  = lu%front_ptr(front + 1) - 1
      do pos = first_pos, last_pos
        local_pos(lu%front_ind(pos)) = pos - first_pos + 1
      end do

      do pos = lu%orig_ptr(front), lu%orig_ptr(front + 1) - 1
        row_local = local_pos(lu%orig_row_pos(pos))
        col_local = local_pos(lu%orig_col_pos(pos))
        if (row_local <= 0 .or. col_local <= 0) then
          call monolis_std_error_string("monolis_fact_analysis: local position not found")
          call monolis_std_error_stop()
        end if
        lu%orig_row_pos(pos) = row_local
        lu%orig_col_pos(pos) = col_local
      end do

      do pos = first_pos, last_pos
        local_pos(lu%front_ind(pos)) = 0
      end do
    end do

    call monolis_dealloc_I_1d(counts)
    call monolis_dealloc_I_1d(next_pos)
    call monolis_dealloc_I_1d(local_pos)
  end subroutine build_original_front_entries

  !> 子フロント寄与の親への対応表を構築
  subroutine build_child_contribution_maps(lu)
    implicit none
    type(monolis_mat_lu), intent(inout) :: lu

    integer(kint) :: nfronts, front, parent, npiv, nupd, i, run, pos
    integer(kint) :: child_pos, run_count, max_update_size
    integer(kint) :: total_update, total_runs
    integer(kint), allocatable :: parent_pos(:)

    nfronts = lu%nfronts
    total_update = 0
    total_runs = 0
    max_update_size = 0

    do front = 1, nfronts
      nupd = lu%front_update_size(front)
      max_update_size = max(max_update_size, nupd)
      total_update = total_update + nupd
    end do

    call monolis_alloc_I_1d(lu%contrib_pos_ptr, nfronts + 1)
    call monolis_alloc_I_1d(lu%contrib_run_ptr, nfronts + 1)
    call monolis_alloc_I_1d(parent_pos, max(1, max_update_size))

    lu%contrib_pos_ptr(1) = 1
    lu%contrib_run_ptr(1) = 1
    do front = 1, nfronts
      parent = lu%front_parent(front)
      npiv = lu%front_pivot_size(front)
      nupd = lu%front_update_size(front)
      lu%contrib_pos_ptr(front + 1) = lu%contrib_pos_ptr(front) + nupd

      run_count = 0
      if (parent > 0 .and. nupd > 0) then
        do i = 1, nupd
          child_pos = lu%front_ptr(front) + npiv + i - 1
          parent_pos(i) = front_local_position(lu, parent, lu%front_ind(child_pos))
          if (parent_pos(i) <= 0) then
            call monolis_std_error_string("monolis_fact_analysis: parent_pos not found")
            call monolis_std_error_stop()
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

      total_runs = total_runs + run_count
      lu%contrib_run_ptr(front + 1) = lu%contrib_run_ptr(front) + run_count
    end do

    call monolis_alloc_I_1d(lu%contrib_parent_pos,       max(1, total_update))
    call monolis_alloc_I_1d(lu%contrib_run_first,        max(1, total_runs))
    call monolis_alloc_I_1d(lu%contrib_run_len,          max(1, total_runs))
    call monolis_alloc_I_1d(lu%contrib_run_parent_first, max(1, total_runs))

    do front = 1, nfronts
      parent = lu%front_parent(front)
      npiv = lu%front_pivot_size(front)
      nupd = lu%front_update_size(front)
      pos = lu%contrib_pos_ptr(front)
      run = lu%contrib_run_ptr(front)

      if (parent > 0 .and. nupd > 0) then
        do i = 1, nupd
          child_pos = lu%front_ptr(front) + npiv + i - 1
          parent_pos(i) = front_local_position(lu, parent, lu%front_ind(child_pos))
          lu%contrib_parent_pos(pos + i - 1) = parent_pos(i)
        end do

        i = 1
        do while (i <= nupd)
          lu%contrib_run_first(run) = i
          lu%contrib_run_parent_first(run) = parent_pos(i)
          do while (i < nupd)
            if (parent_pos(i + 1) /= parent_pos(i) + 1) exit
            i = i + 1
          end do
          lu%contrib_run_len(run) = i - lu%contrib_run_first(run) + 1
          run = run + 1
          i = i + 1
        end do
      end if
    end do

    call monolis_dealloc_I_1d(parent_pos)
  end subroutine build_child_contribution_maps

  !> フロント木の後順走査
  subroutine build_front_postorder(lu, post_count)
    implicit none
    type(monolis_mat_lu), intent(inout) :: lu
    integer(kint), intent(out) :: post_count

    integer(kint) :: nfronts, root, top, child
    integer(kint), allocatable :: stack_node(:), stack_next_child(:)

    nfronts = lu%nfronts
    post_count = 0
    if (nfronts <= 0) return

    call monolis_alloc_I_1d(stack_node, nfronts)
    call monolis_alloc_I_1d(stack_next_child, nfronts)

    do root = 1, nfronts
      if (lu%front_parent(root) /= 0) cycle

      top = 1
      stack_node(top) = root
      stack_next_child(top) = lu%front_first_child(root)
      do while (top > 0)
        child = stack_next_child(top)
        if (child /= 0) then
          stack_next_child(top) = lu%front_next_sibling(child)
          top = top + 1
          stack_node(top) = child
          stack_next_child(top) = lu%front_first_child(child)
        else
          post_count = post_count + 1
          lu%front_postorder(post_count) = stack_node(top)
          top = top - 1
        end if
      end do
    end do

    call monolis_dealloc_I_1d(stack_node)
    call monolis_dealloc_I_1d(stack_next_child)
  end subroutine build_front_postorder

  !> フロント内の global col → local position（1-based）を二分探索で求める
  integer(kint) function front_local_position(lu, front, col) result(local_pos)
    implicit none
    type(monolis_mat_lu), intent(in) :: lu
    integer(kint), intent(in) :: front, col

    integer(kint) :: left, right, mid, base, val

    local_pos = 0
    left  = lu%front_ptr(front)
    right = lu%front_ptr(front + 1) - 1
    base = left - 1
    do while (left <= right)
      mid = (left + right) / 2
      val = lu%front_ind(mid)
      if (val == col) then
        local_pos = mid - base
        return
      else if (val < col) then
        left = mid + 1
      else
        right = mid - 1
      end if
    end do
  end function front_local_position

  !> Shell sort（整数昇順）
  subroutine sort_int_range(values, first, last)
    implicit none
    integer(kint), intent(inout) :: values(:)
    integer(kint), intent(in) :: first, last

    integer(kint) :: gap, i, j, temp, len

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

end module mod_monolis_fact_analysis
