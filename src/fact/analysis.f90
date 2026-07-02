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

  !> amalgamation（スーパーノード融合）: relaxed 融合で対象とする子フロントの
  !> ピボット数上限。これ以下の小フロントを列方向に隣接する親へ融合する。
  integer(kint), parameter :: amalg_nemin = 32
  !> amalgamation: relaxed 融合で許容する padding（親と子のフロントサイズ差）上限。
  !> 余分なゼロ行の演算量を抑えるための制限。
  integer(kint), parameter :: amalg_npad = 64

contains

  !> @ingroup fact
  !> 解析フェーズのエントリ
  subroutine monolis_fact_analysis(monoMAT, lu)
    implicit none
    !> [in] 行列構造体
    type(monolis_mat),    intent(in)    :: monoMAT
    !> [in,out] LU 分解構造体
    type(monolis_mat_lu), intent(inout) :: lu

    integer(kint) :: n, ndof
    integer(kint), allocatable :: row_ptr(:)
    integer(kint), allocatable :: col_ind(:)
    !> スカラー化 CSR の各非ゼロ → monoMAT%R%A の位置（1-based）
    integer(kint), allocatable :: val_pos(:)
    !> 対称化パターン CSR（1-based、N+1）
    integer(kint), allocatable :: sym_row_ptr(:)
    !> 対称化パターン CSR（1-based、sym_nnz）
    integer(kint), allocatable :: sym_col_ind(:)
    !> 置換後列 → フロント
    integer(kint), allocatable :: column_to_super(:)

    call monolis_std_debug_log_header("monolis_fact_analysis")

    call monolis_mat_finalize_LU(lu)

    ndof = monoMAT%NDOF
    if (ndof <= 0) ndof = 1
    n = monoMAT%N * ndof
    lu%N = n
    if (n <= 0) return

    call build_local_csr(monoMAT, ndof, n, row_ptr, col_ind, val_pos)

    call build_symmetric_pattern(n, row_ptr, col_ind, sym_row_ptr, sym_col_ind)

    call build_pord_front_structure(lu, sym_row_ptr, sym_col_ind, column_to_super)

    call build_original_front_entries(lu, n, row_ptr, col_ind, val_pos, column_to_super)

    call build_child_contribution_maps(lu)

    call monolis_dealloc_I_1d(row_ptr)
    call monolis_dealloc_I_1d(col_ind)
    call monolis_dealloc_I_1d(val_pos)
    call monolis_dealloc_I_1d(sym_row_ptr)
    call monolis_dealloc_I_1d(sym_col_ind)
    call monolis_dealloc_I_1d(column_to_super)

    lu%analyzed = .true.
    lu%factorized = .false.
  end subroutine monolis_fact_analysis

  !> monoMAT のブロック CSR を NDOF 倍に展開したスカラー CSR を作る
  !>
  !> 各 NDOF×NDOF ブロックを NDOF^2 個のスカラー非ゼロに展開し、
  !> スカラー行 NDOF*(bi-1)+r、スカラー列 NDOF*(bj-1)+c の順で並べる。
  !> val_pos(scalar_entry) = NDOF^2*(block_entry-1) + NDOF*(r-1) + c
  !> （monoMAT%R%A 内での 1-based 位置）。
  subroutine build_local_csr(monoMAT, ndof, n, row_ptr, col_ind, val_pos)
    implicit none
    type(monolis_mat), intent(in) :: monoMAT
    integer(kint),     intent(in) :: ndof
    integer(kint),     intent(in) :: n
    integer(kint), allocatable, intent(out) :: row_ptr(:)
    integer(kint), allocatable, intent(out) :: col_ind(:)
    integer(kint), allocatable, intent(out) :: val_pos(:)

    integer(kint) :: nb, nz_block, ndof2, bi, bj, ii, jS, jE, r, c
    integer(kint) :: sr, pos, nb_in_row

    nb = monoMAT%N
    nz_block = monoMAT%CSR%index(nb + 1)
    ndof2 = ndof * ndof

    call monolis_alloc_I_1d(row_ptr, n + 1)
    call monolis_alloc_I_1d(col_ind, max(1, nz_block * ndof2))
    call monolis_alloc_I_1d(val_pos, max(1, nz_block * ndof2))

    row_ptr(1) = 1
    do bi = 1, nb
      nb_in_row = monoMAT%CSR%index(bi + 1) - monoMAT%CSR%index(bi)
      do r = 1, ndof
        sr = (bi - 1) * ndof + r
        row_ptr(sr + 1) = row_ptr(sr) + nb_in_row * ndof
      end do
    end do

    pos = 1
    do bi = 1, nb
      jS = monoMAT%CSR%index(bi) + 1
      jE = monoMAT%CSR%index(bi + 1)
      do r = 1, ndof
        do ii = jS, jE
          bj = monoMAT%CSR%item(ii)
          do c = 1, ndof
            col_ind(pos) = (bj - 1) * ndof + c
            val_pos(pos) = ndof2 * (ii - 1) + ndof * (r - 1) + c
            pos = pos + 1
          end do
        end do
      end do
    end do
  end subroutine build_local_csr

  !> 対称化パターンを CSR で構築（1-based）
  subroutine build_symmetric_pattern(n, row_ptr, col_ind, sym_row_ptr, sym_col_ind)
    implicit none
    integer(kint), intent(in) :: n
    integer(kint), intent(in) :: row_ptr(:)
    integer(kint), intent(in) :: col_ind(:)
    integer(kint), allocatable, intent(out) :: sym_row_ptr(:)
    integer(kint), allocatable, intent(out) :: sym_col_ind(:)

    integer(kint) :: row, entry, col, pos, write_pos, cnt, prev, raw_nnz, sym_nnz
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

    call monolis_alloc_I_1d(sym_row_ptr, n + 1)
    call monolis_alloc_I_1d(sym_col_ind, max(1, raw_nnz))

    write_pos = 1
    sym_row_ptr(1) = 1
    do row = 1, n
      call sort_int_range(raw_cols, raw_ptr(row), raw_ptr(row + 1) - 1)
      prev = 0
      cnt = 0
      do entry = raw_ptr(row), raw_ptr(row + 1) - 1
        if (cnt == 0 .or. raw_cols(entry) /= prev) then
          sym_col_ind(write_pos) = raw_cols(entry)
          write_pos = write_pos + 1
          cnt = cnt + 1
          prev = raw_cols(entry)
        end if
      end do
      sym_row_ptr(row + 1) = write_pos
    end do

    sym_nnz = write_pos - 1
    call monolis_alloc_I_1d(compact_cols, sym_nnz)
    compact_cols(1:sym_nnz) = sym_col_ind(1:sym_nnz)
    call move_alloc(compact_cols, sym_col_ind)

    call monolis_dealloc_I_1d(counts)
    call monolis_dealloc_I_1d(raw_ptr)
    call monolis_dealloc_I_1d(next_pos)
    call monolis_dealloc_I_1d(raw_cols)
  end subroutine build_symmetric_pattern

  !> PORD を呼び出し、perm/iperm/front_*/super_* を構築
  subroutine build_pord_front_structure(lu, sym_row_ptr, sym_col_ind, column_to_super)
    implicit none
    type(monolis_mat_lu), intent(inout) :: lu
    integer(kint), intent(in) :: sym_row_ptr(:)
    integer(kint), intent(in) :: sym_col_ind(:)
    integer(kint), allocatable, intent(out) :: column_to_super(:)

    type(graph_t)    :: G
    type(elimtree_t) :: T
    integer(kint) :: opts_dummy(1)
    integer(kint) :: n, nfronts, nz_off, row, entry, col, cnt, pos
    integer(kint) :: front, child, total_front_vars, pivot, first_col
    integer(kint) :: post_pos, new_col, i
    integer(kint) :: nf0, p, r0, size_f, size_p
    logical :: do_merge
    integer(kint), allocatable :: marker(:), work(:)
    integer(kint), allocatable :: vertex_front(:)
    !> amalgamation 用作業配列（元フロント単位）
    integer(kint), allocatable :: piv0(:), upd0(:), parent0(:), sstart0(:)
    integer(kint), allocatable :: merged_into(:), grp_rep(:), newid(:)

    n = lu%N

    nz_off = 0
    do row = 1, n
      do entry = sym_row_ptr(row), sym_row_ptr(row + 1) - 1
        if (sym_col_ind(entry) /= row) nz_off = nz_off + 1
      end do
    end do

    call newGraph(G, n, nz_off)
    G%xadj(0) = 0
    do row = 1, n
      cnt = 0
      do entry = sym_row_ptr(row), sym_row_ptr(row + 1) - 1
        if (sym_col_ind(entry) /= row) cnt = cnt + 1
      end do
      G%xadj(row) = G%xadj(row - 1) + cnt
    end do
    do row = 1, n
      pos = G%xadj(row - 1)
      do entry = sym_row_ptr(row), sym_row_ptr(row + 1) - 1
        col = sym_col_ind(entry)
        if (col /= row) then
          G%adjncy(pos) = col - 1
          pos = pos + 1
        end if
      end do
    end do

    opts_dummy = 0
    call monolis_pord_ordering(G, opts_dummy, .true., T)

    nf0 = T%nfronts

    call monolis_alloc_I_1d(lu%perm,  n)
    call monolis_alloc_I_1d(lu%iperm, n)
    !> permFromElimTree は old -> new を返す（1-based）。
    !> perm として保存し、iperm はその逆。
    call monolis_pord_perm_from_elimtree(T, n, lu%perm)
    do i = 1, n
      lu%iperm(lu%perm(i)) = i
    end do

    !> 元フロント（PORD 出力）の情報を作業配列へ取り出す
    call monolis_alloc_I_1d(vertex_front, n)
    call monolis_alloc_I_1d(piv0,    nf0)
    call monolis_alloc_I_1d(upd0,    nf0)
    call monolis_alloc_I_1d(parent0, nf0)
    call monolis_alloc_I_1d(sstart0, nf0)

    do row = 1, n
      vertex_front(row) = T%vtx2front(row - 1) + 1
    end do

    do front = 1, nf0
      piv0(front) = T%ncolfactor(front - 1)
      upd0(front) = T%ncolupdate(front - 1)
      if (T%parent(front - 1) == -1) then
        parent0(front) = 0
      else
        parent0(front) = T%parent(front - 1) + 1
      end if
    end do

    sstart0(:) = n + 1
    do row = 1, n
      new_col = lu%perm(row)
      front = vertex_front(row)
      sstart0(front) = min(sstart0(front), new_col)
    end do

    !> --- amalgamation（スーパーノード融合） ---
    !> 列方向に親と隣接する子フロント（sstart0(c)+piv0(c)==sstart0(p)）を、
    !> fundamental supernode（更新列が親フロント全体に一致、padding なし）か
    !> relaxed 条件（小フロントかつ padding 上限以下）を満たす場合に親へ融合する。
    !> 融合は列ブロックの連続性を保つため、必ず直接の親に対してのみ行う。
    call monolis_alloc_I_1d(merged_into, nf0)
    do front = 1, nf0
      p = parent0(front)
      if (p <= 0) cycle
      if (sstart0(front) + piv0(front) /= sstart0(p)) cycle
      size_f = piv0(front) + upd0(front)
      size_p = piv0(p) + upd0(p)
      do_merge = (upd0(front) == size_p)
      if (.not. do_merge) then
        do_merge = (piv0(front) <= amalg_nemin .and. size_p - size_f <= amalg_npad)
      end if
      if (do_merge) merged_into(front) = p
    end do

    !> 各元フロントの代表（融合チェーンの最上位フロント）を求める
    call monolis_alloc_I_1d(grp_rep, nf0)
    do front = 1, nf0
      r0 = front
      do while (merged_into(r0) /= 0)
        r0 = merged_into(r0)
      end do
      grp_rep(front) = r0
    end do

    !> 代表フロントへ新番号を付与
    call monolis_alloc_I_1d(newid, nf0)
    nfronts = 0
    do front = 1, nf0
      if (merged_into(front) == 0) then
        nfronts = nfronts + 1
        newid(front) = nfronts
      end if
    end do
    lu%nfronts = nfronts

    call monolis_alloc_I_1d(column_to_super,          n)
    call monolis_alloc_I_1d(lu%super_start,           nfronts)
    call monolis_alloc_I_1d(lu%front_parent,          nfronts)
    call monolis_alloc_I_1d(lu%front_first_child,     nfronts)
    call monolis_alloc_I_1d(lu%front_next_sibling,    nfronts)
    call monolis_alloc_I_1d(lu%front_postorder,       nfronts)
    call monolis_alloc_I_1d(lu%front_size,            nfronts)
    call monolis_alloc_I_1d(lu%front_pivot_size,      nfronts)
    call monolis_alloc_I_1d(lu%front_update_size,     nfronts)
    call monolis_alloc_I_1d(lu%front_ptr,             nfronts + 1)

    lu%super_start(:) = n + 1
    lu%front_pivot_size(:) = 0
    lu%front_first_child(:) = 0
    lu%front_next_sibling(:) = 0
    lu%front_postorder(:) = 0

    !> 融合後フロントのピボット数（融合元の総和）と super_start（最小）を集計
    do front = 1, nf0
      r0 = newid(grp_rep(front))
      lu%front_pivot_size(r0) = lu%front_pivot_size(r0) + piv0(front)
      lu%super_start(r0) = min(lu%super_start(r0), sstart0(front))
    end do

    !> 更新サイズと親は代表（最上位）フロントのものを採用
    do front = 1, nf0
      if (merged_into(front) /= 0) cycle
      r0 = newid(front)
      lu%front_update_size(r0) = upd0(front)
      if (parent0(front) == 0) then
        lu%front_parent(r0) = 0
      else
        lu%front_parent(r0) = newid(grp_rep(parent0(front)))
      end if
    end do

    lu%max_front_size = 0
    lu%front_ptr(1) = 1
    total_front_vars = 0
    do front = 1, nfronts
      lu%front_size(front) = lu%front_pivot_size(front) + lu%front_update_size(front)
      lu%max_front_size = max(lu%max_front_size, lu%front_size(front))
      total_front_vars = total_front_vars + lu%front_size(front)
      lu%front_ptr(front + 1) = lu%front_ptr(front) + lu%front_size(front)
    end do

    call monolis_alloc_I_1d(lu%front_ind, max(1, total_front_vars))

    column_to_super(:) = 0
    do row = 1, n
      new_col = lu%perm(row)
      front = newid(grp_rep(vertex_front(row)))
      column_to_super(new_col) = front
    end do

    !> firstchild / next_sibling
    do front = nfronts, 1, -1
      child = lu%front_parent(front)
      if (child > 0) then
        lu%front_next_sibling(front) = lu%front_first_child(child)
        lu%front_first_child(child) = front
      end if
    end do

    call monolis_dealloc_I_1d(piv0)
    call monolis_dealloc_I_1d(upd0)
    call monolis_dealloc_I_1d(parent0)
    call monolis_dealloc_I_1d(sstart0)
    call monolis_dealloc_I_1d(merged_into)
    call monolis_dealloc_I_1d(grp_rep)
    call monolis_dealloc_I_1d(newid)

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
        do entry = sym_row_ptr(row), sym_row_ptr(row + 1) - 1
          col = lu%perm(sym_col_ind(entry))
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
    call monolis_dealloc_I_1d(vertex_front)
    call freeElimTree(T)
    call freeGraph(G)
  end subroutine build_pord_front_structure

  !> 各非ゼロ要素を所属フロントに対応付け、orig_* を構築
  subroutine build_original_front_entries(lu, n, row_ptr, col_ind, val_pos, column_to_super)
    implicit none
    type(monolis_mat_lu), intent(inout) :: lu
    integer(kint), intent(in) :: n
    integer(kint), intent(in) :: row_ptr(:)
    integer(kint), intent(in) :: col_ind(:)
    integer(kint), intent(in) :: val_pos(:)
    integer(kint), intent(in) :: column_to_super(:)

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
        front = column_to_super(first_col)
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
        front = column_to_super(first_col)
        pos = next_pos(front)

        lu%orig_row_pos(pos) = new_row
        lu%orig_col_pos(pos) = new_col
        lu%orig_entry(pos)   = val_pos(entry)
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
