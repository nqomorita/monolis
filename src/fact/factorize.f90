module mod_monolis_fact_factorize
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc

  implicit none

contains

!==============================================================================
! Build elimination tree using Liu's algorithm (column-by-column)
!
! For each column k (in permuted order 1..N), find all edges to columns j < k
! and set parent(root_of_subtree(j)) = k using Union-Find with path compression.
!
! This requires building a CSR adjacency list from COO first.
!==============================================================================
subroutine build_elimination_tree(mat, lu)
  type(matrix_data), intent(in) :: mat
  type(monolis_mat_lu), intent(inout) :: lu

  integer(kint) :: n, nz
  integer(kint), allocatable :: ancestor(:)
  integer(kint), allocatable :: adj_ptr(:), adj_list(:)
  integer(kint), allocatable :: adj_count(:)
  integer(kint), allocatable :: flag(:)       ! MUMPS ANA_K-style O(1) dedup
  integer(kint) :: k, pi, pj, tmp, r, c, pos, i, nadj

  n = mat%n; nz = mat%nz
  allocate(lu%parent(n))

  associate(row_ptr => mat%row_ptr, col_ind => mat%col_ind, &
            invperm => mat%invperm, perm => mat%perm, parent => lu%parent)

  ! ---- Step 1: Build CSC adjacency for the permuted lower triangle ----
  ! Single-pass counting + FLAG-based O(1) duplicate elimination (MUMPS ANA_K style)
  allocate(adj_count(n), flag(n))
  adj_count = 0
  flag = 0

  ! Count phase with deduplication: for each column pj, count distinct pi < pj
  do r = 1, n
    pi = invperm(r)
    do k = row_ptr(r)+1, row_ptr(r+1)
      c = col_ind(k)
      if (r == c) cycle
      pj = invperm(c)
      ! Ensure pi_local < pj_local (lower triangle in permuted order)
      if (pi < pj) then
        if (flag(pi) /= pj) then   ! O(1) duplicate check
          flag(pi) = pj
          adj_count(pj) = adj_count(pj) + 1
        end if
      else
        if (flag(pj) /= pi) then
          flag(pj) = pi
          adj_count(pi) = adj_count(pi) + 1
        end if
      end if
    end do
  end do

  ! Build pointer array
  allocate(adj_ptr(n+1))
  adj_ptr(1) = 1
  do k = 2, n+1
    adj_ptr(k) = adj_ptr(k-1) + adj_count(k-1)
  end do
  nadj = adj_ptr(n+1) - 1

  ! Fill phase with deduplication
  allocate(adj_list(max(nadj, 1)))
  adj_count = 0
  flag = 0

  do r = 1, n
    pi = invperm(r)
    do k = row_ptr(r)+1, row_ptr(r+1)
      c = col_ind(k)
      if (r == c) cycle
      pj = invperm(c)
      if (pi < pj) then
        if (flag(pi) /= pj) then
          flag(pi) = pj
          adj_count(pj) = adj_count(pj) + 1
          adj_list(adj_ptr(pj) + adj_count(pj) - 1) = pi
        end if
      else
        if (flag(pj) /= pi) then
          flag(pj) = pi
          adj_count(pi) = adj_count(pi) + 1
          adj_list(adj_ptr(pi) + adj_count(pi) - 1) = pj
        end if
      end if
    end do
  end do

  deallocate(flag)

  ! ---- Step 2: Liu's elimination tree algorithm with path compression ----
  allocate(ancestor(n))
  parent = 0

  do k = 1, n
    ancestor(k) = k

    do pos = adj_ptr(k), adj_ptr(k+1) - 1
      pj = adj_list(pos)   ! pj < k

      ! Find root with path halving (faster than full path compression)
      r = pj
      do while (ancestor(r) /= r)
        ancestor(r) = ancestor(ancestor(r))  ! path halving
        r = ancestor(r)
      end do

      if (r /= k) then
        parent(r) = k
        ancestor(r) = k
      end if
    end do
  end do

  ! Ensure connected tree: orphans connect up
  do k = 1, n - 1
    if (parent(k) == 0) parent(k) = k + 1
  end do
  parent(n) = 0

  deallocate(ancestor, adj_ptr, adj_list, adj_count)
  end associate
end subroutine
!
! A supernode is a maximal set of contiguous columns in the factor
! with the same sparsity structure (fundamental supernode).
!
! Detection rule:  column j can be merged with column j-1 if
!   parent(j-1) == j   AND   nonzero-count(j) == nonzero-count(j-1) - 1
!
! We also apply relaxed merging (nemin) to allow a small amount of fill.
!==============================================================================
subroutine identify_supernodes(mat, lu)
  type(matrix_data), intent(in) :: mat
  type(monolis_mat_lu), intent(inout) :: lu

  integer(kint) :: n, nz, nsuper
  integer(kint), allocatable :: sbelong(:), sstart(:), ssize(:), sparent(:), sfsize(:)
  integer(kint), allocatable :: col_count(:)  ! nonzero count in each column of L (permuted)
  integer(kint), allocatable :: children_count(:)
  integer(kint), allocatable :: flag(:)       ! FLAG array for O(1) dedup (MUMPS ANA_K style)
  integer(kint) :: i, k, r, c, pi, pj
  integer(kint) :: nemin
  logical :: can_merge

  n = mat%n; nz = mat%nz
  nemin = 16  ! relaxation parameter

  allocate(col_count(n), children_count(n), sbelong(n), flag(n))
  col_count = 1   ! diagonal always present
  children_count = 0
  flag = 0

  associate(row_ptr => mat%row_ptr, col_ind => mat%col_ind, &
            perm => mat%perm, invperm => mat%invperm, parent => lu%parent)

  ! ---- Step 1: Estimate column counts of L ----
  ! FLAG-based O(1) duplicate elimination (MUMPS ANA_K style)
  ! For each off-diagonal entry, the column with the smaller permuted index gets +1.
  do r = 1, n
    pi = invperm(r)
    do k = row_ptr(r)+1, row_ptr(r+1)
      c = col_ind(k)
      if (r == c) cycle
      pj = invperm(c)
      if (pi < pj) then
        if (flag(pi) /= pj) then  ! O(1) dedup: avoid counting same (pi,pj) twice
          flag(pi) = pj
          col_count(pi) = col_count(pi) + 1
        end if
      else
        if (flag(pj) /= pi) then
          flag(pj) = pi
          col_count(pj) = col_count(pj) + 1
        end if
      end if
    end do
  end do

  ! Propagate fill through the tree (MUMPS ANA_LNEW style):
  ! col_count(parent(i)) gets contribution from child i
  do i = 1, n
    if (parent(i) > 0 .and. parent(i) <= n) then
      if (col_count(i) > 1) then
        col_count(parent(i)) = max(col_count(parent(i)), col_count(i) - 1)
      end if
    end if
  end do

  ! Count children
  do i = 1, n
    if (parent(i) > 0 .and. parent(i) <= n) then
      children_count(parent(i)) = children_count(parent(i)) + 1
    end if
  end do

  ! ---- Step 2: Fundamental supernode detection ----
  sbelong(1) = 1
  nsuper = 1
  do i = 2, n
    can_merge = .true.
    if (parent(i-1) /= i) can_merge = .false.
    if (children_count(i) > 1) can_merge = .false.
    if (can_merge .and. col_count(i) > 0 .and. col_count(i-1) > 0) then
      if (col_count(i) /= col_count(i-1) - 1) then
        if (col_count(i-1) > nemin .and. col_count(i) > nemin) then
          can_merge = .false.
        end if
      end if
    end if

    if (can_merge) then
      sbelong(i) = nsuper
    else
      nsuper = nsuper + 1
      sbelong(i) = nsuper
    end if
  end do

  ! ---- Step 3: Build supernode arrays ----
  allocate(sstart(nsuper), ssize(nsuper), sparent(nsuper), sfsize(nsuper))
  sstart = 0; ssize = 0; sparent = 0; sfsize = 0

  do i = 1, n
    k = sbelong(i)
    ssize(k) = ssize(k) + 1
    if (sstart(k) == 0) sstart(k) = i
  end do

  do i = 1, nsuper
    k = sstart(i) + ssize(i) - 1
    if (parent(k) > 0 .and. parent(k) <= n) then
      sparent(i) = sbelong(parent(k))
    else
      sparent(i) = 0
    end if
  end do

  ! ---- Step 4: Relaxed merging of small supernodes (nemin) ----
  call relax_supernodes(n, nsuper, sbelong, sstart, ssize, sparent, sfsize, &
                        parent, col_count, nemin)

  ! ---- Step 5: Compute frontal sizes (O(nz) CSC-based approach) ----
  call compute_frontal_sizes(n, nz, row_ptr, col_ind, invperm, &
                             nsuper, sstart, ssize, sparent, sbelong, sfsize)

  end associate

  lu%nsuper = nsuper
  lu%snode_belong = sbelong
  lu%snode_start = sstart
  lu%snode_size = ssize
  lu%snode_parent = sparent
  lu%snode_fsize = sfsize

  deallocate(col_count, children_count, flag)
end subroutine

!==============================================================================
! Relaxed supernode merging: merge parent-child pairs when both are small
! Optimized: O(n) sbelong update using sstart/ssize ranges instead of O(n*nsuper)
!==============================================================================
subroutine relax_supernodes(n, nsuper, sbelong, sstart, ssize, sparent, sfsize, &
                            parent, col_count, nemin)
  integer(kint), intent(in)    :: n, nemin
  integer(kint), intent(inout) :: nsuper
  integer(kint), intent(inout) :: sbelong(n)
  integer(kint), intent(inout), allocatable :: sstart(:), ssize(:), sparent(:), sfsize(:)
  integer(kint), intent(in)    :: parent(n), col_count(n)

  logical, allocatable :: merged(:)
  integer(kint) :: i, s, sp, new_nsuper, j
  integer(kint), allocatable :: new_sstart(:), new_ssize(:), new_sparent(:), new_sfsize(:)
  integer(kint), allocatable :: map(:)

  allocate(merged(nsuper))
  merged = .false.

  ! Mark child supernodes for merging into parent when both are small
  do s = 1, nsuper
    sp = sparent(s)
    if (sp > 0 .and. sp <= nsuper .and. sp /= s) then
      if (ssize(s) <= nemin .and. ssize(sp) <= nemin) then
        if (sstart(sp) == sstart(s) + ssize(s)) then
          merged(s) = .true.
        end if
      end if
    end if
  end do

  ! Apply merges using direct range update (O(total_merged_size) instead of O(n*nmerged))
  do s = 1, nsuper
    if (merged(s)) then
      sp = sparent(s)
      sstart(sp) = min(sstart(sp), sstart(s))
      ssize(sp) = ssize(sp) + ssize(s)
      ! Direct range update: only touch variables belonging to supernode s
      do j = sstart(s), sstart(s) + ssize(s) - 1
        if (j >= 1 .and. j <= n) sbelong(j) = sp
      end do
    end if
  end do

  ! Renumber supernodes contiguously
  allocate(map(nsuper))
  map = 0
  new_nsuper = 0
  do s = 1, nsuper
    if (.not. merged(s)) then
      new_nsuper = new_nsuper + 1
      map(s) = new_nsuper
    end if
  end do

  allocate(new_sstart(new_nsuper), new_ssize(new_nsuper))
  allocate(new_sparent(new_nsuper), new_sfsize(new_nsuper))

  do s = 1, nsuper
    if (.not. merged(s)) then
      i = map(s)
      new_sstart(i) = sstart(s)
      new_ssize(i) = ssize(s)
      new_sfsize(i) = sfsize(s)
      if (sparent(s) > 0 .and. sparent(s) <= nsuper) then
        sp = sparent(s)
        do while (merged(sp) .and. sp > 0 .and. sp <= nsuper)
          sp = sparent(sp)
        end do
        if (sp > 0 .and. sp <= nsuper) then
          new_sparent(i) = map(sp)
        else
          new_sparent(i) = 0
        end if
      else
        new_sparent(i) = 0
      end if
    end if
  end do

  ! Update sbelong with map
  do i = 1, n
    sbelong(i) = map(sbelong(i))
  end do

  deallocate(sstart, ssize, sparent, sfsize)
  nsuper = new_nsuper
  allocate(sstart(nsuper), ssize(nsuper), sparent(nsuper), sfsize(nsuper))
  sstart  = new_sstart
  ssize   = new_ssize
  sparent = new_sparent
  sfsize  = new_sfsize

  deallocate(merged, map, new_sstart, new_ssize, new_sparent, new_sfsize)
end subroutine

!==============================================================================
! Compute frontal sizes for each supernode
!
! Optimized O(nz) approach: build permuted CSC, then scan only columns
! belonging to each supernode. Uses marker array with supernode-based
! timestamp to avoid O(n) resets per supernode.
!==============================================================================
subroutine compute_frontal_sizes(n, nz, row_ptr, col_ind, invperm, &
                                 nsuper, sstart, ssize, sparent, sbelong, sfsize)
  integer(kint), intent(in) :: n, nz
  integer(kint), intent(in) :: row_ptr(n+1), col_ind(nz)
  integer(kint), intent(in) :: invperm(n)
  integer(kint), intent(in) :: nsuper
  integer(kint), intent(in) :: sstart(nsuper), ssize(nsuper)
  integer(kint), intent(in) :: sparent(nsuper), sbelong(n)
  integer(kint), intent(inout) :: sfsize(nsuper)

  ! Permuted CSC representation
  integer(kint), allocatable :: csc_ptr(:), csc_row(:)
  integer(kint), allocatable :: csc_count(:)
  ! Marker with supernode-id timestamp (avoids reset)
  integer(kint), allocatable :: marker(:)
  integer(kint) :: s, k, pi, pj, r, c, cnt, nadj
  integer(kint) :: first_col, last_col, col, pos
  integer(kint) :: child_cnt, child_fsize, child_npiv

  ! ---- Build permuted CSC: for each permuted column pj, list of rows pi ----
  allocate(csc_count(n))
  csc_count = 0

  ! Count entries per permuted column (lower triangle: pi > pj stored in col pj)
  do r = 1, n
    pi = invperm(r)
    do k = row_ptr(r)+1, row_ptr(r+1)
      c = col_ind(k)
      if (r == c) cycle
      pj = invperm(c)
      if (pi > pj) then
        csc_count(pj) = csc_count(pj) + 1
      else
        csc_count(pi) = csc_count(pi) + 1
      end if
    end do
  end do

  allocate(csc_ptr(n+1))
  csc_ptr(1) = 1
  do k = 2, n+1
    csc_ptr(k) = csc_ptr(k-1) + csc_count(k-1)
  end do
  nadj = csc_ptr(n+1) - 1

  allocate(csc_row(max(nadj, 1)))
  csc_count = 0

  do r = 1, n
    pi = invperm(r)
    do k = row_ptr(r)+1, row_ptr(r+1)
      c = col_ind(k)
      if (r == c) cycle
      pj = invperm(c)
      if (pi > pj) then
        csc_count(pj) = csc_count(pj) + 1
        csc_row(csc_ptr(pj) + csc_count(pj) - 1) = pi
      else
        csc_count(pi) = csc_count(pi) + 1
        csc_row(csc_ptr(pi) + csc_count(pi) - 1) = pj
      end if
    end do
  end do

  deallocate(csc_count)

  ! ---- Compute frontal sizes using CSC: O(nz) total ----
  ! marker(pi) = s means pi is already counted for supernode s (avoids O(n) reset)
  allocate(marker(n))
  marker = 0

  ! Process supernodes bottom-up (leaves first) to propagate CB indices
  ! Simple approach: process in order 1..nsuper (postorder approximation for chains)
  do s = 1, nsuper
    first_col = sstart(s)
    last_col  = sstart(s) + ssize(s) - 1
    cnt = 0

    ! Scan CSC entries for columns in [first_col..last_col]
    do col = first_col, last_col
      do pos = csc_ptr(col), csc_ptr(col+1) - 1
        pi = csc_row(pos)
        if (pi > last_col .and. marker(pi) /= s) then
          marker(pi) = s
          cnt = cnt + 1
        end if
      end do
    end do

    sfsize(s) = ssize(s) + cnt
    if (sfsize(s) < ssize(s)) sfsize(s) = ssize(s)
  end do

  deallocate(marker, csc_ptr, csc_row)
end subroutine

!==============================================================================
! Build child-sibling representation for the supernode tree
! Optimized: O(nsuper) head-prepend instead of O(nsuper^2) tail-append
!==============================================================================
subroutine build_frontal_tree(mat, lu)
  type(matrix_data), intent(in) :: mat
  type(monolis_mat_lu), intent(inout) :: lu

  integer(kint) :: nsuper, s, p

  nsuper = lu%nsuper
  allocate(lu%sfils(nsuper), lu%sfrere(nsuper))
  lu%sfils  = 0
  lu%sfrere = 0

  associate(sparent => lu%snode_parent, sfils => lu%sfils, sfrere => lu%sfrere)

  ! Head-prepend: O(1) per supernode, O(nsuper) total
  ! Each child is inserted at the head of its parent's child list
  do s = nsuper, 1, -1
    p = sparent(s)
    if (p > 0 .and. p <= nsuper) then
      sfrere(s) = sfils(p)   ! current child becomes sibling of new head
      sfils(p) = s            ! new child becomes head
    end if
  end do

  end associate
end subroutine

!==============================================================================
! Multifrontal lu factorization
!
! Two-pass approach for minimal overhead:
!   Pass 1 (Symbolic): compute index sets, pre-allocate all frontal matrices
!   Pass 2 (Numeric):  assemble + extend-add + factor (ZERO allocations)
!
! Optimized with:
!   - Pre-built permuted CSC with 3 separate arrays (cache-friendly)
!   - All-at-once memory allocation (batch, not per-supernode)
!   - Pre-allocated workspaces for index sets and CB maps
!   - Inline extend-add (no DAXPY overhead for small ndof)
!   - Adaptive LU kernel with small-front inlining
!==============================================================================
subroutine multifrontal_factorize(mat, lu)
  type(matrix_data), intent(in) :: mat
  type(monolis_mat_lu), intent(inout) :: lu

  integer(kint) :: n, nz, nsuper, ndof, info
  integer(kint), allocatable :: postorder(:)
  integer(kint) :: s, ip, nfront, npiv, k, i, j
  integer(kint) :: first_col, last_col
  integer(kint) :: pi, pj, cnt, col, pos
  integer(kint), allocatable :: idx_work(:)   ! reusable workspace
  integer(kint) :: child
  integer(kint), allocatable :: imap(:)
  integer(kint), allocatable :: csc_ptr(:), csc_pi(:), csc_pj(:), csc_origk(:)
  integer(kint), allocatable :: cb_work(:)    ! reusable CB map workspace
  integer(kint) :: max_cb

  n = mat%n; nz = mat%nz; nsuper = lu%nsuper; ndof = mat%ndof
  info = 0

  allocate(lu%factors(nsuper))

  associate(row_ptr => mat%row_ptr, col_ind => mat%col_ind, a_elt => mat%a_elt, &
            perm => mat%perm, invperm => mat%invperm, &
            sstart => lu%snode_start, ssize => lu%snode_size, &
            sparent => lu%snode_parent, sfsize => lu%snode_fsize, &
            sfils => lu%sfils, sfrere => lu%sfrere, factors => lu%factors)

  ! ================================================================
  ! Phase 1: Build permuted CSC (3 separate arrays, cache-friendly)
  ! ================================================================
  call build_permuted_csc(n, nz, row_ptr, col_ind, invperm, &
                          csc_ptr, csc_pi, csc_pj, csc_origk)
  call build_postorder(nsuper, sfils, sfrere, sparent, postorder)

  ! ================================================================
  ! Phase 2: Symbolic factorization — compute exact index sets
  ! Pre-allocate idx_work ONCE (eliminates N per-supernode allocations)
  ! ================================================================
  allocate(imap(n), idx_work(n))
  imap = 0

  do ip = 1, nsuper
    s = postorder(ip)
    npiv = ssize(s)
    first_col = sstart(s)
    last_col  = first_col + npiv - 1

    ! Pivot columns — mark via imap as temporary flag
    do i = first_col, last_col
      imap(i) = -1
    end do
    cnt = npiv
    do i = 1, npiv
      idx_work(i) = first_col + i - 1
    end do

    ! Row indices from CSC columns [first_col..last_col]
    do col = first_col, last_col
      do pos = csc_ptr(col), csc_ptr(col+1) - 1
        pi = csc_pi(pos)
        pj = csc_pj(pos)
        if (pi > last_col .and. imap(pi) == 0) then
          cnt = cnt + 1; idx_work(cnt) = pi; imap(pi) = -1
        end if
        if (pj > last_col .and. imap(pj) == 0) then
          cnt = cnt + 1; idx_work(cnt) = pj; imap(pj) = -1
        end if
      end do
    end do

    ! Children's CB indices
    child = sfils(s)
    do while (child /= 0)
      do i = factors(child)%npiv + 1, factors(child)%nfront
        j = factors(child)%indices(i)
        if (j >= first_col .and. imap(j) == 0) then
          cnt = cnt + 1; idx_work(cnt) = j; imap(j) = -1
        end if
      end do
      child = sfrere(child)
    end do

    nfront = cnt
    if (nfront - npiv > 1) call quicksort(idx_work(npiv+1:nfront), nfront - npiv)

    ! Store symbolic structure
    factors(s)%nfront = nfront
    factors(s)%npiv   = npiv
    allocate(factors(s)%indices(nfront))
    factors(s)%indices(1:nfront) = idx_work(1:nfront)

    ! Clear imap
    do i = 1, nfront
      imap(idx_work(i)) = 0
    end do
  end do

  ! ================================================================
  ! Phase 3: Batch-allocate all frontal matrices + CB workspace
  ! ================================================================
  max_cb = 0
  do s = 1, nsuper
    nfront = factors(s)%nfront
    npiv   = factors(s)%npiv
    allocate(factors(s)%front(nfront*ndof, nfront*ndof))
    allocate(factors(s)%pivorder(npiv))
    do i = 1, npiv
      factors(s)%pivorder(i) = i
    end do
    if (nfront - npiv > max_cb) max_cb = nfront - npiv
  end do
  allocate(cb_work(max(max_cb, 1)))

  ! ================================================================
  ! Phase 4: Numeric factorization (ZERO allocations in this loop)
  ! ================================================================
  do ip = 1, nsuper
    s = postorder(ip)
    nfront = factors(s)%nfront
    npiv   = factors(s)%npiv
    first_col = sstart(s)
    last_col  = first_col + npiv - 1

    ! Build imap
    do i = 1, nfront
      imap(factors(s)%indices(i)) = i
    end do

    ! Zero front
    factors(s)%front = 0.0d0

    ! Assemble original entries from CSC
    call assemble_original_csc(n, ndof, csc_ptr, csc_pi, csc_pj, csc_origk, a_elt, &
                               nfront, factors(s)%front, first_col, last_col, imap)

    ! Extend-add from children (allocation-free using cb_work)
    child = sfils(s)
    do while (child /= 0)
      call extend_add_fast(factors(child), factors(s), imap, ndof, cb_work)
      child = sfrere(child)
    end do

    ! LU factorization with adaptive kernel
    call frontal_lu_factor_opt(factors(s)%front, nfront*ndof, npiv*ndof, info)
    if (info /= 0) info = 0

    ! Clear imap
    do i = 1, nfront
      imap(factors(s)%indices(i)) = 0
    end do
  end do

  deallocate(imap, idx_work, cb_work, postorder)
  deallocate(csc_ptr, csc_pi, csc_pj, csc_origk)
  end associate
end subroutine

!==============================================================================
! Build permuted CSC: 3 separate arrays (csc_pi, csc_pj, csc_origk)
! For each CSR entry, store in column min(pi,pj) with original direction.
! Cache-friendly: no interleaving, stride-1 access per array.
!==============================================================================
subroutine build_permuted_csc(n, nz, row_ptr, col_ind, invperm, &
                              csc_ptr, csc_pi, csc_pj, csc_origk)
  integer(kint), intent(in) :: n, nz
  integer(kint), intent(in) :: row_ptr(n+1), col_ind(nz)
  integer(kint), intent(in) :: invperm(n)
  integer(kint), allocatable, intent(out) :: csc_ptr(:), csc_pi(:), csc_pj(:), csc_origk(:)

  integer(kint), allocatable :: csc_count(:)
  integer(kint) :: r, k, c, pi, pj, nadj, mn

  allocate(csc_count(n))
  csc_count = 0

  do r = 1, n
    pi = invperm(r)
    do k = row_ptr(r)+1, row_ptr(r+1)
      c = col_ind(k)
      pj = invperm(c)
      mn = min(pi, pj)
      csc_count(mn) = csc_count(mn) + 1
    end do
  end do

  allocate(csc_ptr(n+1))
  csc_ptr(1) = 1
  do k = 2, n+1
    csc_ptr(k) = csc_ptr(k-1) + csc_count(k-1)
  end do
  nadj = csc_ptr(n+1) - 1

  allocate(csc_pi(nadj), csc_pj(nadj), csc_origk(nadj))
  csc_count = 0

  do r = 1, n
    pi = invperm(r)
    do k = row_ptr(r)+1, row_ptr(r+1)
      c = col_ind(k)
      pj = invperm(c)
      mn = min(pi, pj)
      csc_count(mn) = csc_count(mn) + 1
      csc_pi(csc_ptr(mn) + csc_count(mn) - 1) = pi
      csc_pj(csc_ptr(mn) + csc_count(mn) - 1) = pj
      csc_origk(csc_ptr(mn) + csc_count(mn) - 1) = k
    end do
  end do

  deallocate(csc_count)
end subroutine

!==============================================================================
! Assemble original entries from 3-array CSC (O(nnz_local) per supernode)
!==============================================================================
subroutine assemble_original_csc(n, ndof, csc_ptr, csc_pi, csc_pj, csc_origk, a_elt, &
                                  nfront, front, first_col, last_col, imap)
  integer(kint), intent(in) :: n, ndof
  integer(kint), intent(in) :: csc_ptr(n+1)
  integer(kint), intent(in) :: csc_pi(*), csc_pj(*), csc_origk(*)
  real(kdouble), intent(in) :: a_elt(*)
  integer(kint), intent(in) :: nfront
  real(kdouble), intent(inout) :: front(nfront*ndof, nfront*ndof)
  integer(kint), intent(in) :: first_col, last_col
  integer(kint), intent(in) :: imap(n)

  integer(kint) :: col, pos, pi, pj, ii, jj, id, jd, base, origk

  do col = first_col, last_col
    do pos = csc_ptr(col), csc_ptr(col+1) - 1
      pi = csc_pi(pos)
      pj = csc_pj(pos)
      ii = imap(pi)
      jj = imap(pj)
      if (ii > 0 .and. jj > 0) then
        origk = csc_origk(pos)
        base = (origk-1)*ndof*ndof
        do jd = 1, ndof
          do id = 1, ndof
            front((ii-1)*ndof+id, (jj-1)*ndof+jd) = &
              front((ii-1)*ndof+id, (jj-1)*ndof+jd) + a_elt(base + (id-1)*ndof + jd)
          end do
        end do
      end if
    end do
  end do
end subroutine

!==============================================================================
! Quicksort for integer arrays (replaces O(n^2) insertion sort)
!==============================================================================
subroutine quicksort(arr, n)
  integer(kint), intent(inout) :: arr(n)
  integer(kint), intent(in) :: n

  integer(kint) :: stack(2, 64)  ! log2(max_n) levels suffice
  integer(kint) :: sp, lo, hi, i, j, pivot, tmp, mid

  if (n <= 1) return

  ! Fall back to insertion sort for small arrays
  if (n <= 32) then
    call isort(arr, n)
    return
  end if

  sp = 1
  stack(1, 1) = 1
  stack(2, 1) = n

  do while (sp > 0)
    lo = stack(1, sp)
    hi = stack(2, sp)
    sp = sp - 1

    if (hi - lo < 16) then
      ! Insertion sort for small partitions
      do i = lo + 1, hi
        tmp = arr(i)
        j = i - 1
        do while (j >= lo)
          if (arr(j) <= tmp) exit
          arr(j+1) = arr(j)
          j = j - 1
        end do
        arr(j+1) = tmp
      end do
      cycle
    end if

    ! Median-of-three pivot selection
    mid = (lo + hi) / 2
    if (arr(lo) > arr(mid)) then; tmp = arr(lo); arr(lo) = arr(mid); arr(mid) = tmp; end if
    if (arr(lo) > arr(hi))  then; tmp = arr(lo); arr(lo) = arr(hi);  arr(hi)  = tmp; end if
    if (arr(mid) > arr(hi)) then; tmp = arr(mid); arr(mid) = arr(hi); arr(hi) = tmp; end if
    pivot = arr(mid)
    arr(mid) = arr(hi - 1)
    arr(hi - 1) = pivot

    i = lo
    j = hi - 1
    do
      do
        i = i + 1
        if (arr(i) >= pivot) exit
      end do
      do
        j = j - 1
        if (arr(j) <= pivot) exit
      end do
      if (i >= j) exit
      tmp = arr(i); arr(i) = arr(j); arr(j) = tmp
    end do
    arr(hi - 1) = arr(i)
    arr(i) = pivot

    ! Push larger partition first (ensures O(log n) stack depth)
    if (i - lo > hi - i) then
      if (lo < i - 1) then; sp = sp + 1; stack(1, sp) = lo; stack(2, sp) = i - 1; end if
      if (i + 1 < hi) then; sp = sp + 1; stack(1, sp) = i + 1; stack(2, sp) = hi; end if
    else
      if (i + 1 < hi) then; sp = sp + 1; stack(1, sp) = i + 1; stack(2, sp) = hi; end if
      if (lo < i - 1) then; sp = sp + 1; stack(1, sp) = lo; stack(2, sp) = i - 1; end if
    end if
  end do
end subroutine

!==============================================================================
! Build postorder traversal of the supernode tree
!==============================================================================
subroutine build_postorder(nsuper, sfils, sfrere, sparent, postorder)
  integer(kint), intent(in)  :: nsuper
  integer(kint), intent(in)  :: sfils(nsuper), sfrere(nsuper), sparent(nsuper)
  integer(kint), allocatable, intent(out) :: postorder(:)

  integer(kint), allocatable :: stack(:)
  integer(kint) :: sp, s, cnt, child
  logical, allocatable :: visited(:)

  allocate(postorder(nsuper), stack(nsuper), visited(nsuper))
  visited = .false.
  cnt = 0
  sp = 0

  ! Push all roots
  do s = 1, nsuper
    if (sparent(s) == 0) then
      sp = sp + 1
      stack(sp) = s
    end if
  end do

  do while (sp > 0)
    s = stack(sp)

    ! Check if all children are visited
    child = sfils(s)
    if (child /= 0) then
      if (.not. visited(child)) then
        ! Push unvisited children
        do while (child /= 0)
          if (.not. visited(child)) then
            sp = sp + 1
            stack(sp) = child
          end if
          child = sfrere(child)
        end do
        cycle
      end if
    end if

    ! All children visited (or leaf): pop and record
    sp = sp - 1
    if (.not. visited(s)) then
      visited(s) = .true.
      cnt = cnt + 1
      postorder(cnt) = s
    end if
  end do

  deallocate(stack, visited)
end subroutine

!==============================================================================
! Simple insertion sort for small integer(kint) arrays (used by quicksort)
!==============================================================================
subroutine isort(arr, n)
  integer(kint), intent(inout) :: arr(n)
  integer(kint), intent(in) :: n
  integer(kint) :: i, j, key

  do i = 2, n
    key = arr(i)
    j = i - 1
    do while (j >= 1)
      if (arr(j) <= key) exit
      arr(j+1) = arr(j)
      j = j - 1
    end do
    arr(j+1) = key
  end do
end subroutine

!==============================================================================
! Extend-add: allocation-free, using pre-allocated cb_work workspace.
! Inline scatter-add (no DAXPY call overhead for small ndof).
!==============================================================================
subroutine extend_add_fast(child_fac, parent_fac, imap, ndof, cb_work)
  type(monolis_mat_frontal), intent(in)    :: child_fac
  type(monolis_mat_frontal), intent(inout) :: parent_fac
  integer(kint), intent(in) :: imap(*)
  integer(kint), intent(in) :: ndof
  integer(kint), intent(inout) :: cb_work(*)

  integer(kint) :: ic, jc, ip, jp, id, jd
  integer(kint) :: npiv_c, nfront_c, ncb

  npiv_c  = child_fac%npiv
  nfront_c = child_fac%nfront
  ncb = nfront_c - npiv_c
  if (ncb <= 0) return

  ! Build CB map using pre-allocated workspace (ZERO allocation)
  do ic = 1, ncb
    cb_work(ic) = imap(child_fac%indices(npiv_c + ic))
  end do

  ! Inline scatter-add (avoids DAXPY function call overhead for small ndof)
  do jc = 1, ncb
    jp = cb_work(jc)
    if (jp == 0) cycle
    do ic = 1, ncb
      ip = cb_work(ic)
      if (ip == 0) cycle
      do jd = 1, ndof
        do id = 1, ndof
          parent_fac%front((ip-1)*ndof+id, (jp-1)*ndof+jd) = &
            parent_fac%front((ip-1)*ndof+id, (jp-1)*ndof+jd) + &
            child_fac%front((npiv_c+ic-1)*ndof+id, (npiv_c+jc-1)*ndof+jd)
        end do
      end do
    end do
  end do
end subroutine

!==============================================================================
! Adaptive LU factorization (MUMPS FAC_H / FAC_P style)
!
! Small fronts (<=24): inline LU without BLAS (avoids function call overhead)
! Medium fronts: single panel (DSCAL + DGER only)
! Large fronts: panel-based DGER + DTRSM + DGEMM (Level-3 BLAS)
!==============================================================================
subroutine frontal_lu_factor_opt(front, nfront, npiv, info)
  real(kdouble), intent(inout) :: front(nfront, nfront)
  integer(kint), intent(in)    :: nfront, npiv
  integer(kint), intent(inout) :: info

  integer(kint) :: k, i, j, kb, ke, nb, nel, nupd, lkjib
  real(kdouble) :: pivot, inv_pivot
  real(kdouble), parameter :: eps = 1.0d-14
  real(kdouble), parameter :: one = 1.0d0, alpha = -1.0d0

  info = 0
  if (npiv <= 0) return

  ! ================================================================
  ! Small-front fast path: inline LU without BLAS calls
  ! Avoids DSCAL/DGER function call overhead for tiny matrices
  ! ================================================================
  if (nfront <= 24) then
    do k = 1, npiv
      pivot = front(k, k)
      if (abs(pivot) < eps) then
        front(k, k) = sign(eps, pivot + eps)
        pivot = front(k, k)
        info = k
      end if
      inv_pivot = one / pivot
      do i = k + 1, nfront
        front(i, k) = front(i, k) * inv_pivot
      end do
      do j = k + 1, nfront
        do i = k + 1, nfront
          front(i, j) = front(i, j) - front(i, k) * front(k, j)
        end do
      end do
    end do
    return
  end if

  ! ================================================================
  ! Adaptive panel size (MUMPS KEEP(4)/KEEP(5)/KEEP(6) style)
  ! ================================================================
  if (nfront < 128) then
    lkjib = nfront
  else if (nfront < 512) then
    lkjib = 48
  else
    lkjib = 64
  end if

  ! ================================================================
  ! Panel-based blocked factorization with Level-3 BLAS
  ! ================================================================
  kb = 1
  do while (kb <= npiv)
    ke = min(kb + lkjib - 1, npiv)
    nb = ke - kb + 1

    do k = kb, ke
      pivot = front(k, k)
      if (abs(pivot) < eps) then
        front(k, k) = sign(eps, pivot + eps)
        pivot = front(k, k)
        info = k
      end if

      nel = nfront - k
      if (nel > 0) then
        call dscal(nel, one/pivot, front(k+1, k), 1)
      end if

      nupd = ke - k
      if (nupd > 0 .and. nel > 0) then
        call dger(nel, nupd, alpha, front(k+1, k), 1, &
                  front(k, k+1), nfront, front(k+1, k+1), nfront)
      end if
    end do

    nupd = nfront - ke
    if (nupd > 0 .and. nb > 0) then
      call dtrsm('L', 'L', 'N', 'U', nb, nupd, one, &
                 front(kb, kb), nfront, front(kb, ke+1), nfront)

      call dgemm('N', 'N', nupd, nupd, nb, alpha, &
                 front(ke+1, kb), nfront, &
                 front(kb, ke+1), nfront, &
                 one, front(ke+1, ke+1), nfront)
    end if

    kb = ke + 1
  end do
end subroutine

!==============================================================================
! Generate right-hand side: b = A * [1, 1, ..., 1]
!==============================================================================
subroutine generate_rhs(mat, lu, rhs)
  type(matrix_data), intent(in) :: mat
  type(monolis_mat_lu), intent(in) :: lu
  real(kdouble), intent(out) :: rhs(mat%ndof*mat%n)

  integer(kint) :: r, k, n, ndof, id, jd, base

  n = mat%n; ndof = mat%ndof
  rhs = 0.0d0
  do r = 1, n
    do k = mat%row_ptr(r)+1, mat%row_ptr(r+1)
      base = (k-1)*ndof*ndof
      do id = 1, ndof
        do jd = 1, ndof
          rhs((r-1)*ndof+id) = rhs((r-1)*ndof+id) + mat%a_elt(base + (id-1)*ndof + jd)
        end do
      end do
    end do
  end do
end subroutine

!==============================================================================
! Multifrontal solve: L U x = P b
!
! Phase 1: Forward substitution  (L y = P b)  — bottom-up through tree
! Phase 2: Backward substitution (U x = y)    — top-down through tree
!==============================================================================
subroutine multifrontal_solve(mat, lu, rhs, x)
  type(matrix_data), intent(in) :: mat
  type(monolis_mat_lu), intent(in) :: lu
  real(kdouble), intent(in)  :: rhs(mat%ndof*mat%n)
  real(kdouble), intent(out) :: x(mat%ndof*mat%n)

  integer(kint) :: n, ndof, ntot, nsuper
  integer(kint), allocatable :: postorder(:)
  real(kdouble), allocatable :: work(:)
  integer(kint) :: ip, s, k, i, j, nfront, npiv
  integer(kint) :: kblk, kdof, iblk, idof, gi, gj
  integer(kint) :: nf_dof, np_dof
  real(kdouble) :: sum_val

  n = mat%n; ndof = mat%ndof; nsuper = lu%nsuper
  ntot = n * ndof

  associate(sstart => lu%snode_start, ssize => lu%snode_size, &
            sfsize => lu%snode_fsize, sparent => lu%snode_parent, &
            sfils => lu%sfils, sfrere => lu%sfrere, &
            perm => mat%perm, invperm => mat%invperm, &
            factors => lu%factors)

  allocate(work(ntot), postorder(nsuper))

  ! Apply block permutation to rhs
  do i = 1, n
    do k = 1, ndof
      work((invperm(i)-1)*ndof + k) = rhs((i-1)*ndof + k)
    end do
  end do

  call build_postorder_array(nsuper, sfils, sfrere, sparent, postorder)

  ! ==== Forward substitution (L y = Pb) in postorder ====
  do ip = 1, nsuper
    s = postorder(ip)
    nfront = factors(s)%nfront
    npiv   = factors(s)%npiv
    nf_dof = nfront * ndof
    np_dof = npiv * ndof

    do k = 1, np_dof
      kblk = (k-1)/ndof + 1
      kdof = mod(k-1, ndof) + 1
      gi = (factors(s)%indices(kblk) - 1)*ndof + kdof
      do i = k + 1, nf_dof
        iblk = (i-1)/ndof + 1
        idof = mod(i-1, ndof) + 1
        gj = (factors(s)%indices(iblk) - 1)*ndof + idof
        work(gj) = work(gj) - factors(s)%front(i, k) * work(gi)
      end do
    end do
  end do

  ! ==== Backward substitution (U x = y) in reverse postorder ====
  do ip = nsuper, 1, -1
    s = postorder(ip)
    nfront = factors(s)%nfront
    npiv   = factors(s)%npiv
    nf_dof = nfront * ndof
    np_dof = npiv * ndof

    do k = np_dof, 1, -1
      kblk = (k-1)/ndof + 1
      kdof = mod(k-1, ndof) + 1
      gi = (factors(s)%indices(kblk) - 1)*ndof + kdof
      sum_val = 0.0d0
      do j = k + 1, nf_dof
        iblk = (j-1)/ndof + 1
        idof = mod(j-1, ndof) + 1
        gj = (factors(s)%indices(iblk) - 1)*ndof + idof
        sum_val = sum_val + factors(s)%front(k, j) * work(gj)
      end do
      work(gi) = (work(gi) - sum_val) / factors(s)%front(k, k)
    end do
  end do

  ! Apply inverse permutation to get solution
  do i = 1, n
    do k = 1, ndof
      x((i-1)*ndof + k) = work((invperm(i)-1)*ndof + k)
    end do
  end do

  deallocate(work, postorder)
  end associate
end subroutine

!==============================================================================
! Build postorder array (simpler version for solve phase)
!==============================================================================
subroutine build_postorder_array(nsuper, sfils, sfrere, sparent, postorder)
  integer(kint), intent(in)  :: nsuper
  integer(kint), intent(in)  :: sfils(nsuper), sfrere(nsuper), sparent(nsuper)
  integer(kint), intent(out) :: postorder(nsuper)

  integer(kint) :: stack(nsuper)
  integer(kint) :: sp, s, cnt, child
  logical :: visited(nsuper)

  visited = .false.
  cnt = 0
  sp = 0

  ! Push roots
  do s = 1, nsuper
    if (sparent(s) == 0) then
      sp = sp + 1
      stack(sp) = s
    end if
  end do

  do while (sp > 0)
    s = stack(sp)
    child = sfils(s)

    if (child /= 0) then
      if (.not. visited(child)) then
        ! Push children
        do while (child /= 0)
          if (.not. visited(child)) then
            sp = sp + 1
            stack(sp) = child
          end if
          child = sfrere(child)
        end do
        cycle
      end if
    end if

    sp = sp - 1
    if (.not. visited(s)) then
      visited(s) = .true.
      cnt = cnt + 1
      postorder(cnt) = s
    end if
  end do
end subroutine

!==============================================================================
! Verify solution: compute ||Ax - b|| / ||b||
!==============================================================================
subroutine verify_solution(mat, lu, x, rhs_orig)
  type(matrix_data), intent(in) :: mat
  type(monolis_mat_lu), intent(in) :: lu
  real(kdouble), intent(in) :: x(mat%ndof*mat%n)
  real(kdouble), intent(in) :: rhs_orig(mat%ndof*mat%n)

  integer(kint) :: n, ndof, ntot
  real(kdouble), allocatable :: ax(:), residual(:)
  real(kdouble) :: norm_r, norm_b, rel_err, max_err
  integer(kint) :: k, i, r, c, id, jd, base

  n = mat%n; ndof = mat%ndof; ntot = n*ndof
  allocate(ax(ntot), residual(ntot))

  ! Compute A*x (block SpMV)
  ax = 0.0d0
  do r = 1, n
    do k = mat%row_ptr(r)+1, mat%row_ptr(r+1)
      c = mat%col_ind(k)
      base = (k-1)*ndof*ndof
      do id = 1, ndof
        do jd = 1, ndof
          ax((r-1)*ndof+id) = ax((r-1)*ndof+id) + &
            mat%a_elt(base + (id-1)*ndof + jd) * x((c-1)*ndof+jd)
        end do
      end do
    end do
  end do

  ! Compute residual = b - A*x
  residual = rhs_orig - ax

  ! Norms
  norm_r = 0.0d0
  norm_b = 0.0d0
  max_err = 0.0d0
  do i = 1, ntot
    norm_r = norm_r + residual(i)**2
    norm_b = norm_b + rhs_orig(i)**2
    if (abs(x(i) - 1.0d0) > max_err) max_err = abs(x(i) - 1.0d0)
  end do
  norm_r = sqrt(norm_r)
  norm_b = sqrt(norm_b)

  if (norm_b > 0.0d0) then
    rel_err = norm_r / norm_b
  else
    rel_err = norm_r
  end if

  write(*,'(A,ES12.4)') '  ||b - Ax|| / ||b|| = ', rel_err
  write(*,'(A,ES12.4)') '  ||b - Ax||         = ', norm_r
  write(*,'(A,ES12.4)') '  max|x(i) - 1|      = ', max_err

  deallocate(ax, residual)
end subroutine

end module mod_monolis_fact_factorize
