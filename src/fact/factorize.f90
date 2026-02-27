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
! Process supernodes in postorder (leaves first, root last).
! Optimized with:
!   - Pre-built permuted CSC for O(nnz_local) index set construction
!   - MUMPS-style imap for O(1) position lookup
!   - CSC-based assembly (only scan relevant columns)
!   - Column-wise extend-add with DAXPY
!   - Quicksort for index sorting
!==============================================================================
subroutine multifrontal_factorize(mat, lu)
  type(matrix_data), intent(in) :: mat
  type(monolis_mat_lu), intent(inout) :: lu

  integer(kint) :: n, nz, nsuper, ndof, info
  integer(kint), allocatable :: postorder(:)
  integer(kint) :: s, ip, nfront, npiv, k, i, j
  integer(kint) :: first_col, last_col
  integer(kint) :: pi, pj, cnt, r, c, col, pos
  integer(kint), allocatable :: idx_set(:)
  integer(kint) :: child
  ! ---- MUMPS-style: global index map for O(1) position lookup ----
  integer(kint), allocatable :: imap(:)

  ! ---- Pre-built permuted CSC for O(nnz_local) assembly ----
  integer(kint), allocatable :: csc_ptr(:), csc_row(:), csc_origk(:)

  n = mat%n; nz = mat%nz; nsuper = lu%nsuper; ndof = mat%ndof
  info = 0

  allocate(lu%factors(nsuper))

  associate(row_ptr => mat%row_ptr, col_ind => mat%col_ind, a_elt => mat%a_elt, &
            perm => mat%perm, invperm => mat%invperm, &
            sstart => lu%snode_start, ssize => lu%snode_size, &
            sparent => lu%snode_parent, sfsize => lu%snode_fsize, &
            sfils => lu%sfils, sfrere => lu%sfrere, factors => lu%factors)

  call build_postorder(nsuper, sfils, sfrere, sparent, postorder)

  ! ---- Pre-build permuted CSC: for each permuted column, list (row, orig_k) ----
  ! This allows O(nnz_local) assembly per supernode instead of O(n*nz)
  call build_permuted_csc(n, nz, row_ptr, col_ind, invperm, &
                          csc_ptr, csc_row, csc_origk)

  allocate(imap(n))
  imap = 0

  ! ================================================================
  ! Process each supernode in postorder (leaves first, root last)
  ! ================================================================
  do ip = 1, nsuper
    s = postorder(ip)
    npiv = ssize(s)
    first_col = sstart(s)
    last_col  = first_col + npiv - 1

    ! ---- Determine index set using CSC (O(nnz_local)) ----
    cnt = 0
    allocate(idx_set(n))

    ! Pivot columns - mark via imap as temporary flag (use negative values)
    do i = first_col, last_col
      imap(i) = -1   ! temporary flag: "in pivot set"
    end do
    cnt = npiv
    idx_set(1:npiv) = (/ (i, i = first_col, last_col) /)

    ! Row indices from CSC columns in [first_col..last_col]
    do col = first_col, last_col
      do pos = csc_ptr(col), csc_ptr(col+1) - 1
        ! csc_row stores packed (pi, pj) pairs; get the "other" index (larger one)
        pi = csc_row(2*pos-1)
        pj = csc_row(2*pos)
        ! Check both indices: the one > last_col needs to be in the index set
        if (pi > last_col .and. imap(pi) == 0) then
          cnt = cnt + 1; idx_set(cnt) = pi; imap(pi) = -1
        end if
        if (pj > last_col .and. imap(pj) == 0) then
          cnt = cnt + 1; idx_set(cnt) = pj; imap(pj) = -1
        end if
      end do
    end do

    ! Indices from children's contribution blocks
    child = sfils(s)
    do while (child /= 0)
      do i = factors(child)%npiv + 1, factors(child)%nfront
        j = factors(child)%indices(i)
        if (j >= first_col .and. imap(j) == 0) then
          cnt = cnt + 1; idx_set(cnt) = j; imap(j) = -1
        end if
      end do
      child = sfrere(child)
    end do

    nfront = cnt

    ! Sort: pivots first (already in place), then others by quicksort
    if (nfront - npiv > 1) then
      call quicksort(idx_set(npiv+1:nfront), nfront - npiv)
    end if

    ! ---- Build global index map (MUMPS-style O(1) lookup) ----
    do i = 1, nfront
      imap(idx_set(i)) = i
    end do

    ! ---- Allocate frontal matrix ----
    factors(s)%nfront = nfront
    factors(s)%npiv   = npiv
    allocate(factors(s)%indices(nfront))
    allocate(factors(s)%front(nfront*ndof, nfront*ndof))
    allocate(factors(s)%pivorder(npiv))
    factors(s)%indices = idx_set(1:nfront)
    factors(s)%front   = 0.0d0
    factors(s)%pivorder = (/ (i, i = 1, npiv) /)
    deallocate(idx_set)

    ! ---- Assemble original entries using CSC + imap (O(nnz_local)) ----
    call assemble_original_csc(n, ndof, csc_ptr, csc_row, csc_origk, a_elt, &
                               nfront, factors(s)%front, &
                               first_col, last_col, imap)

    ! ---- Extend-add from children using imap (O(1) per element) ----
    child = sfils(s)
    do while (child /= 0)
      call extend_add_mapped(factors(child), factors(s), imap, ndof)
      child = sfrere(child)
    end do

    ! ---- BLAS Level-3 blocked lu factorization (MUMPS FAC_H/P style) ----
    call frontal_lu_factor_blas(factors(s)%front, nfront*ndof, npiv*ndof, info)
    if (info /= 0) then
      info = 0
    end if

    ! ---- Clear imap entries for this frontal ----
    do i = 1, nfront
      imap(factors(s)%indices(i)) = 0
    end do
  end do

  deallocate(imap, postorder, csc_ptr, csc_row, csc_origk)
  end associate
end subroutine

!==============================================================================
! Build permuted CSC: for each permuted column min(pi,pj), store (pi, pj, orig_k)
! Stores both original permuted indices to preserve assembly direction.
!==============================================================================
subroutine build_permuted_csc(n, nz, row_ptr, col_ind, invperm, &
                              csc_ptr, csc_row, csc_origk)
  integer(kint), intent(in) :: n, nz
  integer(kint), intent(in) :: row_ptr(n+1), col_ind(nz)
  integer(kint), intent(in) :: invperm(n)
  integer(kint), allocatable, intent(out) :: csc_ptr(:), csc_row(:), csc_origk(:)

  integer(kint), allocatable :: csc_count(:), csc_col2(:)
  integer(kint) :: r, k, c, pi, pj, nadj, mn

  allocate(csc_count(n))
  csc_count = 0

  ! Count: store each edge in the column with smaller permuted index
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

  ! csc_row stores the original pi (invperm(r)), csc_origk stores both orig k
  ! and we add csc_col2 for original pj (invperm(c))
  allocate(csc_row(nadj), csc_origk(nadj), csc_col2(nadj))
  csc_count = 0

  do r = 1, n
    pi = invperm(r)
    do k = row_ptr(r)+1, row_ptr(r+1)
      c = col_ind(k)
      pj = invperm(c)
      mn = min(pi, pj)
      csc_count(mn) = csc_count(mn) + 1
      ! Store original (pi, pj) so assembly direction is preserved
      csc_row(csc_ptr(mn) + csc_count(mn) - 1) = pi
      csc_col2(csc_ptr(mn) + csc_count(mn) - 1) = pj
      csc_origk(csc_ptr(mn) + csc_count(mn) - 1) = k
    end do
  end do

  ! Pack csc_col2 into csc_row by interleaving: store as (pi, pj) pairs
  ! Use separate arrays — csc_row = pi, csc_origk encodes k, extend with csc_col2
  ! For simplicity, store pj in a separate allocatable that we'll pass around.
  ! Actually, we'll encode differently: csc_row(pos) = pi, we need pj too.
  ! Solution: double the array and pack (pi, pj) pairs.

  ! Re-allocate csc_row to hold both pi and pj
  ! We'll use: csc_row(pos) = pi * (n+1) + pj as packed encoding won't work for large n
  ! Instead, simply reallocate csc_row to size 2*nadj: odd=pi, even=pj
  block
    integer(kint), allocatable :: temp(:)
    allocate(temp(2*nadj))
    do k = 1, nadj
      temp(2*k-1) = csc_row(k)    ! pi
      temp(2*k)   = csc_col2(k)   ! pj
    end do
    deallocate(csc_row)
    allocate(csc_row(2*nadj))
    csc_row = temp
    deallocate(temp)
  end block

  deallocate(csc_count, csc_col2)
end subroutine

!==============================================================================
! Assemble original entries using pre-built CSC (O(nnz_local) per supernode)
! Preserves original (pi, pj) direction for correct non-symmetric assembly.
!==============================================================================
subroutine assemble_original_csc(n, ndof, csc_ptr, csc_row, csc_origk, a_elt, &
                                  nfront, front, first_col, last_col, imap)
  integer(kint), intent(in) :: n, ndof
  integer(kint), intent(in) :: csc_ptr(n+1), csc_row(*), csc_origk(*)
  real(kdouble), intent(in) :: a_elt(*)
  integer(kint), intent(in) :: nfront
  real(kdouble), intent(inout) :: front(nfront*ndof, nfront*ndof)
  integer(kint), intent(in) :: first_col, last_col
  integer(kint), intent(in) :: imap(n)

  integer(kint) :: col, pos, pi, pj, ii, jj, id, jd, base, origk

  ! Scan only columns [first_col..last_col] in permuted CSC
  do col = first_col, last_col
    do pos = csc_ptr(col), csc_ptr(col+1) - 1
      ! Retrieve original (pi, pj) from packed storage
      pi = csc_row(2*pos-1)
      pj = csc_row(2*pos)
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
! Extend-add: add child's contribution block into parent's frontal matrix
!
! The contribution block is the Schur complement of the child's frontal
! matrix, i.e., the (nfront-npiv) x (nfront-npiv) lower-right block
! after lu elimination of the pivot rows/columns.
!==============================================================================
!==============================================================================
! Extend-add using precomputed imap + column-wise DAXPY (MUMPS-style)
! Parent's imap maps child's CB indices directly to parent positions.
! Column-wise traversal for cache-friendly access + DAXPY for inner loop.
!==============================================================================
subroutine extend_add_mapped(child_fac, parent_fac, imap, ndof)
  type(monolis_mat_frontal), intent(in)    :: child_fac
  type(monolis_mat_frontal), intent(inout) :: parent_fac
  integer(kint), intent(in) :: imap(*)
  integer(kint), intent(in) :: ndof

  integer(kint) :: ic, jc, ip, jp, id, jd
  integer(kint) :: npiv_c, nfront_c, ncb, nf_p
  integer(kint), allocatable :: cb_map(:)
  integer(kint) :: ncb_valid, iv
  integer(kint), allocatable :: valid_ic(:), valid_ip(:)

  npiv_c  = child_fac%npiv
  nfront_c = child_fac%nfront
  ncb = nfront_c - npiv_c
  nf_p = parent_fac%nfront
  if (ncb <= 0) return

  ! Pre-compute CB index mapping and filter valid entries
  allocate(cb_map(ncb), valid_ic(ncb), valid_ip(ncb))
  ncb_valid = 0
  do ic = 1, ncb
    cb_map(ic) = imap(child_fac%indices(npiv_c + ic))
    if (cb_map(ic) > 0) then
      ncb_valid = ncb_valid + 1
      valid_ic(ncb_valid) = ic
      valid_ip(ncb_valid) = cb_map(ic)
    end if
  end do
  if (ncb_valid == 0) then
    deallocate(cb_map, valid_ic, valid_ip)
    return
  end if

  ! Column-wise scatter-add (cache-friendly: access parent_fac%front by columns)
  do jc = 1, ncb
    jp = cb_map(jc)
    if (jp == 0) cycle
    do jd = 1, ndof
      do iv = 1, ncb_valid
        ic = valid_ic(iv)
        ip = valid_ip(iv)
        ! DAXPY-style: add child CB column segment to parent
        call daxpy(ndof, 1.0d0, &
                   child_fac%front((npiv_c+ic-1)*ndof+1, (npiv_c+jc-1)*ndof+jd), 1, &
                   parent_fac%front((ip-1)*ndof+1, (jp-1)*ndof+jd), 1)
      end do
    end do
  end do

  deallocate(cb_map, valid_ic, valid_ip)
end subroutine

!==============================================================================
! Frontal lu factorization with partial pivoting
!
! Performs lu decomposition on the first npiv columns of the nfront x nfront
! frontal matrix, leaving the Schur complement in the lower-right block.
!
! After this routine:
!   front(1:npiv, 1:npiv) contains U in the upper triangle and L below
!   front(npiv+1:nfront, 1:npiv) contains the remaining L entries
!   front(1:npiv, npiv+1:nfront) contains the remaining U entries
!   front(npiv+1:nfront, npiv+1:nfront) is the Schur complement (contribution block)
!==============================================================================
!==============================================================================
! BLAS Level-3 blocked lu factorization (MUMPS FAC_H / FAC_P style)
!
! Uses a panel-based approach:
!   1. Factor a panel of lkjib columns with Level-2 BLAS (DGER)
!   2. Apply panel to remaining columns with DTRSM + DGEMM (Level-3)
!
! After this routine, front(1:npiv, 1:npiv) stores L\U,
! front(npiv+1:nfront, 1:npiv) stores L below pivots,
! front(1:npiv, npiv+1:nfront) stores U to the right,
! front(npiv+1:nfront, npiv+1:nfront) is the Schur complement (CB).
!==============================================================================
subroutine frontal_lu_factor_blas(front, nfront, npiv, info)
  real(kdouble), intent(inout) :: front(nfront, nfront)
  integer(kint), intent(in)    :: nfront, npiv
  integer(kint), intent(inout) :: info

  integer(kint) :: k, i, j, kb, ke, nb, nel, nupd
  real(kdouble) :: pivot
  real(kdouble), parameter :: eps = 1.0d-14
  real(kdouble), parameter :: one = 1.0d0, alpha = -1.0d0
  ! MUMPS-style blocking factor (lkjib)
  integer(kint), parameter :: lkjib = 64

  info = 0

  ! ================================================================
  ! Panel-based blocked factorization
  ! Process NPIV columns in blocks of size lkjib
  ! ================================================================
  kb = 1
  do while (kb <= npiv)
    ke = min(kb + lkjib - 1, npiv)   ! end of this panel
    nb = ke - kb + 1                  ! panel width

    ! ---- Factor this panel with Level-2 operations ----
    ! (MUMPS FAC_H / FAC_M / FAC_N style: column-by-column with DGER)
    do k = kb, ke
      pivot = front(k, k)
      if (abs(pivot) < eps) then
        front(k, k) = sign(eps, pivot + eps)
        pivot = front(k, k)
        info = k
      end if

      ! Scale L column k: front(k+1:nfront, k) /= pivot
      nel = nfront - k
      if (nel > 0) then
        call dscal(nel, one/pivot, front(k+1, k), 1)
      end if

      ! Rank-1 update within the panel: columns k+1..ke only
      ! front(k+1:nfront, k+1:ke) -= L(:,k) * U(k, k+1:ke)
      nupd = ke - k   ! number of remaining columns in panel
      if (nupd > 0 .and. nel > 0) then
        call dger(nel, nupd, alpha, front(k+1, k), 1, &
                  front(k, k+1), nfront, front(k+1, k+1), nfront)
      end if
    end do

    ! ---- Apply panel to trailing columns with Level-3 BLAS ----
    ! (MUMPS FAC_P style: DTRSM + DGEMM)
    nupd = nfront - ke   ! remaining columns/rows to update
    if (nupd > 0 .and. nb > 0) then
      nel = nfront - ke

      ! Step 1: Solve for U block using panel's L factor
      ! L(kb:ke, kb:ke) * U_new(kb:ke, ke+1:nfront) = front(kb:ke, ke+1:nfront)
      ! L is unit-lower-triangular (diagonal of L = 1), stored below diagonal
      call dtrsm('L', 'L', 'N', 'U', nb, nel, one, &
                 front(kb, kb), nfront, front(kb, ke+1), nfront)

      ! Step 2: Schur complement update
      ! front(ke+1:nfront, ke+1:nfront) -= L(ke+1:nfront, kb:ke) * U(kb:ke, ke+1:nfront)
      call dgemm('N', 'N', nel, nel, nb, alpha, &
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
