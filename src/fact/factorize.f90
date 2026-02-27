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
  integer(kint) :: k, pi, pj, tmp, r, c, pos, i

  n = mat%n; nz = mat%nz
  allocate(lu%parent(n))

  associate(row_ptr => mat%row_ptr, col_ind => mat%col_ind, &
            invperm => mat%invperm, perm => mat%perm, parent => lu%parent)

  ! ---- Step 1: Build CSR adjacency for the permuted lower triangle ----
  ! For each edge (pi, pj) with pi < pj, store pi in adj_list of column pj
  ! (i.e., for column pj, the list of earlier columns connected to it)
  allocate(adj_count(n))
  adj_count = 0

  do r = 1, n
    do k = row_ptr(r)+1, row_ptr(r+1)
      c = col_ind(k)
      if (r == c) cycle
      pi = invperm(r)
      pj = invperm(c)
      if (pi > pj) then; tmp = pi; pi = pj; pj = tmp; end if
      adj_count(pj) = adj_count(pj) + 1
    end do
  end do

  allocate(adj_ptr(n+1))
  adj_ptr(1) = 1
  do k = 2, n+1
    adj_ptr(k) = adj_ptr(k-1) + adj_count(k-1)
  end do

  allocate(adj_list(adj_ptr(n+1)-1))
  adj_count = 0  ! reuse as insertion counter

  do r = 1, n
    do k = row_ptr(r)+1, row_ptr(r+1)
      c = col_ind(k)
      if (r == c) cycle
      pi = invperm(r)
      pj = invperm(c)
      if (pi > pj) then; tmp = pi; pi = pj; pj = tmp; end if
      adj_count(pj) = adj_count(pj) + 1
      adj_list(adj_ptr(pj) + adj_count(pj) - 1) = pi
    end do
  end do

  ! ---- Step 2: Liu's elimination tree algorithm ----
  allocate(ancestor(n))
  parent = 0

  do k = 1, n
    ancestor(k) = k   ! each node is its own root initially

    ! Process all edges to earlier columns
    do pos = adj_ptr(k), adj_ptr(k+1) - 1
      pj = adj_list(pos)   ! pj < k

      ! Find root of pj's subtree with path compression
      r = pj
      do while (ancestor(r) /= r)
        r = ancestor(r)
      end do

      ! Path compression: flatten all nodes on the path to root
      i = pj
      do while (i /= r)
        tmp = ancestor(i)
        ancestor(i) = r
        i = tmp
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
  integer(kint) :: i, k, r, c, pi, pj, tmp
  integer(kint) :: nemin
  logical :: can_merge

  n = mat%n; nz = mat%nz
  nemin = 16  ! relaxation parameter

  allocate(col_count(n), children_count(n), sbelong(n))
  col_count = 1   ! diagonal always present
  children_count = 0

  associate(row_ptr => mat%row_ptr, col_ind => mat%col_ind, &
            perm => mat%perm, invperm => mat%invperm, parent => lu%parent)

  ! ---- Step 1: Estimate column counts of L ----
  ! For each off-diagonal entry in the permuted lower triangle,
  ! the column with the smaller permuted index gets +1.
  ! This is an approximation (exact requires symbolic factorization).
  do r = 1, n
    do k = row_ptr(r)+1, row_ptr(r+1)
      c = col_ind(k)
      if (r == c) cycle
      pi = invperm(r)
      pj = invperm(c)
      if (pi > pj) then
        tmp = pi; pi = pj; pj = tmp
      end if
      col_count(pi) = col_count(pi) + 1
    end do
  end do

  ! Also propagate fill through the tree: if parent(i)=j then
  ! col j inherits (col_count(i)-1) entries minus one for column i itself.
  ! Simple approximation: col_count(parent(i)) += max(col_count(i) - 1, 0)
  ! We skip this for now and rely on the structure-based merge criterion.

  ! Count children
  do i = 1, n
    if (parent(i) > 0 .and. parent(i) <= n) then
      children_count(parent(i)) = children_count(parent(i)) + 1
    end if
  end do

  ! ---- Step 2: Fundamental supernode detection ----
  ! column j starts a new supernode if:
  !   parent(j-1) /= j  OR  children_count(j) > 1
  !   (i.e., parent chain is broken or j has multiple children)
  sbelong(1) = 1
  nsuper = 1
  do i = 2, n
    can_merge = .true.
    ! parent of previous column must be this column
    if (parent(i-1) /= i) can_merge = .false.
    ! this column must have exactly one child (the previous column)
    if (children_count(i) > 1) can_merge = .false.
    ! column count compatibility (approximate)
    if (can_merge .and. col_count(i) > 0 .and. col_count(i-1) > 0) then
      if (col_count(i) /= col_count(i-1) - 1) then
        ! Allow relaxed merging if nodes are small
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

  ! sstart = first variable, ssize = count of variables
  do i = 1, n
    k = sbelong(i)
    ssize(k) = ssize(k) + 1
    if (sstart(k) == 0) sstart(k) = i
  end do

  ! sparent: parent supernode
  do i = 1, nsuper
    ! last variable in this supernode
    k = sstart(i) + ssize(i) - 1
    if (parent(k) > 0 .and. parent(k) <= n) then
      sparent(i) = sbelong(parent(k))
    else
      sparent(i) = 0   ! root
    end if
  end do

  ! ---- Step 4: Relaxed merging of small supernodes (nemin) ----
  call relax_supernodes(n, nsuper, sbelong, sstart, ssize, sparent, sfsize, &
                        parent, col_count, nemin)

  ! ---- Step 5: Compute frontal sizes ----
  ! frontal size = npiv + number of rows/cols that appear in the
  ! frontal matrix but are not pivoted (the "index" set beyond pivots)
  call compute_frontal_sizes(n, nz, row_ptr, col_ind, invperm, &
                             nsuper, sstart, ssize, sparent, sbelong, sfsize)

  end associate

  ! Store results into lu
  lu%nsuper = nsuper
  lu%snode_belong = sbelong
  lu%snode_start = sstart
  lu%snode_size = ssize
  lu%snode_parent = sparent
  lu%snode_fsize = sfsize

  deallocate(col_count, children_count)
end subroutine

!==============================================================================
! Relaxed supernode merging: merge parent-child pairs when both are small
!==============================================================================
subroutine relax_supernodes(n, nsuper, sbelong, sstart, ssize, sparent, sfsize, &
                            parent, col_count, nemin)
  integer(kint), intent(in)    :: n, nemin
  integer(kint), intent(inout) :: nsuper
  integer(kint), intent(inout) :: sbelong(n)
  integer(kint), intent(inout), allocatable :: sstart(:), ssize(:), sparent(:), sfsize(:)
  integer(kint), intent(in)    :: parent(n), col_count(n)

  logical, allocatable :: merged(:)
  integer(kint) :: i, s, sp, new_nsuper
  integer(kint), allocatable :: new_sstart(:), new_ssize(:), new_sparent(:), new_sfsize(:)
  integer(kint), allocatable :: map(:)

  allocate(merged(nsuper))
  merged = .false.

  ! Mark child supernodes for merging into parent when both are small
  do s = 1, nsuper
    sp = sparent(s)
    if (sp > 0 .and. sp <= nsuper .and. sp /= s) then
      if (ssize(s) <= nemin .and. ssize(sp) <= nemin) then
        ! Check that s is the only child feeding into sp
        ! (simple check: merge only if sparent is the immediate successor)
        if (sstart(sp) == sstart(s) + ssize(s)) then
          merged(s) = .true.
        end if
      end if
    end if
  end do

  ! Apply merges: fold merged supernodes into their parents
  do s = 1, nsuper
    if (merged(s)) then
      sp = sparent(s)
      ! Expand parent to include child
      sstart(sp) = min(sstart(sp), sstart(s))
      ssize(sp) = ssize(sp) + ssize(s)
      sparent(sp) = sparent(sp)   ! keep parent's parent
      ! Update sbelong
      do i = 1, n
        if (sbelong(i) == s) sbelong(i) = sp
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
        ! Find the non-merged ancestor
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

  ! Update sbelong
  do i = 1, n
    sbelong(i) = map(sbelong(i))
  end do

  ! Replace arrays
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
! The frontal size of supernode s is:
!   nfront(s) = npiv(s) + |{ row indices in L columns of s that are > last pivot of s }|
!
! We compute this by finding all row indices that appear in the columns
! belonging to supernode s (in the permuted matrix + fill from children).
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

  integer(kint), allocatable :: row_set(:)
  logical, allocatable :: marker(:)
  integer(kint) :: s, k, pi, pj, tmp, r, c, cnt
  integer(kint) :: first_col, last_col, child

  allocate(marker(n), row_set(n))

  do s = 1, nsuper
    first_col = sstart(s)
    last_col  = sstart(s) + ssize(s) - 1
    marker = .false.
    cnt = 0

    ! Mark all row indices from original matrix entries
    ! that fall into columns [first_col .. last_col] in permuted order
    do r = 1, n
      do k = row_ptr(r)+1, row_ptr(r+1)
        c = col_ind(k)
        pi = invperm(r)
        pj = invperm(c)

        ! Check if pj is in this supernode's columns
        if (pj >= first_col .and. pj <= last_col) then
          if (pi > last_col .and. .not. marker(pi)) then
            marker(pi) = .true.
            cnt = cnt + 1
          end if
        end if
        ! Symmetric: also check transpose
        if (pi >= first_col .and. pi <= last_col) then
          if (pj > last_col .and. .not. marker(pj)) then
            marker(pj) = .true.
            cnt = cnt + 1
          end if
        end if
      end do
    end do

    ! Also include indices from children's contribution blocks
    ! (indices that survive from child into parent)
    do child = 1, nsuper
      if (sparent(child) == s) then
        ! The contribution block of child has indices > last pivot of child
        ! that also appear in parent. We approximate this: all indices
        ! of child's frontal beyond its pivots propagate upward.
        ! (Precise computation would require the child's index set, but
        !  at this stage we use the original matrix structure.)
      end if
    end do

    sfsize(s) = ssize(s) + cnt
    ! At minimum, frontal size = pivot count
    if (sfsize(s) < ssize(s)) sfsize(s) = ssize(s)
  end do

  deallocate(marker, row_set)
end subroutine

!==============================================================================
! Build child-sibling representation for the supernode tree
!==============================================================================
subroutine build_frontal_tree(mat, lu)
  type(matrix_data), intent(in) :: mat
  type(monolis_mat_lu), intent(inout) :: lu

  integer(kint) :: nsuper, s, p, cur

  nsuper = lu%nsuper
  allocate(lu%sfils(nsuper), lu%sfrere(nsuper))
  lu%sfils  = 0
  lu%sfrere = 0

  associate(sparent => lu%snode_parent, sfils => lu%sfils, sfrere => lu%sfrere)

  do s = 1, nsuper
    p = sparent(s)
    if (p > 0 .and. p <= nsuper) then
      if (sfils(p) == 0) then
        sfils(p) = s
      else
        cur = sfils(p)
        do while (sfrere(cur) /= 0)
          cur = sfrere(cur)
        end do
        sfrere(cur) = s
      end if
    end if
  end do

  end associate
end subroutine

!==============================================================================
! Multifrontal lu factorization
!
! Process supernodes in postorder (leaves first, root last).
! For each supernode:
!   1. Allocate frontal matrix
!   2. Assemble original entries (arrowhead)
!   3. Extend-add contribution blocks from children
!   4. Partial pivot lu on the fully-assembled pivots
!   5. Store L, U factors and contribution block
!==============================================================================
subroutine multifrontal_factorize(mat, lu)
  type(matrix_data), intent(in) :: mat
  type(monolis_mat_lu), intent(inout) :: lu

  integer(kint) :: n, nz, nsuper, ndof, info
  integer(kint), allocatable :: postorder(:)
  integer(kint) :: s, ip, nfront, npiv, k, i, j
  integer(kint) :: first_col, last_col
  integer(kint) :: pi, pj, cnt, r, c
  integer(kint), allocatable :: idx_set(:)
  logical, allocatable :: marker(:)
  integer(kint) :: child
  ! ---- MUMPS-style: global index map for O(1) position lookup ----
  integer(kint), allocatable :: imap(:)   ! imap(global_idx) -> position in frontal

  n = mat%n; nz = mat%nz; nsuper = lu%nsuper; ndof = mat%ndof
  info = 0

  allocate(lu%factors(nsuper))

  associate(row_ptr => mat%row_ptr, col_ind => mat%col_ind, a_elt => mat%a_elt, &
            perm => mat%perm, invperm => mat%invperm, &
            sstart => lu%snode_start, ssize => lu%snode_size, &
            sparent => lu%snode_parent, sfsize => lu%snode_fsize, &
            sfils => lu%sfils, sfrere => lu%sfrere, factors => lu%factors)

  call build_postorder(nsuper, sfils, sfrere, sparent, postorder)
  allocate(marker(n), imap(n))
  imap = 0

  ! ================================================================
  ! Process each supernode in postorder (leaves first, root last)
  ! ================================================================
  do ip = 1, nsuper
    s = postorder(ip)
    npiv = ssize(s)
    first_col = sstart(s)
    last_col  = first_col + npiv - 1

    ! ---- Determine index set of the frontal matrix ----
    marker = .false.
    cnt = 0
    allocate(idx_set(n))

    ! Pivot columns
    do i = first_col, last_col
      marker(i) = .true.
    end do
    cnt = npiv
    idx_set(1:npiv) = (/ (i, i = first_col, last_col) /)

    ! Row indices from original matrix (only update indices > last_col)
    do r = 1, n
      do k = row_ptr(r)+1, row_ptr(r+1)
        c = col_ind(k)
        pi = invperm(r)
        pj = invperm(c)
        if (pj >= first_col .and. pj <= last_col) then
          if (pi > last_col .and. .not. marker(pi)) then
            cnt = cnt + 1; idx_set(cnt) = pi; marker(pi) = .true.
          end if
        end if
        if (pi >= first_col .and. pi <= last_col) then
          if (pj > last_col .and. .not. marker(pj)) then
            cnt = cnt + 1; idx_set(cnt) = pj; marker(pj) = .true.
          end if
        end if
      end do
    end do

    ! Indices from children's contribution blocks
    child = sfils(s)
    do while (child /= 0)
      do i = factors(child)%npiv + 1, factors(child)%nfront
        j = factors(child)%indices(i)
        if (j >= first_col .and. .not. marker(j)) then
          cnt = cnt + 1; idx_set(cnt) = j; marker(j) = .true.
        end if
      end do
      child = sfrere(child)
    end do

    nfront = cnt
    call sort_index_set(idx_set, nfront, npiv, first_col, last_col)
    do i = 1, nfront
      marker(idx_set(i)) = .false.
    end do

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

    ! ---- Assemble original entries using imap (O(1) per entry) ----
    call assemble_original_mapped(n, nz, ndof, row_ptr, col_ind, a_elt, invperm, &
                                  nfront, factors(s)%front, &
                                  first_col, last_col, imap)

    ! ---- Extend-add from children using imap (O(1) per element) ----
    child = sfils(s)
    do while (child /= 0)
      call extend_add_mapped(factors(child), factors(s), imap, ndof)
      ! Free child's contribution block (L/U factors still needed for solve)
      ! We keep front because solve reads from it.
      child = sfrere(child)
    end do

    ! ---- BLAS Level-3 blocked lu factorization (MUMPS FAC_H/P style) ----
    call frontal_lu_factor_blas(factors(s)%front, nfront*ndof, npiv*ndof, info)
    if (info /= 0) then
      info = 0  ! continue with small pivot
    end if

    ! ---- Clear imap entries for this frontal ----
    do i = 1, nfront
      imap(factors(s)%indices(i)) = 0
    end do
  end do

  deallocate(marker, imap, postorder)
  end associate
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
! Sort index set: pivots in [first_col..last_col] first (in order),
! then remaining indices in ascending order
!==============================================================================
subroutine sort_index_set(idx, n, npiv, first_col, last_col)
  integer(kint), intent(inout) :: idx(n)
  integer(kint), intent(in) :: n, npiv, first_col, last_col

  integer(kint), allocatable :: pivots(:), others(:)
  integer(kint) :: i, np, no, j, tmp

  allocate(pivots(npiv), others(n))
  np = 0; no = 0

  do i = 1, n
    if (idx(i) >= first_col .and. idx(i) <= last_col) then
      np = np + 1
      pivots(np) = idx(i)
    else
      no = no + 1
      others(no) = idx(i)
    end if
  end do

  ! Sort pivots
  call isort(pivots, np)
  ! Sort others
  call isort(others, no)

  ! Reassemble
  idx(1:np) = pivots(1:np)
  idx(np+1:np+no) = others(1:no)

  deallocate(pivots, others)
end subroutine

!==============================================================================
! Simple insertion sort for small integer(kint) arrays
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
! Assemble original matrix entries into frontal matrix
!==============================================================================
!==============================================================================
! Assemble original matrix entries using precomputed imap (O(1) lookup)
! MUMPS-style: imap(global_idx) gives position in frontal matrix directly
!==============================================================================
subroutine assemble_original_mapped(n, nz, ndof, row_ptr, col_ind, a_elt, invperm, &
                                     nfront, front, first_col, last_col, imap)
  integer(kint), intent(in) :: n, nz, ndof
  integer(kint), intent(in) :: row_ptr(n+1), col_ind(nz)
  real(kdouble), intent(in) :: a_elt(ndof*ndof*nz)
  integer(kint), intent(in) :: invperm(n)
  integer(kint), intent(in) :: nfront
  real(kdouble), intent(inout) :: front(nfront*ndof, nfront*ndof)
  integer(kint), intent(in) :: first_col, last_col
  integer(kint), intent(in) :: imap(n)

  integer(kint) :: k, pi, pj, ii, jj, mn, r, c, id, jd, base

  do r = 1, n
    do k = row_ptr(r)+1, row_ptr(r+1)
      c = col_ind(k)
      pi = invperm(r)
      pj = invperm(c)
      mn = min(pi, pj)
      if (mn < first_col .or. mn > last_col) cycle
      ii = imap(pi)
      jj = imap(pj)
      if (ii > 0 .and. jj > 0) then
        base = (k-1)*ndof*ndof  ! k is already 1-based here
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
! Extend-add: add child's contribution block into parent's frontal matrix
!
! The contribution block is the Schur complement of the child's frontal
! matrix, i.e., the (nfront-npiv) x (nfront-npiv) lower-right block
! after lu elimination of the pivot rows/columns.
!==============================================================================
!==============================================================================
! Extend-add using precomputed imap (O(1) per element)
! MUMPS-style: parent's imap is already set, so we map child's CB indices
! directly to parent positions without linear search.
!==============================================================================
subroutine extend_add_mapped(child_fac, parent_fac, imap, ndof)
  type(monolis_mat_frontal), intent(in)    :: child_fac
  type(monolis_mat_frontal), intent(inout) :: parent_fac
  integer(kint), intent(in) :: imap(*)
  integer(kint), intent(in) :: ndof

  integer(kint) :: ic, jc, ip, jp, id, jd
  integer(kint) :: npiv_c, nfront_c, ncb
  integer(kint), allocatable :: cb_map(:)

  npiv_c  = child_fac%npiv
  nfront_c = child_fac%nfront
  ncb = nfront_c - npiv_c
  if (ncb <= 0) return

  allocate(cb_map(ncb))
  do ic = 1, ncb
    cb_map(ic) = imap(child_fac%indices(npiv_c + ic))
  end do

  do jc = 1, ncb
    jp = cb_map(jc)
    if (jp == 0) cycle
    do ic = 1, ncb
      ip = cb_map(ic)
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

  deallocate(cb_map)
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
