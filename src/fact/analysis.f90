module mod_monolis_fact_analysis
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
  ! Reverse Cuthill-McKee ordering
  ! O(N + nz*logN) — much better fill-reduction for FEM matrices
  ! than the simple greedy degree ordering.
  !
  ! 1. Find pseudo-peripheral starting node (2-round BFS refinement)
  ! 2. BFS with degree-ordered neighbor insertion
  ! 3. Reverse the resulting ordering
  !==============================================================================
  subroutine reverse_cuthill_mckee_ordering(mat, lu)
    type(matrix_data), intent(inout) :: mat
    type(monolis_mat_lu), intent(inout) :: lu

    integer(kint) :: n, i, k, j, r, c
    integer(kint), allocatable :: degree(:), queue(:), nbr_buf(:)
    logical, allocatable :: visited(:)
    integer(kint) :: head, tail, node, nbr, start
    integer(kint) :: min_deg, nbr_count, tmp

    n = mat%n
    allocate(mat%perm(n), mat%invperm(n))
    allocate(degree(n), visited(n), queue(n), nbr_buf(n))

    ! Compute degrees (off-diagonal only)
    degree = 0
    do r = 1, n
      do k = mat%row_ptr(r)+1, mat%row_ptr(r+1)
        c = mat%col_ind(k)
        if (r /= c) degree(r) = degree(r) + 1
      end do
    end do

    ! Find pseudo-peripheral starting node
    ! Start with minimum degree node, then refine with 2 BFS rounds
    min_deg = n + 1
    start = 1
    do i = 1, n
      if (degree(i) < min_deg .and. degree(i) > 0) then
        min_deg = degree(i)
        start = i
      end if
    end do

    ! Refine: 2 rounds of BFS; last node of BFS is pseudo-peripheral
    do j = 1, 2
      visited = .false.
      head = 1; tail = 1
      queue(1) = start
      visited(start) = .true.

      do while (head <= tail)
        node = queue(head); head = head + 1
        do k = mat%row_ptr(node)+1, mat%row_ptr(node+1)
          nbr = mat%col_ind(k)
          if (nbr /= node .and. .not. visited(nbr)) then
            visited(nbr) = .true.
            tail = tail + 1
            queue(tail) = nbr
          end if
        end do
      end do

      start = queue(tail)   ! last node = pseudo-peripheral
    end do

    ! Main BFS with degree-sorted neighbor insertion (Cuthill-McKee)
    visited = .false.
    head = 1; tail = 1
    queue(1) = start
    visited(start) = .true.

    do while (head <= tail)
      node = queue(head); head = head + 1

      ! Collect unvisited neighbors
      nbr_count = 0
      do k = mat%row_ptr(node)+1, mat%row_ptr(node+1)
        nbr = mat%col_ind(k)
        if (nbr /= node .and. .not. visited(nbr)) then
          nbr_count = nbr_count + 1
          nbr_buf(nbr_count) = nbr
        end if
      end do

      ! Sort neighbors by non-decreasing degree (insertion sort — small lists)
      do i = 2, nbr_count
        tmp = nbr_buf(i)
        j = i - 1
        do while (j >= 1)
          if (degree(nbr_buf(j)) <= degree(tmp)) exit
          nbr_buf(j+1) = nbr_buf(j)
          j = j - 1
        end do
        nbr_buf(j+1) = tmp
      end do

      ! Add sorted neighbors to queue
      do i = 1, nbr_count
        nbr = nbr_buf(i)
        if (.not. visited(nbr)) then
          visited(nbr) = .true.
          tail = tail + 1
          queue(tail) = nbr
        end if
      end do
    end do

    ! Handle disconnected components
    do i = 1, n
      if (.not. visited(i)) then
        tail = tail + 1
        queue(tail) = i
      end if
    end do

    ! Reverse: RCM = reverse of CM ordering
    do i = 1, n
      mat%perm(n - i + 1) = queue(i)
      mat%invperm(queue(i)) = n - i + 1
    end do

    deallocate(degree, visited, queue, nbr_buf)
  end subroutine

end module mod_monolis_fact_analysis
