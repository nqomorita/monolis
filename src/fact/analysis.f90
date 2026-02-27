module mod_monolis_fact_analysis
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc

  implicit none

contains

  ! Build elimination tree using Liu's algorithm (column-by-column)
  !
  ! For each column k (in permuted order 1..N), find all edges to columns j < k
  ! and set parent(root_of_subtree(j)) = k using Union-Find with path compression.
  !
  ! This requires building a CSR adjacency list from COO first.
  subroutine build_elimination_tree(mat, lu)
    type(matrix_data), intent(in) :: mat
    type(monolis_mat_lu), intent(inout) :: lu

    integer :: n, nz
    integer, allocatable :: ancestor(:)
    integer, allocatable :: adj_ptr(:), adj_list(:)
    integer, allocatable :: adj_count(:)
    integer :: k, pi, pj, tmp, r, c, pos, i

    n = mat%N; nz = mat%NZ
    allocate(lu%PARENT(n))

    associate(row_ptr => mat%ROW_PTR, col_ind => mat%COL_IND, &
              invperm => mat%INVPERM, perm => mat%PERM, parent => lu%PARENT)

    ! ---- Step 1: Build CSR adjacency for the permuted lower triangle ----
    ! For each edge (pi, pj) with pi < pj, store pi in adj_list of column pj
    ! (i.e., for column pj, the list of earlier columns connected to it)
    allocate(adj_count(n))
    adj_count = 0

    do r = 1, n
      do k = row_ptr(r), row_ptr(r+1)-1
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
      do k = row_ptr(r), row_ptr(r+1)-1
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

  ! A supernode is a maximal set of contiguous columns in the factor
  ! with the same sparsity structure (fundamental supernode).
  !
  ! Detection rule:  column j can be merged with column j-1 if
  !   parent(j-1) == j   AND   nonzero-count(j) == nonzero-count(j-1) - 1
  !
  ! We also apply relaxed merging (NEMIN) to allow a small amount of fill.
  subroutine identify_supernodes(mat, lu)
    type(matrix_data), intent(in) :: mat
    type(monolis_mat_lu), intent(inout) :: lu

    integer :: n, nz, nsuper
    integer, allocatable :: sbelong(:), sstart(:), ssize(:), sparent(:), sfsize(:)
    integer, allocatable :: col_count(:)  ! nonzero count in each column of L (permuted)
    integer, allocatable :: children_count(:)
    integer :: i, k, r, c, pi, pj, tmp
    integer :: NEMIN
    logical :: can_merge

    n = mat%N; nz = mat%NZ
    NEMIN = 16  ! relaxation parameter

    allocate(col_count(n), children_count(n), sbelong(n))
    col_count = 1   ! diagonal always present
    children_count = 0

    associate(row_ptr => mat%ROW_PTR, col_ind => mat%COL_IND, &
              perm => mat%PERM, invperm => mat%INVPERM, parent => lu%PARENT)

    ! ---- Step 1: Estimate column counts of L ----
    ! For each off-diagonal entry in the permuted lower triangle,
    ! the column with the smaller permuted index gets +1.
    ! This is an approximation (exact requires symbolic factorization).
    do r = 1, n
      do k = row_ptr(r), row_ptr(r+1)-1
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
          if (col_count(i-1) > NEMIN .and. col_count(i) > NEMIN) then
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

    ! ---- Step 4: Relaxed merging of small supernodes (NEMIN) ----
    call relax_supernodes(n, nsuper, sbelong, sstart, ssize, sparent, sfsize, &
                          parent, col_count, NEMIN)

    ! ---- Step 5: Compute frontal sizes ----
    ! frontal size = npiv + number of rows/cols that appear in the
    ! frontal matrix but are not pivoted (the "index" set beyond pivots)
    call compute_frontal_sizes(n, nz, row_ptr, col_ind, invperm, &
                               nsuper, sstart, ssize, sparent, sbelong, sfsize)

    end associate

    ! Store results into lu
    lu%NSUPER = nsuper
    lu%SNODE_BELONG = sbelong
    lu%SNODE_START = sstart
    lu%SNODE_SIZE = ssize
    lu%SNODE_PARENT = sparent
    lu%SNODE_FSIZE = sfsize

    deallocate(col_count, children_count)
  end subroutine

  ! Relaxed supernode merging: merge parent-child pairs when both are small
  subroutine relax_supernodes(n, nsuper, sbelong, sstart, ssize, sparent, sfsize, &
                              parent, col_count, nemin)
    integer, intent(in)    :: n, nemin
    integer, intent(inout) :: nsuper
    integer, intent(inout) :: sbelong(n)
    integer, intent(inout), allocatable :: sstart(:), ssize(:), sparent(:), sfsize(:)
    integer, intent(in)    :: parent(n), col_count(n)

    logical, allocatable :: merged(:)
    integer :: i, s, sp, new_nsuper
    integer, allocatable :: new_sstart(:), new_ssize(:), new_sparent(:), new_sfsize(:)
    integer, allocatable :: map(:)

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

  ! Compute frontal sizes for each supernode
  !
  ! The frontal size of supernode s is:
  !   nfront(s) = npiv(s) + |{ row indices in L columns of s that are > last pivot of s }|
  !
  ! We compute this by finding all row indices that appear in the columns
  ! belonging to supernode s (in the permuted matrix + fill from children).
  subroutine compute_frontal_sizes(n, nz, row_ptr, col_ind, invperm, &
                                   nsuper, sstart, ssize, sparent, sbelong, sfsize)
    integer, intent(in) :: n, nz
    integer, intent(in) :: row_ptr(n+1), col_ind(nz)
    integer, intent(in) :: invperm(n)
    integer, intent(in) :: nsuper
    integer, intent(in) :: sstart(nsuper), ssize(nsuper)
    integer, intent(in) :: sparent(nsuper), sbelong(n)
    integer, intent(inout) :: sfsize(nsuper)

    integer, allocatable :: row_set(:)
    logical, allocatable :: marker(:)
    integer :: s, k, pi, pj, tmp, r, c, cnt
    integer :: first_col, last_col, child

    allocate(marker(n), row_set(n))

    do s = 1, nsuper
      first_col = sstart(s)
      last_col  = sstart(s) + ssize(s) - 1
      marker = .false.
      cnt = 0

      ! Mark all row indices from original matrix entries
      ! that fall into columns [first_col .. last_col] in permuted order
      do r = 1, n
        do k = row_ptr(r), row_ptr(r+1)-1
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

  ! Build child-sibling representation for the supernode tree
  subroutine build_frontal_tree(mat, lu)
    type(matrix_data), intent(in) :: mat
    type(monolis_mat_lu), intent(inout) :: lu

    integer :: nsuper, s, p, cur

    nsuper = lu%NSUPER
    allocate(lu%SFILS(nsuper), lu%SFRERE(nsuper))
    lu%SFILS  = 0
    lu%SFRERE = 0

    associate(sparent => lu%SNODE_PARENT, sfils => lu%SFILS, sfrere => lu%SFRERE)

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

end module mod_monolis_fact_analysis
