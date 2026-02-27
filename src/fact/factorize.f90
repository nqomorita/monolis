module mod_monolis_fact_factorize
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc

  implicit none

contains

  ! Multifrontal LU factorization
  !
  ! Process supernodes in postorder (leaves first, root last).
  ! For each supernode:
  !   1. Allocate frontal matrix
  !   2. Assemble original entries (arrowhead)
  !   3. Extend-add contribution blocks from children
  !   4. Partial pivot LU on the fully-assembled pivots
  !   5. Store L, U factors and contribution block
  subroutine multifrontal_factorize(mat, lu)
    type(matrix_data), intent(in) :: mat
    type(monolis_mat_lu), intent(inout) :: lu

    integer :: n, nz, nsuper, ndof, info
    integer, allocatable :: postorder(:)
    integer :: s, ip, nfront, npiv, k, i, j
    integer :: first_col, last_col
    integer :: pi, pj, cnt, r, c
    integer, allocatable :: idx_set(:)
    logical, allocatable :: marker(:)
    integer :: child
    ! ---- MUMPS-style: global index map for O(1) position lookup ----
    integer, allocatable :: imap(:)   ! imap(global_idx) -> position in frontal

    n = mat%N; nz = mat%NZ; nsuper = lu%NSUPER; ndof = mat%NDOF
    info = 0

    allocate(lu%FACTORS(nsuper))

    associate(row_ptr => mat%ROW_PTR, col_ind => mat%COL_IND, a_elt => mat%A_elt, &
              perm => mat%PERM, invperm => mat%INVPERM, &
              sstart => lu%SNODE_START, ssize => lu%SNODE_SIZE, &
              sparent => lu%SNODE_PARENT, sfsize => lu%SNODE_FSIZE, &
              sfils => lu%SFILS, sfrere => lu%SFRERE, factors => lu%FACTORS)

    call build_postorder(nsuper, sfils, sfrere, sparent, postorder)
    allocate(marker(n), imap(n))
    imap = 0

    ! Process each supernode in postorder (leaves first, root last)
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
        do k = row_ptr(r), row_ptr(r+1)-1
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

      ! ---- BLAS Level-3 blocked LU factorization (MUMPS FAC_H/P style) ----
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

  ! Build postorder traversal of the supernode tree
  subroutine build_postorder(nsuper, sfils, sfrere, sparent, postorder)
    integer, intent(in)  :: nsuper
    integer, intent(in)  :: sfils(nsuper), sfrere(nsuper), sparent(nsuper)
    integer, allocatable, intent(out) :: postorder(:)

    integer, allocatable :: stack(:)
    integer :: sp, s, cnt, child
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

  ! Sort index set: pivots in [first_col..last_col] first (in order),
  ! then remaining indices in ascending order
  subroutine sort_index_set(idx, n, npiv, first_col, last_col)
    integer, intent(inout) :: idx(n)
    integer, intent(in) :: n, npiv, first_col, last_col

    integer, allocatable :: pivots(:), others(:)
    integer :: i, np, no, j, tmp

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

  ! Simple insertion sort for small integer arrays
  subroutine isort(arr, n)
    integer, intent(inout) :: arr(n)
    integer, intent(in) :: n
    integer :: i, j, key

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

  ! Assemble original matrix entries into frontal matrix
  ! Assemble original matrix entries using precomputed imap (O(1) lookup)
  ! MUMPS-style: imap(global_idx) gives position in frontal matrix directly
  subroutine assemble_original_mapped(n, nz, ndof, row_ptr, col_ind, a_elt, invperm, &
                                       nfront, front, first_col, last_col, imap)
    integer, intent(in) :: n, nz, ndof
    integer, intent(in) :: row_ptr(n+1), col_ind(nz)
    double precision, intent(in) :: a_elt(ndof*ndof*nz)
    integer, intent(in) :: invperm(n)
    integer, intent(in) :: nfront
    double precision, intent(inout) :: front(nfront*ndof, nfront*ndof)
    integer, intent(in) :: first_col, last_col
    integer, intent(in) :: imap(n)

    integer :: k, pi, pj, ii, jj, mn, r, c, id, jd, base

    do r = 1, n
      do k = row_ptr(r), row_ptr(r+1)-1
        c = col_ind(k)
        pi = invperm(r)
        pj = invperm(c)
        mn = min(pi, pj)
        if (mn < first_col .or. mn > last_col) cycle
        ii = imap(pi)
        jj = imap(pj)
        if (ii > 0 .and. jj > 0) then
          base = (k-1)*ndof*ndof
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

  ! Extend-add: add child's contribution block into parent's frontal matrix
  !
  ! The contribution block is the Schur complement of the child's frontal
  ! matrix, i.e., the (nfront-npiv) x (nfront-npiv) lower-right block
  ! after LU elimination of the pivot rows/columns.
  ! Extend-add using precomputed imap (O(1) per element)
  ! MUMPS-style: parent's imap is already set, so we map child's CB indices
  ! directly to parent positions without linear search.
  subroutine extend_add_mapped(child_fac, parent_fac, imap, ndof)
    type(monolis_mat_frontal), intent(in)    :: child_fac
    type(monolis_mat_frontal), intent(inout) :: parent_fac
    integer, intent(in) :: imap(*)
    integer, intent(in) :: ndof

    integer :: ic, jc, ip, jp, id, jd
    integer :: npiv_c, nfront_c, ncb
    integer, allocatable :: cb_map(:)

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

  ! Frontal LU factorization with partial pivoting
  !
  ! Performs LU decomposition on the first npiv columns of the nfront x nfront
  ! frontal matrix, leaving the Schur complement in the lower-right block.
  !
  ! After this routine:
  !   front(1:npiv, 1:npiv) contains U in the upper triangle and L below
  !   front(npiv+1:nfront, 1:npiv) contains the remaining L entries
  !   front(1:npiv, npiv+1:nfront) contains the remaining U entries
  !   front(npiv+1:nfront, npiv+1:nfront) is the Schur complement (contribution block)
  ! BLAS Level-3 blocked LU factorization (MUMPS FAC_H / FAC_P style)
  !
  ! Uses a panel-based approach:
  !   1. Factor a panel of LKJIB columns with Level-2 BLAS (DGER)
  !   2. Apply panel to remaining columns with DTRSM + DGEMM (Level-3)
  !
  ! After this routine, front(1:npiv, 1:npiv) stores L\U,
  ! front(npiv+1:nfront, 1:npiv) stores L below pivots,
  ! front(1:npiv, npiv+1:nfront) stores U to the right,
  ! front(npiv+1:nfront, npiv+1:nfront) is the Schur complement (CB).
  subroutine frontal_lu_factor_blas(front, nfront, npiv, info)
    double precision, intent(inout) :: front(nfront, nfront)
    integer, intent(in)    :: nfront, npiv
    integer, intent(inout) :: info

    integer :: k, i, j, kb, ke, nb, nel, nupd
    double precision :: pivot
    double precision, parameter :: EPS = 1.0d-14
    double precision, parameter :: ONE = 1.0d0, ALPHA = -1.0d0
    ! MUMPS-style blocking factor (LKJIB)
    integer, parameter :: LKJIB = 64

    info = 0

    ! Panel-based blocked factorization
    ! Process NPIV columns in blocks of size LKJIB
    kb = 1
    do while (kb <= npiv)
      ke = min(kb + LKJIB - 1, npiv)   ! end of this panel
      nb = ke - kb + 1                  ! panel width

      ! ---- Factor this panel with Level-2 operations ----
      ! (MUMPS FAC_H / FAC_M / FAC_N style: column-by-column with DGER)
      do k = kb, ke
        pivot = front(k, k)
        if (abs(pivot) < EPS) then
          front(k, k) = sign(EPS, pivot + EPS)
          pivot = front(k, k)
          info = k
        end if

        ! Scale L column k: front(k+1:nfront, k) /= pivot
        nel = nfront - k
        if (nel > 0) then
          call dscal(nel, ONE/pivot, front(k+1, k), 1)
        end if

        ! Rank-1 update within the panel: columns k+1..ke only
        ! front(k+1:nfront, k+1:ke) -= L(:,k) * U(k, k+1:ke)
        nupd = ke - k   ! number of remaining columns in panel
        if (nupd > 0 .and. nel > 0) then
          call dger(nel, nupd, ALPHA, front(k+1, k), 1, &
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
        call dtrsm('L', 'L', 'N', 'U', nb, nel, ONE, &
                   front(kb, kb), nfront, front(kb, ke+1), nfront)

        ! Step 2: Schur complement update
        ! front(ke+1:nfront, ke+1:nfront) -= L(ke+1:nfront, kb:ke) * U(kb:ke, ke+1:nfront)
        call dgemm('N', 'N', nel, nel, nb, ALPHA, &
                   front(ke+1, kb), nfront, &
                   front(kb, ke+1), nfront, &
                   ONE, front(ke+1, ke+1), nfront)
      end if

      kb = ke + 1
    end do
  end subroutine

  ! Multifrontal solve: L U x = P b
  !
  ! Phase 1: Forward substitution  (L y = P b)  — bottom-up through tree
  ! Phase 2: Backward substitution (U x = y)    — top-down through tree
  subroutine multifrontal_solve(mat, lu, rhs, x)
    type(matrix_data), intent(in) :: mat
    type(monolis_mat_lu), intent(in) :: lu
    double precision, intent(in)  :: rhs(mat%NDOF*mat%N)
    double precision, intent(out) :: x(mat%NDOF*mat%N)

    integer :: n, ndof, ntot, nsuper
    integer, allocatable :: postorder(:)
    double precision, allocatable :: work(:)
    integer :: ip, s, k, i, j, nfront, npiv
    integer :: kblk, kdof, iblk, idof, gi, gj
    integer :: nf_dof, np_dof
    double precision :: sum_val

    n = mat%N; ndof = mat%NDOF; nsuper = lu%NSUPER
    ntot = n * ndof

    associate(sstart => lu%SNODE_START, ssize => lu%SNODE_SIZE, &
              sfsize => lu%SNODE_FSIZE, sparent => lu%SNODE_PARENT, &
              sfils => lu%SFILS, sfrere => lu%SFRERE, &
              perm => mat%PERM, invperm => mat%INVPERM, &
              factors => lu%FACTORS)

    allocate(work(ntot), postorder(nsuper))

    ! Apply block permutation to RHS
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

  ! Build postorder array (simpler version for solve phase)
  subroutine build_postorder_array(nsuper, sfils, sfrere, sparent, postorder)
    integer, intent(in)  :: nsuper
    integer, intent(in)  :: sfils(nsuper), sfrere(nsuper), sparent(nsuper)
    integer, intent(out) :: postorder(nsuper)

    integer :: stack(nsuper)
    integer :: sp, s, cnt, child
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

end module mod_monolis_fact_factorize
