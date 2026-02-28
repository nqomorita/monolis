module mod_monolis_fact_factorize
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc

  implicit none

contains

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
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: mat
    type(monolis_mat), intent(inout) :: lu

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

    n = mat%N
    nz = mat%CSR%index(mat%NP + 1)
    ndof = mat%NDOF
    nsuper = lu%lu%nsuper
    info = 0

    allocate(lu%lu%factors(nsuper))

    associate(row_ptr => mat%CSR%index, col_ind => mat%CSR%item, &
              a_elt => lu%R%A, &
              invperm => mat%REORDER%iperm, perm => mat%REORDER%perm, &
              sstart => lu%lu%snode_start, ssize => lu%lu%snode_size, &
              sparent => lu%lu%snode_parent, sfsize => lu%lu%snode_fsize, &
              sfils => lu%lu%sfils, sfrere => lu%lu%sfrere, factors => lu%lu%factors)

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
      if (nfront - npiv > 1)then
        call monolis_qsort_I_1d(idx_work, npiv + 1, nfront)
      endif

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
  ! Multifrontal solve: L U x = P b
  !
  ! Phase 1: Forward substitution  (L y = P b)  — bottom-up through tree
  ! Phase 2: Backward substitution (U x = y)    — top-down through tree
  !==============================================================================
  subroutine multifrontal_solve(mat, lu, rhs, x)
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: mat
    type(monolis_mat), intent(in) :: lu
    real(kdouble), intent(in)  :: rhs(:)
    real(kdouble), intent(out) :: x(:)

    integer(kint) :: n, ndof, ntot, nsuper
    integer(kint), allocatable :: postorder(:)
    real(kdouble), allocatable :: work(:)
    integer(kint) :: ip, s, k, i, j, nfront, npiv
    integer(kint) :: kblk, kdof, iblk, idof, gi, gj
    integer(kint) :: nf_dof, np_dof
    real(kdouble) :: sum_val

    n = mat%N
    ndof = mat%NDOF
    nsuper = lu%lu%nsuper
    ntot = n * ndof

    associate(sstart => lu%lu%snode_start, ssize => lu%lu%snode_size, &
              sfsize => lu%lu%snode_fsize, sparent => lu%lu%snode_parent, &
              sfils => lu%lu%sfils, sfrere => lu%lu%sfrere, &
              invperm => mat%REORDER%iperm, perm => mat%REORDER%perm, &
              factors => lu%lu%factors)

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

end module mod_monolis_fact_factorize
