!> 直接法向けフィルイン判定モジュール
module mod_monolis_spmat_fillin
  use mod_monolis_utils
  use mod_monolis_def_mat

  implicit none

  !> フィルイン構造体
  type monolis_fillin
    !> 祖先の数
    integer(kint) :: n_ancestor
    !> 祖先のインデックス
    integer(kint), pointer :: ancestor(:)
  endtype monolis_fillin

contains

  !> フィルインの決定（シンボリック解析）
  subroutine monolis_matrix_get_fillin(CSR, SCSR, is_asym)
    implicit none
    !>
    type(monolis_mat_CSR) :: CSR
    !>
    type(monolis_mat_separated_CSR) :: SCSR
    !>
    logical :: is_asym
    type(monolis_fillin), allocatable :: tree(:)
    integer(kint), pointer :: array(:)
    integer(kint), allocatable :: fillin_mask(:)
    integer(kint), allocatable :: child_mask(:)
    integer(kint), allocatable :: parent_mask(:)
    integer(kint), allocatable :: count(:)
    integer(kint) :: N, NPU, NPL
    integer(kint) :: i, j, k, jS, jE, in, c
    integer(kint) :: nbytes
    integer(kint) :: is, ie
    integer(kint) :: range, parent
    integer(kint) :: bit = kint*8

    N = CSR%N
    allocate(tree(N))

    do i = 1, N
      tree(i)%n_ancestor = 0

      in = 0
      jS = CSR%index(i-1) + 1
      jE = CSR%index(i  )
      do j = jS, jE
        if(i < CSR%item(j) .and. CSR%item(j) <= N)then
          in = in + 1
        endif
      enddo

      tree(i)%n_ancestor = in
      allocate(tree(i)%ancestor(in))

      in = 0
      do j = jS, jE
        if(i < CSR%item(j) .and. CSR%item(j) <= N)then
          in = in + 1
          tree(i)%ancestor(in) = CSR%item(j)
        endif
      enddo
    enddo

    nbytes = N/bit+1

    call monolis_alloc_I_1d(child_mask, nbytes)
    call monolis_alloc_I_1d(parent_mask, nbytes)
    call monolis_alloc_I_1d(fillin_mask, nbytes)

    do i = 1, N
      if(tree(i)%n_ancestor < 2) cycle
      is = i/bit + 1
      child_mask(is:nbytes) = 0
      parent_mask(is:nbytes) = 0

      parent = tree(i)%ancestor(1)
      range = 0
      do j = 2, tree(i)%n_ancestor
        in = tree(i)%ancestor(j)
        ie = in/bit + 1
        child_mask(ie) = ibset(child_mask(ie), mod(in,bit))
        range = in
      enddo

      k = tree(parent)%n_ancestor
      do j = 1, k
        in = tree(parent)%ancestor(j)
        ie = in/bit + 1
        parent_mask(ie) = ibset(parent_mask(ie), mod(in,bit))
        range = max(range,in)
      enddo

      ie = range/bit + 1
      fillin_mask(is:ie) = ior(child_mask(is:ie), parent_mask(is:ie))

      c = 0
      do j = is, ie
        c = c + popcnt(fillin_mask(j))
      enddo

      if(0 < c)then
        allocate(array(c))
        tree(parent)%n_ancestor = c
        in = 0
        do j = is, ie
          do k = 1, popcnt(fillin_mask(j))
            in = in + 1
            c = popcnt( iand(fillin_mask(j), - fillin_mask(j)) -1 )
            fillin_mask(j) = ibclr(fillin_mask(j),c)
            array(in) = bit*(j-1)+c
          enddo
        enddo
        deallocate(tree(parent)%ancestor)
        tree(parent)%ancestor => array
      endif
    enddo

    call monolis_dealloc_I_1d(child_mask )
    call monolis_dealloc_I_1d(parent_mask)
    call monolis_dealloc_I_1d(fillin_mask)

    !> upper part
    call monolis_alloc_I_1d(SCSR%indexU, N + 1)

    in = 0
    SCSR%indexU(0) = 0
    do i = 1, N
      SCSR%indexU(i) = SCSR%indexU(i-1) + tree(i)%n_ancestor
      in = in + tree(i)%n_ancestor
    enddo

    NPU = in
    SCSR%NPU = in

    call monolis_alloc_I_1d(SCSR%itemU, NPU)

    in = 0
    do i = 1, N
      do j = 1, tree(i)%n_ancestor
        in = in + 1
        SCSR%itemU(in) = tree(i)%ancestor(j)
      enddo
    enddo

    if(is_asym)then
      !> lower part
      call monolis_alloc_I_1d(count, N)
      call monolis_alloc_I_1d(SCSR%indexL, N + 1)
      count = 0
      SCSR%indexL(0) = 0

      do i = 1, SCSR%indexU(N)
        in = SCSR%itemU(i)
        count(in) = count(in) + 1
      enddo

      do i = 1, N
        SCSR%indexL(i) = SCSR%indexL(i-1) + count(i)
      enddo

      NPL = NPU
      c = 1
      call monolis_alloc_I_1d(SCSR%itemL, NPL)
      do i = 1, N
        aa:do k = 1, i
          jS = SCSR%indexU(k-1) + 1
          jE = SCSR%indexU(k)
          do j = jS, jE
            in = SCSR%itemU(j)
            if(i < in) cycle aa
            if(i == in)then
              SCSR%itemL(c) = k
              c = c + 1
              cycle aa
            endif
          enddo
        enddo aa
      enddo
      call monolis_dealloc_I_1d(count)
    endif
  end subroutine monolis_matrix_get_fillin

!  subroutine monolis_matrix_copy_with_fillin(monoPRM, monoCOM, monoMAT, monoTREE, is_asym)
!    implicit none
!    type(monolis_prm) :: monoPRM
!    type(monolis_com) :: monoCOM
!    type(monolis_mat) :: monoMAT
!    type(monolis_mat_LDU) :: monoTREE
!    integer(kint) :: N, NPU, NPL, NDOF, NDOF2
!    integer(kint) :: i, j, k, l, iS, iE, jS, jE, lS, lE
!    integer(kint) :: in, jn
!    logical :: is_asym
!
!    N = monoTREE%N
!    NPU = monoTREE%NPU
!    NDOF = monoTREE%NDOF
!    NDOF2 = NDOF*NDOF
!
!    !value
!    allocate(monoTREE%D(NDOF2*N), source = 0.0d0)
!    do i = 1, N
!      jS = monoMAT%index(i-1) + 1
!      jE = monoMAT%index(i)
!      do j = jS, jE
!        in = monoMAT%item(j)
!        if(i == in)then
!          do l = 1, NDOF2
!            monoTREE%D(NDOF2*(i-1) + l) = monoMAT%A(NDOF2*(j-1) + l)
!          enddo
!        endif
!      enddo
!    enddo
!
!    allocate(monoTREE%U(NDOF2*NPU), source = 0.0d0)
!
!    do k = 1, N
!      iS = monoTREE%indexU(k-1) + 1
!      iE = monoTREE%indexU(k)
!      jS = monoMAT%index(k-1) + 1
!      jE = monoMAT%index(k)
!      aa:do j = jS, jE
!        jn = monoMAT%item(j)
!        if(k < jn)then
!          do i = iS, iE
!            in = monoTREE%itemU(i)
!            if(jn == in)then
!              lS = NDOF2*(i-1)
!              lE = NDOF2*(j-1)
!              do l = 1, NDOF2
!                monoTREE%U(lS + l) = monoMAT%A(lE + l)
!              enddo
!              iS = iS + 1
!              cycle aa
!            endif
!          enddo
!        endif
!      enddo aa
!    enddo
!
!    if(is_asym)then
!      monoTREE%NPL = NPU
!      NPL = NPU
!      allocate(monoTREE%L(NDOF2*NPL), source = 0.0d0)
!
!      do k = 1, N
!        iS = monoTREE%indexL(k-1) + 1
!        iE = monoTREE%indexL(k)
!        jS = monoMAT%index(k-1) + 1
!        jE = monoMAT%index(k)
!        bb:do j = jS, jE
!          jn = monoMAT%item(j)
!          if(jn < k)then
!            do i = iS, iE
!              in = monoTREE%itemL(i)
!              if(jn == in)then
!                lS = NDOF2*(i-1)
!                lE = NDOF2*(j-1)
!                do l = 1, NDOF2
!                  monoTREE%L(lS + l) = monoMAT%A(lE + l)
!                enddo
!                iS = iS + 1
!                cycle bb
!              endif
!            enddo
!          endif
!        enddo bb
!      enddo
!    endif
!  end subroutine monolis_matrix_copy_with_fillin
end module mod_monolis_spmat_fillin