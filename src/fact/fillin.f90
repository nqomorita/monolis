module mod_monolis_fact_fillin
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc

  implicit none

  type monolis_fillin
    integer(kint) :: n_ancestor
    integer(kint), pointer :: ancestor(:)
  end type monolis_fillin

contains

  subroutine monolis_matrix_get_fillin(monoMAT, monoTREE, is_asym, is_fillin)
    implicit none
    type(monolis_mat) :: monoMAT
    type(monolis_mat) :: monoTREE
    type(monolis_fillin), pointer:: tree(:)
    integer(kint), pointer :: array(:)
    integer(kint), pointer :: fillin_mask(:)
    integer(kint), pointer :: child_mask(:)
    integer(kint), pointer :: parent_mask(:)
    integer(kint) :: N, NPU, NPL
    integer(kint) :: i, j, k, jS, jE, in, c
    integer(kint) :: nbytes
    integer(kint) :: is, ie
    integer(kint) :: range, parent
    integer(kint) :: bit = kint*8
    integer(kint), allocatable :: count(:)
    logical :: is_asym, is_fillin

    N = monoMAT%N
    allocate(tree(N))

    do i = 1, N
      tree(i)%n_ancestor = 0

      in = 0
      jS = monoMAT%CSR%index(i) + 1
      jE = monoMAT%CSR%index(i + 1)
      do j = jS, jE
        if(i < monoMAT%CSR%item(j) .and. monoMAT%CSR%item(j) <= N)then
          in = in + 1
        endif
      enddo

      tree(i)%n_ancestor = in
      allocate(tree(i)%ancestor(in))

      in = 0
      do j = jS, jE
        if(i < monoMAT%CSR%item(j) .and. monoMAT%CSR%item(j) <= N)then
          in = in + 1
          tree(i)%ancestor(in) = monoMAT%CSR%item(j)
        endif
      enddo
    enddo

    if(is_fillin)then
      nbytes = N/bit+1
      allocate(child_mask (nbytes), source = 0)
      allocate(parent_mask(nbytes), source = 0)
      allocate(fillin_mask(nbytes), source = 0)

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

      deallocate(child_mask )
      deallocate(parent_mask)
      deallocate(fillin_mask)
    endif

    !> SCSR allocation part
    monoTREE%N = monoMAT%N
    monoTREE%NDOF = monoMAT%NDOF

    allocate(monoTREE%SCSR%indexU(N+1))
    in = 0
    monoTREE%SCSR%indexU(1) = 0
    do i = 1, N
      monoTREE%SCSR%indexU(i+1) = monoTREE%SCSR%indexU(i) + tree(i)%n_ancestor + 1
      in = in + tree(i)%n_ancestor + 1
    enddo

    NPU = in

    allocate(monoTREE%SCSR%itemU(NPU))
    in = 0
    do i = 1, N
      in = in + 1
      monoTREE%SCSR%itemU(in) = i
      do j = 1, tree(i)%n_ancestor
        in = in + 1
        monoTREE%SCSR%itemU(in) = tree(i)%ancestor(j)
      enddo
    enddo

!    if(is_asym)then
!      !> lower part
!      allocate(count(N))
!      allocate(monoTREE%SCSR%indexL(0:N))
!      count = 0
!      monoTREE%SCSR%indexL(0) = 0
!
!      do i = 1, monoTREE%SCSR%indexU(N)
!        in = monoTREE%SCSR%itemU(i)
!        count(in) = count(in) + 1
!      enddo
!
!      do i = 1, N
!        monoTREE%SCSR%indexL(i) = monoTREE%SCSR%indexL(i-1) + count(i)
!      enddo
!
!      NPL = NPU
!      c = 1
!      allocate(monoTREE%SCSR%itemL(NPL))
!      do i = 1, N
!        aa:do k = 1, i
!          jS = monoTREE%SCSR%indexU(k-1) + 1
!          jE = monoTREE%SCSR%indexU(k)
!          do j = jS, jE
!            in = monoTREE%SCSR%itemU(j)
!            if(i < in) cycle aa
!            if(i == in)then
!              monoTREE%SCSR%itemL(c) = k
!              c = c + 1
!              cycle aa
!            endif
!          enddo
!        enddo aa
!      enddo
!      deallocate(count)
!    endif
  end subroutine monolis_matrix_get_fillin

  subroutine monolis_matrix_alloc_with_fillin(monoTREE, is_asym)
    implicit none
    type(monolis_mat) :: monoTREE
    logical :: is_asym
    integer(kint) :: N, NDOF, NZ

    NDOF = monoTREE%NDOF
    N = monoTREE%N
    NZ = monoTREE%SCSR%indexU(N + 1)

    call monolis_palloc_R_1d(monoTREE%R%A, NZ*NDOF*NDOF)
    !call monolis_palloc_R_1d(MAT%R%B, NP*NDOF)
    !call monolis_palloc_R_1d(MAT%R%X, NP*NDOF)
  end subroutine monolis_matrix_alloc_with_fillin

!  subroutine monolis_matrix_copy_with_fillin(monoMAT, monoTREE, is_asym)
!    implicit none
!    type(monolis_mat) :: monoMAT
!    type(monolis_mat) :: monoTREE
!    integer(kint) :: N, NPU, NPL, NDOF, NDOF2
!    integer(kint) :: i, j, k, l, iS, iE, jS, jE, lS, lE
!    integer(kint) :: in, jn
!    logical :: is_asym
!
!    N = monoTREE%N
!    NPU = monoTREE%SCSR%indexU(N + 1)
!    NDOF = monoTREE%NDOF
!    NDOF2 = NDOF*NDOF
!
!    !value
!    allocate(monoTREE%R%D(NDOF2*N), source = 0.0d0)
!    do i = 1, N
!      jS = monoMAT%CSR%index(i-1) + 1
!      jE = monoMAT%CSR%index(i)
!      do j = jS, jE
!        in = monoMAT%CSR%item(j)
!        if(i == in)then
!          do l = 1, NDOF2
!            monoTREE%R%D(NDOF2*(i-1) + l) = monoMAT%R%A(NDOF2*(j-1) + l)
!          enddo
!        endif
!      enddo
!    enddo
!
!    allocate(monoTREE%R%U(NDOF2*NPU), source = 0.0d0)
!
!    do k = 1, N
!      iS = monoTREE%SCSR%indexU(k-1) + 1
!      iE = monoTREE%SCSR%indexU(k)
!      jS = monoMAT%CSR%index(k-1) + 1
!      jE = monoMAT%CSR%index(k)
!      aa:do j = jS, jE
!        jn = monoMAT%CSR%item(j)
!        !if(k < jn)then
!        if(k <= jn)then
!          do i = iS, iE
!            in = monoTREE%SCSR%itemU(i)
!            if(jn == in)then
!              lS = NDOF2*(i-1)
!              lE = NDOF2*(j-1)
!              do l = 1, NDOF2
!                monoTREE%R%U(lS + l) = monoMAT%R%A(lE + l)
!              enddo
!              iS = iS + 1
!              cycle aa
!            endif
!          enddo
!        endif
!      enddo aa
!    enddo
!
!!    if(is_asym)then
!!      monoTREE%NPL = NPU
!!      NPL = NPU
!!      allocate(monoTREE%L(NDOF2*NPL), source = 0.0d0)
!!
!!      do k = 1, N
!!        iS = monoTREE%indexL(k-1) + 1
!!        iE = monoTREE%indexL(k)
!!        jS = monoMAT%index(k-1) + 1
!!        jE = monoMAT%index(k)
!!        bb:do j = jS, jE
!!          jn = monoMAT%item(j)
!!          if(jn < k)then
!!            do i = iS, iE
!!              in = monoTREE%itemL(i)
!!              if(jn == in)then
!!                lS = NDOF2*(i-1)
!!                lE = NDOF2*(j-1)
!!                do l = 1, NDOF2
!!                  monoTREE%L(lS + l) = monoMAT%A(lE + l)
!!                enddo
!!                iS = iS + 1
!!                cycle bb
!!              endif
!!            enddo
!!          endif
!!        enddo bb
!!      enddo
!!    endif
!  end subroutine monolis_matrix_copy_with_fillin
end module mod_monolis_fact_fillin
