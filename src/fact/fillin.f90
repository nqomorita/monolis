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
    type(monolis_fillin), allocatable:: tree(:)
    integer(kint), pointer :: array(:)
    integer(kint), allocatable :: fillin_mask(:)
    integer(kint), allocatable :: child_mask(:)
    integer(kint), allocatable :: parent_mask(:)
    integer(kint), allocatable :: parent_id(:)
    integer(kint), allocatable :: perm(:)
    integer(kint), allocatable :: temp(:)
    integer(kint), allocatable :: index(:)
    integer(kint), allocatable :: is_used(:)
    integer(kint) :: N, M, NPU
    integer(kint) :: i, j, k, p, jS, jE, in, c, d
    integer(kint) :: nbytes
    integer(kint) :: is, ie
    integer(kint) :: range, parent, child
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

    !> relax supernode part
    if(.true.)then
      allocate(is_used(N), source = 0)
      allocate(temp(N), source = 0)
      allocate(perm(N), source = 0)

      do i = 1, N - 1
        temp(i) = tree(i)%ancestor(1)
        perm(i) = (i)
      enddo
      temp(N) = N + 1
      perm(N) = N + 1

      call monolis_qsort_I_2d(temp, perm, 1, N)

      allocate(parent_id(N), source = 0)
      parent_id = temp
      call monolis_get_uniq_array_I(parent_id, N, M)

      allocate(index(M + 1), source = 0)
      in = 2
      do i = 1, N - 1
        if(temp(i) == temp(i + 1))then
          index(in) = index(in) + 1
        else 
          index(in) = index(in) + 1
          in = in + 1
        endif
      enddo

      do i = 2, M + 1
        index(i) = index(i) + index(i-1)
      enddo

      do i = 1, M - 1
        jS = index(i) + 1
        jE = index(i + 1)
        call monolis_qsort_I_1d(perm, jS, jE)
      enddo

write(*,*)"temp"
write(*,"(20i4)")temp
write(*,*)"perm"
write(*,"(20i4)")perm
write(*,*)"index"
write(*,"(20i4)")index

      nbytes = N/bit+1
      allocate(child_mask (nbytes), source = 0)
      allocate(parent_mask(nbytes), source = 0)
      allocate(fillin_mask(nbytes), source = 0)

      do i = M - 1, 1, -1
        jS = index(i) + 1
        jE = index(i + 1)

        is = minval(perm(jS:jE))/bit + 1
        parent_mask(is:nbytes) = 0

        range = 0
        parent = parent_id(i)
write(*,*)"parent", parent
        call set_bit(parent, parent_mask, range, bit)
        do j = 1, tree(parent)%n_ancestor
          in = tree(parent)%ancestor(j)
          call set_bit(in, parent_mask, range, bit)
        enddo

        do m = jE, jS, -1
          child = perm(m)
          if(child == N - tree(child)%n_ancestor) cycle
write(*,*)"child", child
          child_mask(is:nbytes) = 0

          do j = 1, tree(child)%n_ancestor
            in = tree(child)%ancestor(j)
            call set_bit(in, child_mask, range, bit)
          enddo

          ie = range/bit + 1
          fillin_mask(is:ie) = ieor(child_mask(is:ie), parent_mask(is:ie))

          c = 0
          do j = is, ie
            c = c + popcnt(fillin_mask(j))
          enddo
          write(*,*)"C", c

          if(c == 0) cycle
          call set_bit(child, child_mask, range, bit)
          fillin_mask(is:ie) = ior(child_mask(is:ie), parent_mask(is:ie))
          
          allocate(array(c))
          tree(child)%n_ancestor = c
          in = 0
          do j = is, ie
            do k = 1, popcnt(fillin_mask(j))
              in = in + 1
              c = popcnt( iand(fillin_mask(j), - fillin_mask(j)) -1 )
              fillin_mask(j) = ibclr(fillin_mask(j),c)
              array(in) = bit*(j-1)+c
            enddo
          enddo
          deallocate(tree(child)%ancestor)
          tree(parent)%ancestor => array
        enddo
      enddo

      deallocate(child_mask )
      deallocate(parent_mask)
      deallocate(fillin_mask)
    endif

    !> SCSR allocation part
    monoTREE%N = monoMAT%N
    monoTREE%NDOF = monoMAT%NDOF

    allocate(monoTREE%SCSR%indexU(N+1))
    NPU = 0
    monoTREE%SCSR%indexU(1) = 0
    do i = 1, N
      monoTREE%SCSR%indexU(i+1) = monoTREE%SCSR%indexU(i) + tree(i)%n_ancestor + 1
      NPU = NPU + tree(i)%n_ancestor + 1
    enddo

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

  subroutine set_bit(in, parent_mask, range, bit)
    implicit none
    integer(kint) :: in, range, ie
    integer(kint) :: parent_mask(:)
    integer(kint) :: bit 
    ie = in/bit + 1
    parent_mask(ie) = ibset(parent_mask(ie), mod(in, bit))
    range = max(range, in)
  end subroutine set_bit

  subroutine monolis_matrix_alloc_with_fillin(monoTREE, is_asym)
    implicit none
    type(monolis_mat) :: monoTREE
    logical :: is_asym
    integer(kint) :: N, NDOF, NZ

    NDOF = monoTREE%NDOF
    N = monoTREE%N
    NZ = monoTREE%SCSR%indexU(N + 1)

    call monolis_palloc_R_1d(monoTREE%R%A, NZ*NDOF*NDOF)
  end subroutine monolis_matrix_alloc_with_fillin
end module mod_monolis_fact_fillin
