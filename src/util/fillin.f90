module mod_monolis_matrix_fillin
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat

  implicit none

  type monolis_fillin
    integer(kind=kint) :: n_descendant
    integer(kind=kint) :: n_ancestor
    integer(kind=kint), pointer :: descendant(:)
    integer(kind=kint), pointer :: ancestor(:)
    integer(kind=kint), pointer :: update_index(:)
    real(kind=kdouble), pointer :: update(:,:,:)
    logical :: factorized
    logical :: updated
  endtype monolis_fillin

contains

  subroutine monolis_matrix_get_fillin(monoPRM, monoCOM, monoMAT, monoTREE, is_fillin, is_asym)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    type(monolis_mat) :: monoTREE
    type(monolis_fillin), pointer:: T(:)
    integer(kind=kint), pointer :: idxU(:)
    integer(kind=kint), pointer :: itemU(:)
    integer(kind=kint), pointer :: idxL(:)
    integer(kind=kint), pointer :: itemL(:)
    integer(kind=kint) :: N, NPU, NPL
    integer(kind=kint) :: i, j, k, jS, jE, in, c
    integer(kind=kint) :: Nbytes
    integer(kind=kint) :: Start, End
    integer(kind=kint) :: range, parent
    integer(kind=kint) :: ZERO = 0
    integer(kind=kint), pointer :: array(:)
    integer(kind=kint), pointer :: fillin_mask(:)
    integer(kind=kint), pointer :: child_mask(:)
    integer(kind=kint), pointer :: parent_mask(:)
    integer(kind=kint) :: bit = kint*8
    integer(kind=kint), allocatable :: count(:), diff(:)
    logical :: is_fillin, is_asym

    N = monoMAT%N
    allocate(T(N))

    do i = 1, N
      T(i)%factorized = .false.
      T(i)%updated = .false.
      T(i)%n_descendant = 0
      T(i)%n_ancestor = 0

      jS = monoMAT%indexU(i-1) + 1
      jE = monoMAT%indexU(i  )
      in = 0
      do j = jS, jE
        if(monoMAT%itemU(j) <= N)then
          in = in + 1
        endif
      enddo
      T(i)%n_ancestor = in
      allocate(T(i)%ancestor(in))

      in = 0
      do j = jS, jE
        if(monoMAT%itemU(j) <= N)then
          in = in + 1
          T(i)%ancestor(in) = monoMAT%itemU(j)
        endif
      enddo
    enddo

    if(is_fillin)then
      Nbytes = N/bit+1
      allocate(child_mask (Nbytes))
      allocate(parent_mask(Nbytes))
      allocate(fillin_mask(Nbytes))

      do i = 1, N
        if(T(i)%n_ancestor < 2) cycle
        Start = i/bit + 1
        child_mask(Start:Nbytes) = 0
        parent_mask(Start:Nbytes) = 0

        parent = T(i)%ancestor(1)
        range = 0
        do j = 2, T(i)%n_ancestor
          in = T(i)%ancestor(j)
          End = in/bit + 1
          child_mask(End) = ibset(child_mask(End),mod(in,bit))
          range = in
        enddo
        k = T(parent)%n_ancestor
        do j = 1, k
          in = T(parent)%ancestor(j)
          End = in/bit + 1
          parent_mask(End) = ibset(parent_mask(End),mod(in,bit))
          range = max(range,in)
        enddo
        End = range/bit + 1

        fillin_mask(Start:End) = ior(child_mask(Start:End), parent_mask(Start:End))

        c = 0
        do j = Start, End
          c = c + popcnt(fillin_mask(j))
        enddo

        if(0 < c)then
          allocate(array(c))
          T(parent)%n_ancestor=c
          in = 0
          do j = Start, End
            do k = 1, popcnt(fillin_mask(j))
              in = in + 1
              c = popcnt( iand(fillin_mask(j), - fillin_mask(j)) -1 )
              fillin_mask(j) = ibclr(fillin_mask(j),c)
              array(in) = bit*(j-1)+c
            enddo
          enddo
          deallocate(T(parent)%ancestor)
          T(parent)%ancestor => array
        endif
      enddo
      deallocate(child_mask )
      deallocate(parent_mask)
      deallocate(fillin_mask)

      !do i = 1, N
      !  do j = 1, T(i)%n_ancestor
      !    T(T(i)%ancestor(j))%n_descendant = T(T(i)%ancestor(j))%n_descendant + 1
      !  enddo
      !enddo
      !do i = 1, N
      !  allocate(T(i)%descendant(T(i)%n_descendant))
      !enddo

      !allocate(diff(N))
      !diff(:) = 0
      !do i = 1, N
      !  do j = 1, T(i)%n_ancestor
      !    diff(T(i)%ancestor(j)) = diff(T(i)%ancestor(j)) + 1
      !    T(T(i)%ancestor(j))%descendant(diff(T(i)%ancestor(j))) = i
      !  enddo
      !enddo
    endif

    idxU => monoTREE%indexU
    itemU => monoTREE%itemU

    allocate(idxU(0:N))
    in = 0
    idxU(0) = 0
    do i = 1, N
      idxU(i) = idxU(i-1) + T(i)%n_ancestor + 1
      in = in + T(i)%n_ancestor + 1
    enddo
    NPU = in

    allocate(itemU(NPU))
    in = 0
    do i = 1, N
      in = in + 1
      itemU(in) = i
      do j = 1, T(i)%n_ancestor
        in = in + 1
        itemU(in) = T(i)%ancestor(j)
      enddo
    enddo

    if(is_asym)then
      !lower part
      idxL => monoTREE%indexL
      itemL => monoTREE%itemL

      allocate(count(N))
      allocate(idxL(0:N))
      count = 0
      idxL(0) = 0

      do i = 1, idxU(N)
        in = itemU(i)
        count(in) = count(in) + 1
      enddo

      do i = 1, N
        idxL(i) = idxL(i-1) + count(i)
      enddo

      NPL = NPU
      c = 1
      allocate(itemL(NPL))
      do i = 1, N
        aa:do k = 1, i
          jS = idxU(k-1) + 1
          jE = idxU(k)
          do j = jS, jE
            in = itemU(j)
            if(i < in) cycle aa
            if(i == in)then
              itemL(c) = k
              c = c + 1
              cycle aa
            endif
          enddo
        enddo aa
      enddo
    endif
  end subroutine monolis_matrix_get_fillin

  subroutine monolis_matrix_copy_with_fillin(monoPRM, monoCOM, monoMAT, monoTREE, is_asym)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    type(monolis_mat) :: monoTREE
    integer(kind=kint), pointer :: idxU(:)
    integer(kind=kint), pointer :: itemU(:)
    integer(kind=kint), pointer :: idxL(:)
    integer(kind=kint), pointer :: itemL(:)
    real(kind=kdouble), pointer :: AU(:)
    real(kind=kdouble), pointer :: AL(:)
    integer(kind=kint) :: N, NPU, NPL
    integer(kind=kint) :: i, j, k, iS, iE, jS, jE
    integer(kind=kint) :: in, jn ,kn, nn
    logical :: is_asym

    N = monoTREE%N
    NPU = monoTREE%NPU
    idxU => monoTREE%indexU
    itemU => monoTREE%itemU
    idxL => monoTREE%indexL
    itemL => monoTREE%itemL

    !value
    allocate(monoTREE%AU(9*NPU))
    AU => monoTREE%AU
    AU = 0.0d0

    do k = 1, N
      in = idxU(k-1)+1
      AU(9*in-8) = monoMAT%D(9*k-8)
      AU(9*in-7) = monoMAT%D(9*k-7)
      AU(9*in-6) = monoMAT%D(9*k-6)
      AU(9*in-5) = monoMAT%D(9*k-5)
      AU(9*in-4) = monoMAT%D(9*k-4)
      AU(9*in-3) = monoMAT%D(9*k-3)
      AU(9*in-2) = monoMAT%D(9*k-2)
      AU(9*in-1) = monoMAT%D(9*k-1)
      AU(9*in  ) = monoMAT%D(9*k  )
    enddo

    do k = 1, N
      iS = idxU(k-1) + 1
      iE = idxU(k)
      jS = monoMAT%indexU(k-1) + 1
      jE = monoMAT%indexU(k)
      aa:do j = jS, jE
        jn = monoMAT%itemU(j)
        do i = iS, iE
          in = itemU(i)
          if(jn == in)then
            AU(9*i-8) = monoMAT%AU(9*j-8)
            AU(9*i-7) = monoMAT%AU(9*j-7)
            AU(9*i-6) = monoMAT%AU(9*j-6)
            AU(9*i-5) = monoMAT%AU(9*j-5)
            AU(9*i-4) = monoMAT%AU(9*j-4)
            AU(9*i-3) = monoMAT%AU(9*j-3)
            AU(9*i-2) = monoMAT%AU(9*j-2)
            AU(9*i-1) = monoMAT%AU(9*j-1)
            AU(9*i  ) = monoMAT%AU(9*j  )
            iS = iS + 1
            cycle aa
          endif
        enddo
      enddo aa
    enddo

    if(is_asym)then
      allocate(AL(9*NPL))
      AL = 0.0d0

      do k = 1, N
        iS = idxL(k-1) + 1
        iE = idxL(k)
        jS = monoMAT%indexL(k-1) + 1
        jE = monoMAT%indexL(k)
        bb:do j = jS, jE
          jn = monoMAT%itemL(j)
          do i = iS, iE
            in = itemL(i)
            if(jn == in)then
              AL(9*i-8) = monoMAT%AL(9*j-8)
              AL(9*i-7) = monoMAT%AL(9*j-7)
              AL(9*i-6) = monoMAT%AL(9*j-6)
              AL(9*i-5) = monoMAT%AL(9*j-5)
              AL(9*i-4) = monoMAT%AL(9*j-4)
              AL(9*i-3) = monoMAT%AL(9*j-3)
              AL(9*i-2) = monoMAT%AL(9*j-2)
              AL(9*i-1) = monoMAT%AL(9*j-1)
              AL(9*i  ) = monoMAT%AL(9*j  )
              iS = iS + 1
              cycle bb
            endif
          enddo
        enddo bb
      enddo
    endif
  end subroutine monolis_matrix_copy_with_fillin
end module mod_monolis_matrix_fillin