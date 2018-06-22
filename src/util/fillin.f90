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

  subroutine monolis_matrix_get_fillin(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    type(monolis_fillin), pointer:: T(:)
    integer(kind=kint), pointer :: idxU(:)
    integer(kind=kint), pointer :: itemU(:)
    integer(kind=kint) :: N, NPU
    integer(kind=kint) :: i, j, k, jS, jE, in, c
    integer(kind=kint) :: Nbytes
    integer(kind=kint) :: Start, End
    integer(kind=kint) :: range, parent
    integer(kind=kint) :: ZERO = 0
    integer(kind=kint), pointer :: Array(:)
    integer(kind=kint), pointer :: FillinMask(:)
    integer(kind=kint), pointer :: ChildMask(:)
    integer(kind=kint), pointer :: ParentMask(:)
    integer(kind=kint) :: bit = kint*8
    integer(kind=kint), allocatable :: count(:), diff(:)
    logical :: is_asym = .false.

    N = monoMAT%N
    allocate(T(N))

    do i = 1, N
      T(i)%factorized = .false.
      T(i)%updated = .false.
      T(i)%n_descendant = 0
      T(i)%n_ancestor = 0
      jS = monoMAT%indexU(i-1) + 1
      jE = monoMAT%indexU(i  )
      T(i)%n_ancestor = jE - jS + 1
      allocate(T(i)%ancestor(T(i)%n_ancestor))
      in = 0
      do j = jS, jE
        in = in + 1
        T(i)%ancestor(in) = monoMAT%itemU(j)
      enddo
    enddo

    Nbytes = N/bit+1
    allocate(ChildMask (Nbytes))
    allocate(ParentMask(Nbytes))
    allocate(FillinMask(Nbytes))

    do i = 1, N
      if(T(i)%n_ancestor < 2) cycle
      Start = i/bit + 1
      ChildMask(Start:Nbytes) = 0
      ParentMask(Start:Nbytes) = 0

      parent = T(i)%ancestor(1)
      range = 0
      do j = 2, T(i)%n_ancestor
        in = T(i)%ancestor(j)
        End = in/bit + 1
        ChildMask(End) = ibset(ChildMask(End),mod(in,bit))
        range = in
      enddo
      k = T(parent)%n_ancestor
      do j = 1, k
        in = T(parent)%ancestor(j)
        End = in/bit + 1
        ParentMask(End) = ibset(ParentMask(End),mod(in,bit))
        range = max(range,in)
      enddo
      End = range/bit + 1

      FillinMask(Start:End) = ior(ChildMask(Start:End), ParentMask(Start:End))

      c = 0
      do j = Start, End
        c = c + popcnt(FillinMask(j))
      enddo

      if(0 < c)then
        allocate(Array(c))
        T(parent)%n_ancestor=c
        in = 0
        do j = Start, End
          do k = 1, popcnt(FillinMask(j))
            in = in + 1
            c = popcnt( iand(FillinMask(j), - FillinMask(j)) -1 )
            FillinMask(j) = ibclr(FillinMask(j),c)
            Array(in) = bit*(j-1)+c
          enddo
        enddo
        deallocate(T(parent)%ancestor)
        T(parent)%ancestor => Array
      endif
    enddo
    deallocate(ChildMask )
    deallocate(ParentMask)
    deallocate(FillinMask)

    do i = 1, N
      do j = 1, T(i)%n_ancestor
        T(T(i)%ancestor(j))%n_descendant = T(T(i)%ancestor(j))%n_descendant + 1
      enddo
    enddo
    do i = 1, N
      allocate(T(i)%descendant(T(i)%n_descendant))
    enddo

    allocate(diff(N))
    diff(:) = 0
    do i = 1, N
      do j = 1, T(i)%n_ancestor
        diff(T(i)%ancestor(j)) = diff(T(i)%ancestor(j)) + 1
        T(T(i)%ancestor(j))%descendant(diff(T(i)%ancestor(j))) = i
      enddo
    enddo

    allocate(idxU(0:N))
    in = 0
    idxU(0) = 0
    do i=1,N
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

  subroutine monolis_matrix_copy_with_fillin(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint), pointer :: idxU(:)
    integer(kind=kint), pointer :: itemU(:)
    integer(kind=kint), pointer :: idxL(:)
    integer(kind=kint), pointer :: itemL(:)
    real(kind=kreal), pointer :: AU(:)
    real(kind=kreal), pointer :: AL(:)
    integer(kind=kint) :: N, NPU, NPL
    integer(kind=kint) :: i, j, k, iS, iE, jS, jE
    integer(kind=kint) :: in, jn ,kn, nn
    logical :: is_asym = .false.

    N   = hecT%N

    !value
    allocate(AU(9*NPU))
    AU = 0.0d0

    do k=1,N
      in = idxU(k-1)+1
      AU(9*in-8) = hecT%D(9*k-8)
      AU(9*in-7) = hecT%D(9*k-7)
      AU(9*in-6) = hecT%D(9*k-6)
      AU(9*in-5) = hecT%D(9*k-5)
      AU(9*in-4) = hecT%D(9*k-4)
      AU(9*in-3) = hecT%D(9*k-3)
      AU(9*in-2) = hecT%D(9*k-2)
      AU(9*in-1) = hecT%D(9*k-1)
      AU(9*in  ) = hecT%D(9*k  )
    enddo

    do k=1,N
      iS = idxU(k-1) + 1
      iE = idxU(k)
      jS = hecT%indexU(k-1) + 1
      jE = hecT%indexU(k)
      aa:do j=jS,jE
        jn = hecT%itemU(j)
        do i=iS,iE
          in = itemU(i)
          if(jn == in)then
            AU(9*i-8) = hecT%AU(9*j-8)
            AU(9*i-7) = hecT%AU(9*j-7)
            AU(9*i-6) = hecT%AU(9*j-6)
            AU(9*i-5) = hecT%AU(9*j-5)
            AU(9*i-4) = hecT%AU(9*j-4)
            AU(9*i-3) = hecT%AU(9*j-3)
            AU(9*i-2) = hecT%AU(9*j-2)
            AU(9*i-1) = hecT%AU(9*j-1)
            AU(9*i  ) = hecT%AU(9*j  )
            iS = iS + 1
            cycle aa
          endif
        enddo
      enddo aa
    enddo

    if(is_asym)then
      allocate(AL(9*NPL))
      AL = 0.0d0

      do k=1,N
        iS = idxL(k-1) + 1
        iE = idxL(k)
        jS = hecT%indexL(k-1) + 1
        jE = hecT%indexL(k)
        bb:do j=jS,jE
          jn = hecT%itemL(j)
          do i=iS,iE
            in = itemL(i)
            if(jn == in)then
              AL(9*i-8) = hecT%AL(9*j-8)
              AL(9*i-7) = hecT%AL(9*j-7)
              AL(9*i-6) = hecT%AL(9*j-6)
              AL(9*i-5) = hecT%AL(9*j-5)
              AL(9*i-4) = hecT%AL(9*j-4)
              AL(9*i-3) = hecT%AL(9*j-3)
              AL(9*i-2) = hecT%AL(9*j-2)
              AL(9*i-1) = hecT%AL(9*j-1)
              AL(9*i  ) = hecT%AL(9*j  )
              iS = iS + 1
              cycle bb
            endif
          enddo
        enddo bb
      enddo
    endif
  end subroutine monolis_matrix_copy_with_fillin
end module mod_monolis_matrix_fillin