module mod_monolis_restruct
  use mod_monolis_prm
  use mod_monolis_mat

  implicit none
  private
  public :: monolis_restruct_comm
  public :: monolis_restruct_matrix

contains

  subroutine monolis_restruct_comm(monoCOM, monoCOM_reorder, perm)
    implicit none
    type(monolis_com) :: monoCOM
    type(monolis_com) :: monoCOM_reorder
    integer(kind=kint) :: perm(:)
    integer(kind=kint) :: i, in, N

!    monoCOM_reorder%myrank   = monoCOM%myrank
!    monoCOM_reorder%comm     = monoCOM%comm
!    monoCOM_reorder%commsize = monoCOM%commsize
!    monoCOM_reorder%send_n_neib = monoCOM%send_n_neib
!    monoCOM_reorder%recv_n_neib = monoCOM%recv_n_neib
!
!    if(monoCOM%send_n_neib /= 0)then
!      monoCOM_reorder%recv_neib_pe => monoCOM%recv_neib_pe
!      monoCOM_reorder%recv_index   => monoCOM%recv_index
!      monoCOM_reorder%recv_item    => monoCOM%recv_item
!
!      monoCOM_reorder%send_neib_pe => monoCOM%send_neib_pe
!      monoCOM_reorder%send_index   => monoCOM%send_index
!
!      N = monoCOM%send_index(monoCOM%send_n_neib)
!      allocate(monoCOM_reorder%send_item(N))
!      do i = 1, N
!        in = perm(monoCOM%send_item(i))
!        monoCOM_reorder%send_item(i) = in
!      enddo
!    endif
  end subroutine monolis_restruct_comm

  subroutine monolis_restruct_matrix(monoMAT, monoMAT_reorder, perm, iperm)
    implicit none
    type(monolis_mat) :: monoMAT
    type(monolis_mat) :: monoMAT_reorder
    integer(kind=kint) :: perm(:), iperm(:)
    integer(kind=kint) :: N, NP, NZ, NDOF, NDOF2

!    N = monoMAT%N
!    NP = monoMAT%NP
!    NZ = monoMAT%index(NP)
!    NDOF = monoMAT%NDOF
!    NDOF2 = NDOF*NDOF
!
!    monoMAT_reorder%N = N
!    monoMAT_reorder%NP = NP
!    monoMAT_reorder%NZ = NZ
!    monoMAT_reorder%NDOF = NDOF
!    allocate(monoMAT_reorder%index(0:NP))
!    allocate(monoMAT_reorder%item(NZ))
!    call monolis_restruct_matrix_profile(NP, perm, iperm, &
!       & monoMAT%index, monoMAT%item, monoMAT_reorder%index, monoMAT_reorder%item)
!
!    allocate(monoMAT_reorder%A(NDOF2*NZ))
!    call monolis_restruct_matrix_values(NP, NDOF, perm, iperm, &
!       & monoMAT%index, monoMAT%item, monoMAT%A, &
!       & monoMAT_reorder%index, monoMAT_reorder%item, monoMAT_reorder%A)
!
!    allocate(monoMAT_reorder%X(NDOF*NP))
!    allocate(monoMAT_reorder%B(NDOF*NP))
  end subroutine monolis_restruct_matrix

  subroutine monolis_restruct_matrix_profile(N, perm, iperm, &
    & index, item, indexp, itemp)
    implicit none
    integer(kind=kint) :: N
    integer(kind=kint) :: perm(:), iperm(:)
    integer(kind=kint) :: index(0:), item(:)
    integer(kind=kint) :: indexp(0:), itemp(:)
    integer(kind=kint) :: cnt, i, in, j, jo, jn

!    cnt = 0
!    indexp(0) = 0
!    do i = 1, N
!      in = perm(i)
!      do j = index(in-1)+1, index(in)
!        jo = item(j)
!        jn = iperm(jo)
!        cnt = cnt + 1
!        itemp(cnt) = jn
!      enddo
!      indexp(i) = cnt
!      call sort_int_array(itemp, indexp(i-1)+1, indexp(i))
!    enddo
  end subroutine monolis_restruct_matrix_profile

  subroutine monolis_restruct_matrix_values(N, NDOF, perm, iperm, index, item, A, &
      & indexp, itemp, Ap)
    implicit none
    integer(kind=kint) :: N, NDOF
    integer(kind=kint) :: perm(:), iperm(:)
    integer(kind=kint) :: index(0:), item(:)
    real(kind=kdouble) :: A(:)
    integer(kind=kint) :: indexp(0:), itemp(:)
    real(kind=kdouble) :: Ap(:)
    integer(kind=kint) :: NDOF2, in, i
    integer(kind=kint) :: jSn, jEn
    integer(kind=kint) :: jo, ko, kn, jn, lo, ln, l

!    NDOF2 = NDOF*NDOF
!    do i = 1, N
!      in = iperm(i)
!      jSn = indexp(in-1)+1
!      jEn = indexp(in)
!      do jo = index(i-1)+1, index(i)
!        ko = item(jo)
!        kn = iperm(ko)
!        call bsearch_int_array(itemp, jSn, jEn, kn, jn)
!        lo = (jo-1)*NDOF2
!        ln = (jn-1)*NDOF2
!        do l = 1, NDOF2
!          Ap(ln + l) = A(lo + l)
!        enddo
!      enddo
!    enddo
  end subroutine monolis_restruct_matrix_values
end module mod_monolis_restruct