module mod_monolis_reorder
  use, intrinsic :: iso_c_binding
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_metis
  implicit none

  type monolis_edge_info
    integer(kind=kint) :: N = 0
    integer(kind=kint), pointer :: node(:) => null()
  endtype monolis_edge_info

contains

  subroutine monolis_reorde_matrix_metis(N, NP, indexL, itemL, indexU, itemU, perm, iperm)
    implicit none
    type(monolis_edge_info), allocatable :: edge(:)
    integer(kind=kint), intent(in)  :: N, NP
    integer(kind=kint), intent(in)  :: indexL(0:), indexU(0:)
    integer(kind=kint), intent(in)  :: itemL(:), itemU(:)
    integer(kind=kint), intent(out) :: perm(:), iperm(:)
    integer(kind=kint) :: i, j, k, iS, jS, jE, nu, nl
    integer(kind=kint) :: icel, nedge, ic_type, in, jn, kn, nn, ne
    integer(kind=kint) :: imax, imin
    integer(kind=kint) :: nlocal(20)
    integer(kind=kint), allocatable :: check(:), nozero(:)
    integer(idx_t) :: ierr
    integer(idx_t) :: nvtxs
    integer(idx_t), pointer :: xadj(:)        => null()
    integer(idx_t), pointer :: adjncy(:)      => null()
    integer(idx_t), pointer :: vwgt(:)        => null()
    integer(idx_t), pointer :: options(:)     => null()
    integer(idx_t), pointer :: metis_perm(:)  => null()
    integer(idx_t), pointer :: metis_iperm(:) => null()

    nvtxs = N
    allocate(edge(N))

    do i = 1, N
      nl = indexL(i) - indexL(i-1)
      nu = indexU(i) - indexU(i-1)
      in = nl + nu

      if(0 < in)then
        allocate(nozero(in))
        nozero = 0

        jn = 1
        jS = indexL(i-1) + 1
        jE = indexL(i  )
        do j = jS, jE
          nozero(jn) = itemL(j)
          jn = jn + 1
        enddo
        jS = indexU(i-1) + 1
        jE = indexU(i  )
        do j = jS, jE
          if(N < itemU(j))then
            nozero(jn) = i
          else
            nozero(jn) = itemU(j)
          endif
          jn = jn + 1
        enddo

        call reallocate_array(edge(i)%N, in, edge(i)%node)
        edge(i)%N = in

        do j = 1, in
          edge(i)%node(j) = nozero(j)
        enddo
        deallocate(nozero)
      endif
    enddo

    !buget sort (delete multiple entries)
    do i = 1, N
      in = edge(i)%N
      if(0 < in)then
        imax = maxval(edge(i)%node)
        imin = minval(edge(i)%node)
        allocate(check(imin:imax))
        check(imin:imax) = 0
        do j = 1, in
          check(edge(i)%node(j)) = 1
        enddo
        nn = 0
        do j = imin, imax
          if(check(j) == 1) nn = nn + 1
        enddo
        edge(i)%N = nn
        deallocate(edge(i)%node)
        allocate(edge(i)%node(nn))
        in = 1
        do j = imin, imax
          if(check(j) == 1)then
            edge(i)%node(in) = j
            in = in + 1
          endif
        enddo
        deallocate(check)
      endif
    enddo

    allocate(xadj(N+1))
    xadj(1) = 0
    do i = 1, N
      xadj(i+1) = xadj(i) + edge(i)%N
    enddo

    nedge = xadj(N+1)
    allocate(adjncy(nedge))
    in = 1
    do i = 1, N
      do j = 1, edge(i)%N
        adjncy(in) = edge(i)%node(j) - 1
        in = in + 1
      enddo
    enddo

    allocate(metis_perm(N))
    allocate(metis_iperm(N))

    ierr = monolis_metis_nodeND(nvtxs, xadj, adjncy, vwgt, options, metis_perm, metis_iperm)
    do i = 1, N
      perm (i) = metis_perm(i)  + 1
      iperm(i) = metis_iperm(i) + 1
    enddo

    deallocate(edge)
    deallocate(check)
    do i = 1, N
      deallocate(edge(i)%node)
    enddo
    deallocate(edge)
    deallocate(xadj)
    deallocate(adjncy)
    deallocate(metis_perm)
    deallocate(metis_iperm)
  end subroutine monolis_reorde_matrix_metis

  subroutine reallocate_array(iold, inew, x)
    implicit none
    integer(kind=kint), intent(in) :: iold, inew
    integer(kind=kint), pointer :: x(:), t(:)
    integer(kind=kint) :: i

    if(.not. associated(x))then
      allocate(x(inew))
    else
      t => x
      x => null()
      allocate(x(inew))
      do i=1,iold
        x(i) = t(i)
      enddo
      deallocate(t)
    endif
  end subroutine reallocate_array
end module mod_monolis_reorder