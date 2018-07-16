module mod_monolis_linalg_com
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  implicit none

  integer(kind=kint), parameter :: tagSum = 1
  integer(kind=kint), parameter :: tagMax = 2
  integer(kind=kint), parameter :: tagMin = 3

contains

  subroutine monolis_allreduce_R1(val, tag, comm)
    implicit none
    integer(kind=kint), intent(in) :: tag, comm
    integer(kind=kint) :: N, ierr
    real(kind=kdouble) :: val, temp


#ifdef WITH_MPI
    temp = 0.0d0
    N = 1
    if(tag == tagSum)then
      call MPI_allreduce(val, temp, N, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
    elseif(tag == tagMax)then
      call MPI_allreduce(val, temp, N, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)
    elseif(tag == tagMin)then
      call MPI_allreduce(val, temp, N, MPI_DOUBLE_PRECISION, MPI_MIN, comm, ierr)
    endif
    val = temp
#endif
  end subroutine monolis_allreduce_R1

  subroutine monolis_allreduce_R(N, val, tag, comm)
    implicit none
    integer(kind=kint), intent(in) :: tag, comm
    integer(kind=kint) :: N, ierr
    real(kind=kdouble) :: val(N), temp(N)

#ifdef WITH_MPI
    temp = 0.0d0
    if(tag == tagSum)then
      call MPI_allreduce(val, temp, N, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
    elseif(tag == tagMax)then
      call MPI_allreduce(val, temp, N, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)
    elseif(tag == tagMin)then
      call MPI_allreduce(val, temp, N, MPI_DOUBLE_PRECISION, MPI_MIN, comm, ierr)
    endif
    val = temp
#endif
  end subroutine monolis_allreduce_R

  subroutine monolis_allreduce_I1(val, tag, comm)
    implicit none
    integer(kind=kint), intent(in) :: tag, comm
    integer(kind=kint)  :: N, ierr, val, temp

#ifdef WITH_MPI
    temp = 0
    N = 1
    if(tag == tagSum)then
      call MPI_allreduce(val, temp, N, MPI_INTEGER, MPI_SUM, comm, ierr)
    elseif(tag == tagMax)then
      call MPI_allreduce(val, temp, N, MPI_INTEGER, MPI_MAX, comm, ierr)
    elseif(tag == tagMin)then
      call MPI_allreduce(val, temp, N, MPI_INTEGER, MPI_MIN, comm, ierr)
    endif
    val = temp
#endif
  end subroutine monolis_allreduce_I1

  subroutine monolis_allreduce_I(N, val, tag, comm)
    implicit none
    integer(kind=kint), intent(in) :: tag, comm
    integer(kind=kint)  :: N, ierr, val(N), temp(N)

#ifdef WITH_MPI
    temp = 0
    if(tag == tagSum)then
      call MPI_allreduce(val, temp, N, MPI_INTEGER, MPI_SUM, comm, ierr)
    elseif(tag == tagMax)then
      call MPI_allreduce(val, temp, N, MPI_INTEGER, MPI_MAX, comm, ierr)
    elseif(tag == tagMin)then
      call MPI_allreduce(val, temp, N, MPI_INTEGER, MPI_MIN, comm, ierr)
    endif
    val = temp
#endif
  end subroutine monolis_allreduce_I

  subroutine monolis_SendRecv_R(n_neib, neib_pe, send_index, send_item, recv_index, recv_item, &
  & ws, wr, val, ndof, comm)
    implicit none
    integer(kind=kint) :: n_neib
    integer(kind=kint) :: iS, in, j, k, ierr
    integer(kind=kint) :: i, ndof, comm
    integer(kind=kint), pointer :: neib_pe(:)
    integer(kind=kint), pointer :: send_index(:)
    integer(kind=kint), pointer :: send_item (:)
    integer(kind=kint), pointer :: recv_index(:)
    integer(kind=kint), pointer :: recv_item (:)
    integer(kind=kint) :: sta1(monolis_status_size, n_neib)
    integer(kind=kint) :: sta2(monolis_status_size, n_neib)
    integer(kind=kint) :: req1(n_neib)
    integer(kind=kint) :: req2(n_neib)
    real(kind=kdouble) :: val(:), ws(:), wr(:)

#ifdef WITH_MPI
    do i = 1, n_neib
      iS = send_index(i-1)
      in = send_index(i  ) - iS
      if(in == 0) cycle
      do j = iS+1, iS+in
        do k = 1, ndof
          ws(ndof*(j-1)+k) = val(ndof*(send_item(j)-1)+k)
        enddo
      enddo
      call MPI_Isend(ws(ndof*iS+1), ndof*in, MPI_DOUBLE_PRECISION, neib_pe(i), 0, comm, req1(i), ierr)
    enddo

    do i = 1, n_neib
      iS = recv_index(i-1)
      in = recv_index(i  ) - iS
      if(in == 0) cycle
      call MPI_Irecv(wr(ndof*iS+1), ndof*in, MPI_DOUBLE_PRECISION, neib_pe(i), 0, comm, req2(i), ierr)
    enddo

    call MPI_waitall(n_neib, req2, sta2, ierr)

    do i = 1, n_neib
      iS = recv_index(i-1)
      in = recv_index(i  ) - iS
      do j = iS+1, iS+in
        do k = 1, ndof
          val(ndof*(recv_item(j)-1)+k) = wr(ndof*(j-1)+k)
        enddo
      enddo
    enddo

    call MPI_waitall(n_neib, req1, sta1, ierr)
#endif
  end subroutine monolis_SendRecv_R

  subroutine monolis_SendRecv_I(n_neib, neib_pe, send_index, send_item, recv_index, recv_item, &
  & ws, wr, val, ndof, comm)
    implicit none
    integer(kind=kint) :: n_neib
    integer(kind=kint) :: iS, in, j, k, ierr
    integer(kind=kint) :: i, ndof, comm
    integer(kind=kint) :: val(:), ws(:), wr(:)
    integer(kind=kint), pointer :: neib_pe(:)
    integer(kind=kint), pointer :: send_index(:)
    integer(kind=kint), pointer :: send_item (:)
    integer(kind=kint), pointer :: recv_index(:)
    integer(kind=kint), pointer :: recv_item (:)
    integer(kind=kint) :: sta1(monolis_status_size, n_neib)
    integer(kind=kint) :: sta2(monolis_status_size, n_neib)
    integer(kind=kint) :: req1(n_neib)
    integer(kind=kint) :: req2(n_neib)

#ifdef WITH_MPI
    do i = 1, n_neib
      iS = send_index(i-1)
      in = send_index(i  ) - iS
      if(in == 0) cycle
      do j = iS+1, iS+in
        do k = 1, ndof
          ws(ndof*(j-1)+k) = val(ndof*(send_item(j)-1)+k)
        enddo
      enddo
      call MPI_Isend(ws(ndof*iS+1), ndof*in, MPI_INTEGER, neib_pe(i), 0, comm, req1(i), ierr)
    enddo

    do i = 1, n_neib
      iS = recv_index(i-1)
      in = recv_index(i  ) - iS
      if(in == 0) cycle
      call MPI_Irecv(wr(ndof*iS+1), ndof*in, MPI_INTEGER, neib_pe(i), 0, comm, req2(i), ierr)
    enddo

    call MPI_waitall(n_neib, req2, sta2, ierr)

    do i = 1, n_neib
      iS = recv_index(i-1)
      in = recv_index(i  ) - iS
      do j = iS+1, iS+in
        do k = 1, ndof
          val(ndof*(recv_item(j)-1)+k) = wr(ndof*(j-1)+k)
        enddo
      enddo
    enddo

    call MPI_waitall(n_neib, req1, sta1, ierr)
#endif
  end subroutine monolis_SendRecv_I

end module mod_monolis_linalg_com