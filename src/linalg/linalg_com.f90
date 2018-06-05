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

#ifdef WITHMPI
    temp = 0.0d0
    N    = 1
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

#ifdef WITHMPI
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

#ifdef WITHMPI
    temp = 0
    N    = 1
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

#ifdef WITHMPI
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
  & ws, wr, val, comm)
    implicit none
    integer(kind=kint) :: n_neib
    integer(kind=kint) :: istart, inum, k, ierr
    integer(kind=kint) :: neib, comm
    integer(kind=kint), pointer :: neib_pe(:)
    integer(kind=kint), pointer :: send_index(:)
    integer(kind=kint), pointer :: send_item (:)
    integer(kind=kint), pointer :: recv_index(:)
    integer(kind=kint), pointer :: recv_item (:)
    integer(kind=kint) :: sta1(MPI_STATUS_SIZE, n_neib)
    integer(kind=kint) :: sta2(MPI_STATUS_SIZE, n_neib)
    integer(kind=kint) :: req1(n_neib)
    integer(kind=kint) :: req2(n_neib)
    real(kind=kdouble) :: val(:), ws(:), wr(:)

#ifdef WITHMPI
    do neib = 1, n_neib
      istart= send_index(neib-1)
      inum  = send_index(neib  ) - istart
      do k = istart+1, istart+inum
        ws(k) = val(send_item(k))
      enddo
      call MPI_Isend(ws(istart+1), inum, MPI_DOUBLE_PRECISION, neib_pe(neib), 0, comm, req1(neib), ierr)
    enddo

    do neib = 1, n_neib
      istart= recv_index(neib-1)
      inum  = recv_index(neib  ) - istart
      call MPI_Irecv(wr(istart+1), inum, MPI_DOUBLE_PRECISION, neib_pe(neib), 0, comm, req2(neib), ierr)
    enddo

    call MPI_WAITALL(n_neib, req2, sta2, ierr)

    do neib= 1, n_neib
      istart= recv_index(neib-1)
      inum  = recv_index(neib  ) - istart
      do k= istart+1, istart+inum
        val(recv_item(k))= wr(k)
      enddo
    enddo

    call MPI_waitall(n_neib, req1, sta1, ierr)
#endif
  end subroutine monolis_SendRecv_R

  subroutine monolis_SendRecv_I(n_neib, neib_pe, send_index, send_item, recv_index, recv_item, &
  & ws, wr, val, comm)
    implicit none
    integer(kind=kint) :: n_neib
    integer(kind=kint) :: istart, inum, k, ierr
    integer(kind=kint) :: neib, comm
    integer(kind=kint) :: val(:), ws(:), wr(:)
    integer(kind=kint), pointer :: neib_pe(:)
    integer(kind=kint), pointer :: send_index(:)
    integer(kind=kint), pointer :: send_item (:)
    integer(kind=kint), pointer :: recv_index(:)
    integer(kind=kint), pointer :: recv_item (:)
    integer(kind=kint) :: sta1(MPI_STATUS_SIZE, n_neib)
    integer(kind=kint) :: sta2(MPI_STATUS_SIZE, n_neib)
    integer(kind=kint) :: req1(n_neib)
    integer(kind=kint) :: req2(n_neib)

#ifdef WITHMPI
    do neib = 1, n_neib
      istart= send_index(neib-1)
      inum  = send_index(neib  ) - istart
      do k = istart+1, istart+inum
        ws(k) = val(send_item(k))
      enddo
      call MPI_Isend(ws(istart+1), inum, MPI_INTEGER, neib_pe(neib), 0, comm, req1(neib), ierr)
    enddo

    do neib = 1, n_neib
      istart= recv_index(neib-1)
      inum  = recv_index(neib  ) - istart
      call MPI_Irecv(wr(istart+1), inum, MPI_INTEGER, neib_pe(neib), 0, comm, req2(neib), ierr)
    enddo

    call MPI_WAITALL(n_neib, req2, sta2, ierr)

    do neib= 1, n_neib
      istart= recv_index(neib-1)
      inum  = recv_index(neib  ) - istart
      do k= istart+1, istart+inum
        val(recv_item(k))= wr(k)
      enddo
    enddo

    call MPI_waitall(n_neib, req1, sta1, ierr)
#endif
  end subroutine monolis_SendRecv_I

end module mod_monolis_linalg_com