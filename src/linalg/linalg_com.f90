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
      call MPI_allREDUCE(val, temp, N, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
    elseif(tag == tagMax)then
      call MPI_allREDUCE(val, temp, N, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)
    elseif(tag == tagMin)then
      call MPI_allREDUCE(val, temp, N, MPI_DOUBLE_PRECISION, MPI_MIN, comm, ierr)
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
      call MPI_allREDUCE(val, temp, N, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
    elseif(tag == tagMax)then
      call MPI_allREDUCE(val, temp, N, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)
    elseif(tag == tagMin)then
      call MPI_allREDUCE(val, temp, N, MPI_DOUBLE_PRECISION, MPI_MIN, comm, ierr)
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
      call MPI_allREDUCE(val, temp, N, MPI_INTEGER, MPI_SUM, comm, ierr)
    elseif(tag == tagMax)then
      call MPI_allREDUCE(val, temp, N, MPI_INTEGER, MPI_MAX, comm, ierr)
    elseif(tag == tagMin)then
      call MPI_allREDUCE(val, temp, N, MPI_INTEGER, MPI_MIN, comm, ierr)
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
      call MPI_allREDUCE(val, temp, N, MPI_INTEGER, MPI_SUM, comm, ierr)
    elseif(tag == tagMax)then
      call MPI_allREDUCE(val, temp, N, MPI_INTEGER, MPI_MAX, comm, ierr)
    elseif(tag == tagMin)then
      call MPI_allREDUCE(val, temp, N, MPI_INTEGER, MPI_MIN, comm, ierr)
    endif
    val = temp
#endif
  end subroutine monolis_allreduce_I

  subroutine monolis_SendRecv_R &
     &       ( NEIBPETOT, NEIBPE, SendIndex, SendItem, RecvIndex, RecvItem, &
     &         WS, WR, VAL, comm)
    implicit none
    integer(kind=kint) :: NEIBPETOT
    integer(kind=kint) :: istart, inum, k, ierr
    integer(kind=kint) :: neib, comm
    integer(kind=kint), pointer :: NEIBPE   (:)
    integer(kind=kint), pointer :: SendIndex(:)
    integer(kind=kint), pointer :: SendItem (:)
    integer(kind=kint), pointer :: RecvIndex(:)
    integer(kind=kint), pointer :: RecvItem (:)
    integer(kind=kint), allocatable :: sta1(:,:)
    integer(kind=kint), allocatable :: sta2(:,:)
    integer(kind=kint), allocatable :: req1(:)
    integer(kind=kint), allocatable :: req2(:)
    real(kind=kdouble) :: VAL(:)
    real(kind=kdouble) :: WS(:), WR(:)

#ifdef WITHMPI
    allocate(sta1(MPI_STATUS_SIZE, NEIBPETOT))
    allocate(sta2(MPI_STATUS_SIZE, NEIBPETOT))
    allocate(req1(NEIBPETOT))
    allocate(req2(NEIBPETOT))

    do neib = 1, NEIBPETOT
      istart= SendIndex(neib-1)
      inum  = SendIndex(neib  ) - istart
      do k = istart+1, istart+inum
        WS(k) = VAL(SendItem(k))
      enddo
      call MPI_ISEND(WS(istart+1), inum, MPI_DOUBLE_PRECISION, NEIBPE(neib), 0, comm, req1(neib), ierr)
    enddo

    do neib = 1, NEIBPETOT
      istart= RecvIndex(neib-1)
      inum  = RecvIndex(neib  ) - istart
      call MPI_IRECV(WR(istart+1), inum, MPI_DOUBLE_PRECISION, NEIBPE(neib), 0, comm, req2(neib), ierr)
    enddo

    call MPI_WAITALL(NEIBPETOT, req2, sta2, ierr)

    do neib= 1, NEIBPETOT
      istart= RecvIndex(neib-1)
      inum  = RecvIndex(neib  ) - istart
      do k= istart+1, istart+inum
        VAL(RecvItem(k))= WR(k)
      enddo
    enddo

    call MPI_WAITALL(NEIBPETOT, req1, sta1, ierr)

    deallocate(sta1, sta2, req1, req2)
#endif

  end subroutine monolis_SendRecv_R

  subroutine monolis_SendRecv_I &
     &       ( NEIBPETOT, NEIBPE, SendIndex, SendItem, RecvIndex, RecvItem, &
     &         WS, WR, VAL, comm)
    implicit none
    integer(kind=kint) :: NEIBPETOT
    integer(kind=kint) :: istart, inum, k, ierr
    integer(kind=kint) :: neib, comm
    integer(kind=kint) :: VAL(:)
    integer(kind=kint) :: WS(:), WR(:)
    integer(kind=kint), pointer :: NEIBPE   (:)
    integer(kind=kint), pointer :: SendIndex(:)
    integer(kind=kint), pointer :: SendItem (:)
    integer(kind=kint), pointer :: RecvIndex(:)
    integer(kind=kint), pointer :: RecvItem (:)
    integer(kind=kint), allocatable :: sta1(:,:)
    integer(kind=kint), allocatable :: sta2(:,:)
    integer(kind=kint), allocatable :: req1(:)
    integer(kind=kint), allocatable :: req2(:)

#ifdef WITHMPI
    allocate(sta1(MPI_STATUS_SIZE, NEIBPETOT))
    allocate(sta2(MPI_STATUS_SIZE, NEIBPETOT))
    allocate(req1(NEIBPETOT))
    allocate(req2(NEIBPETOT))

    do neib = 1, NEIBPETOT
      istart= SendIndex(neib-1)
      inum  = SendIndex(neib  ) - istart
      do k = istart+1, istart+inum
        WS(k) = VAL(SendItem(k))
      enddo
      call MPI_ISEND(WS(istart+1), inum, MPI_INTEGER, NEIBPE(neib), 0, comm, req1(neib), ierr)
    enddo

    do neib = 1, NEIBPETOT
      istart= RecvIndex(neib-1)
      inum  = RecvIndex(neib  ) - istart
      call MPI_IRECV(WR(istart+1), inum, MPI_INTEGER, NEIBPE(neib), 0, comm, req2(neib), ierr)
    enddo

    call MPI_WAITALL(NEIBPETOT, req2, sta2, ierr)

    do neib= 1, NEIBPETOT
      istart= RecvIndex(neib-1)
      inum  = RecvIndex(neib  ) - istart
      do k= istart+1, istart+inum
        VAL(RecvItem(k))= WR(k)
      enddo
    enddo

    call MPI_WAITALL(NEIBPETOT, req1, sta1, ierr)

    deallocate (sta1, sta2, req1, req2)
#endif
  end subroutine monolis_SendRecv_I

end module mod_monolis_linalg_com