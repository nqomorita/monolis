module mod_monolis_linalg_com
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  implicit none

  integer(kint), parameter :: tagSum = 1
  integer(kint), parameter :: tagMax = 2
  integer(kint), parameter :: tagMin = 3

contains

  subroutine monolis_allreduce_R1(val, tag, comm)
    implicit none
    integer(kint), intent(in) :: tag, comm
    integer(kint) :: N, ierr
    real(kdouble) :: val, temp

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
    integer(kint), intent(in) :: tag, comm
    integer(kint) :: N, ierr
    real(kdouble) :: val(N), temp(N)

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
    integer(kint), intent(in) :: tag, comm
    integer(kint)  :: N, ierr, val, temp

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
    integer(kint), intent(in) :: tag, comm
    integer(kint)  :: N, ierr, val(N), temp(N)

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

  subroutine monolis_SendRecv_pre_R(send_n_neib, send_neib_pe, recv_n_neib, recv_neib_pe, &
    & send_index, send_item, recv_index, recv_item, &
    & ws, wr, val, ndof, comm)
    implicit none
    integer(kint) :: send_n_neib, recv_n_neib
    integer(kint) :: iS, in, j, k, ierr
    integer(kint) :: i, ndof, comm
    integer(kint), pointer :: send_neib_pe(:)
    integer(kint), pointer :: send_index(:)
    integer(kint), pointer :: send_item (:)
    integer(kint), pointer :: recv_neib_pe(:)
    integer(kint), pointer :: recv_index(:)
    integer(kint), pointer :: recv_item (:)
    integer(kint) :: sta1(monolis_status_size, send_n_neib)
    integer(kint) :: sta2(monolis_status_size, recv_n_neib)
    integer(kint) :: req1(send_n_neib)
    integer(kint) :: req2(recv_n_neib)
    real(kdouble) :: val(:), ws(:), wr(:)

#ifdef WITH_MPI
    do i = 1, send_n_neib
      iS = send_index(i-1)
      in = send_index(i  ) - iS
      if(in == 0) cycle
      do j = iS+1, iS+in
        do k = 1, ndof
          ws(ndof*(j-1)+k) = val(ndof*(send_item(j)-1)+k)
        enddo
      enddo
      call MPI_Isend(ws(ndof*iS+1), ndof*in, MPI_DOUBLE_PRECISION, send_neib_pe(i), 0, comm, req1(i), ierr)
    enddo

    do i = 1, recv_n_neib
      iS = recv_index(i-1)
      in = recv_index(i  ) - iS
      if(in == 0) cycle
      call MPI_Irecv(wr(ndof*iS+1), ndof*in, MPI_DOUBLE_PRECISION, recv_neib_pe(i), 0, comm, req2(i), ierr)
    enddo

    call MPI_waitall(recv_n_neib, req2, sta2, ierr)

    do i = 1, recv_n_neib
      iS = recv_index(i-1)
      in = recv_index(i  ) - iS
      do j = iS+1, iS+in
        do k = 1, ndof
          val(ndof*(recv_item(j)-1)+k) = wr(ndof*(j-1)+k)
        enddo
      enddo
    enddo

    call MPI_waitall(send_n_neib, req1, sta1, ierr)
#endif
  end subroutine monolis_SendRecv_pre_R

  subroutine monolis_SendRecv_post_R(n_neib, neib_pe, send_index, send_item, recv_index, recv_item, &
  & ws, wr, val, ndof, comm)
    implicit none
    integer(kint) :: n_neib
    integer(kint) :: iS, in, j, k, ierr
    integer(kint) :: i, ndof, comm
    integer(kint), pointer :: neib_pe(:)
    integer(kint), pointer :: send_index(:)
    integer(kint), pointer :: send_item (:)
    integer(kint), pointer :: recv_index(:)
    integer(kint), pointer :: recv_item (:)
    integer(kint) :: sta1(monolis_status_size, n_neib)
    integer(kint) :: sta2(monolis_status_size, n_neib)
    integer(kint) :: req1(n_neib)
    integer(kint) :: req2(n_neib)
    real(kdouble) :: val(:), ws(:), wr(:)

#ifdef WITH_MPI

#endif
  end subroutine monolis_SendRecv_post_R

  subroutine monolis_SendRecv_pre_I(send_n_neib, send_neib_pe, recv_n_neib, recv_neib_pe, &
    & send_index, send_item, recv_index, recv_item, &
    & ws, wr, val, ndof, comm)
    implicit none
    integer(kint) :: send_n_neib, recv_n_neib
    integer(kint) :: iS, in, j, k, ierr
    integer(kint) :: i, ndof, comm
    integer(kint), pointer :: send_neib_pe(:)
    integer(kint), pointer :: send_index(:)
    integer(kint), pointer :: send_item (:)
    integer(kint), pointer :: recv_neib_pe(:)
    integer(kint), pointer :: recv_index(:)
    integer(kint), pointer :: recv_item (:)
    integer(kint) :: sta1(monolis_status_size, send_n_neib)
    integer(kint) :: sta2(monolis_status_size, recv_n_neib)
    integer(kint) :: req1(send_n_neib)
    integer(kint) :: req2(recv_n_neib)
    integer(kint) :: val(:), ws(:), wr(:)

#ifdef WITH_MPI
    do i = 1, send_n_neib
      iS = send_index(i-1)
      in = send_index(i  ) - iS
      if(in == 0) cycle
      do j = iS+1, iS+in
        do k = 1, ndof
          ws(ndof*(j-1)+k) = val(ndof*(send_item(j)-1)+k)
        enddo
      enddo
      call MPI_Isend(ws(ndof*iS+1), ndof*in, MPI_DOUBLE_PRECISION, send_neib_pe(i), 0, comm, req1(i), ierr)
    enddo

    do i = 1, recv_n_neib
      iS = recv_index(i-1)
      in = recv_index(i  ) - iS
      if(in == 0) cycle
      call MPI_Irecv(wr(ndof*iS+1), ndof*in, MPI_DOUBLE_PRECISION, recv_neib_pe(i), 0, comm, req2(i), ierr)
    enddo

    call MPI_waitall(recv_n_neib, req2, sta2, ierr)

    do i = 1, recv_n_neib
      iS = recv_index(i-1)
      in = recv_index(i  ) - iS
      do j = iS+1, iS+in
        do k = 1, ndof
          val(ndof*(recv_item(j)-1)+k) = wr(ndof*(j-1)+k)
        enddo
      enddo
    enddo

    call MPI_waitall(send_n_neib, req1, sta1, ierr)
#endif
  end subroutine monolis_SendRecv_pre_I

  subroutine monolis_allgather_I1(sval, rbuf, comm)
    implicit none
    integer(kint) :: sval    !send buffer
    integer(kint) :: rbuf(:) !receive buffer
    integer(kint) :: comm
#ifdef WITH_MPI
    integer(kint) :: ierr
    call MPI_allgather(sval, 1, MPI_INTEGER, rbuf, 1, MPI_INTEGER, comm, ierr )
#else
    rbuf(1) = sval
#endif
  end subroutine monolis_allgather_I1

  subroutine monolis_gatherv_R(sbuf, sc, &
      &     rbuf, rcs, disp, root, comm)
    implicit none
    integer(kint) :: sc      !send counts
    integer(kint) :: rcs(:)  !receive counts
    integer(kint) :: disp(:) !displacement
    integer(kint) :: root
    integer(kint) :: comm
    real(kdouble) :: sbuf(:) !send buffer
    real(kdouble) :: rbuf(:) !receive buffer
#ifdef WITH_MPI
    integer(kint) :: ierr
    call MPI_gatherv( sbuf, sc, MPI_REAL8, &
      &     rbuf, rcs, disp, MPI_REAL8, root, comm, ierr )
#else
    rbuf(1:sc) = sbuf(1:sc)
#endif
  end subroutine monolis_gatherv_R

  subroutine monolis_scatterv_R(sbuf, scs, disp, &
      &     rbuf, rc, root, comm)
    implicit none
    integer(kint) :: scs(:)  !send counts
    integer(kint) :: disp(:) !displacement
    integer(kint) :: rc      !receive counts
    integer(kint) :: root
    integer(kint) :: comm
    real(kdouble) :: rbuf(:) !receive buffer
    real(kdouble) :: sbuf(:) !send buffer
#ifdef WITH_MPI
    integer(kint) :: ierr
    call MPI_scatterv( sbuf, scs, disp, MPI_REAL8, &
      &     rbuf, rc, MPI_REAL8, root, comm, ierr )
#else
    rbuf(1:rc) = sbuf(1:rc)
#endif
  end subroutine monolis_scatterv_R

end module mod_monolis_linalg_com