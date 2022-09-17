module mod_monolis_util_com
  use mod_monolis_util
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_stdlib
  implicit none

contains

  subroutine monolis_barrier(monolis)
    implicit none
    type(monolis_structure) :: monolis

#ifdef WITH_MPI
    call monolis_barrier_(monolis%COM%comm)
#endif
  end subroutine monolis_barrier

  subroutine monolis_com_set_communicator(monolis, comm)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: comm
    monolis%COM%comm = comm
  end subroutine monolis_com_set_communicator

  subroutine monolis_com_set_myrank(monolis, myrank)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: myrank
    monolis%COM%myrank = myrank
  end subroutine monolis_com_set_myrank

  subroutine monolis_com_set_commsize(monolis, commsize)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: commsize
    monolis%COM%commsize = commsize
  end subroutine monolis_com_set_commsize

  subroutine monolis_com_get_comm_table(monolis, N, NP, nid)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: N, NP, M, myrank, commsize, ierr
    integer(kint) :: i, id, idx, j, jS, jE
    integer(kint) :: nid(:)
    integer(kint), allocatable :: counts(:), outer_id_local(:), outer_id_all(:)
    integer(kint), allocatable :: displs(:), internal_id(:), is_neib(:)

    !if(monolis%MAT%NP /= size(nid))then
    !  write(*,*)"*** ERROR: The number of matrix DoF and the number of global column id are different."
    !  stop
    !endif

    write(*,*)size(nid)
    write(*,*)nid

    myrank = monolis_global_myrank()
    commsize = monolis_global_commsize()

    M = NP - N

    !> 外点を送信
    allocate(counts(commsize), source = 0)

    call mpi_allgather(M, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, monolis%COM%comm, ierr)

    allocate(outer_id_local(M), source = 0)
    allocate(displs(commsize+1), source = 0)

    do i = N+1, NP
      outer_id_local(i-N) = nid(i)
    enddo

    do i = 1, commsize
      displs(i+1) = displs(i) + counts(i)
    enddo

    write(100+monolis_global_myrank(),*)"counts"
    write(100+monolis_global_myrank(),*)counts
    write(100+monolis_global_myrank(),*)"displs"
    write(100+monolis_global_myrank(),*)displs

    allocate(outer_id_all(displs(commsize+1)), source = 0)

    call mpi_allgatherv(outer_id_local, M, MPI_INTEGER, &
      outer_id_all, counts, displs, MPI_INTEGER, monolis%COM%comm, ierr)

    write(100+monolis_global_myrank(),*)"outer_id_all"
    write(100+monolis_global_myrank(),*)outer_id_all

    !> 隣接領域チェック
    allocate(internal_id(N), source = 0)
    do i = 1, N
      internal_id(i) = nid(i)
    enddo
    call monolis_qsort_int(internal_id, 1, N)

    allocate(is_neib(commsize), source = 0)

    aa:do i = 1, commsize
      if(i == myrank+1) cycle
      jS = displs(i) + 1
      jE = displs(i+1)
      do j = jS, jE
        id = outer_id_all(j)
        call monolis_bsearch_int(internal_id, 1, N, id, idx)
        is_neib(i) = 1
        if(idx /= -1) cycle aa
      enddo
    enddo aa

    write(100+monolis_global_myrank(),*)"neib"
    write(100+monolis_global_myrank(),*)is_neib

    !> slave 側の情報を master 側に送信（recv を作成）

    !> master 側の情報を整理（send を作成）

  end subroutine monolis_com_get_comm_table

end module mod_monolis_util_com
