module mod_monolis_util_com
  use mod_monolis_util
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_linalg_com
  use mod_monolis_stdlib
  implicit none

  private

  public :: monolis_barrier
  public :: monolis_com_set_communicator
  public :: monolis_com_set_myrank
  public :: monolis_com_set_commsize
  public :: monolis_com_get_comm_table

  type monolis_node_list
    integer(kint) :: nnode = 0
    integer(kint) :: domid = -1
    !integer(kint), allocatable :: domid(:)
    integer(kint), allocatable :: local_nid(:)
  end type monolis_node_list

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
    type(monolis_node_list), allocatable :: send_list(:)
    type(monolis_node_list), allocatable :: recv_list(:)
    integer(kint), intent(in) :: N, NP
    integer(kint), intent(in) :: nid(:)
    integer(kint) :: M, n_outer, myrank, commsize, ierr
    integer(kint) :: i, in, id, idx, j, jS, jE, recv_rank, ns, nr
    integer(kint) :: n_neib_recv, n_neib_send, local_id, global_id, n_data
    integer(kint), allocatable :: counts(:), outer_node_id_local(:), local_nid(:)
    integer(kint), allocatable :: outer_node_id_all(:), outer_dom_id_all(:), temp(:)
    integer(kint), allocatable :: displs(:), internal_node_id(:), is_neib(:), neib_id(:)
    integer(kint), allocatable :: send_n_list(:)
    integer(kint), allocatable :: ws(:), wr(:)
    integer(kint), allocatable :: sta1(:,:)
    integer(kint), allocatable :: sta2(:,:)
    integer(kint), allocatable :: req1(:)
    integer(kint), allocatable :: req2(:)

    !if(monolis%MAT%NP /= size(nid))then
    !  write(*,*)"*** ERROR: The number of matrix DoF and the number of global column id are different."
    !  stop
    !endif

    write(100+monolis_global_myrank(),*)size(nid)
    write(100+monolis_global_myrank(),*)nid

    myrank = monolis_global_myrank()
    commsize = monolis_global_commsize()

    M = NP - N

    !> 外点を全体で共有
    allocate(counts(commsize), source = 0)

    call mpi_allgather(M, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, monolis%COM%comm, ierr)

    allocate(outer_node_id_local(M), source = 0)
    allocate(displs(commsize + 1), source = 0)

    do i = N + 1, NP
      outer_node_id_local(i-N) = nid(i)
    enddo

    do i = 1, commsize
      displs(i + 1) = displs(i) + counts(i)
    enddo

    write(100+monolis_global_myrank(),*)"counts"
    write(100+monolis_global_myrank(),*)counts
    write(100+monolis_global_myrank(),*)"displs"
    write(100+monolis_global_myrank(),*)displs

    n_outer = displs(commsize + 1)
    allocate(outer_node_id_all(n_outer), source = 0)

    call mpi_allgatherv(outer_node_id_local, M, MPI_INTEGER, &
      outer_node_id_all, counts, displs, MPI_INTEGER, monolis%COM%comm, ierr)

    write(100+monolis_global_myrank(),*)"outer_node_id_all"
    write(100+monolis_global_myrank(),*)outer_node_id_all

    !> 外点が属する領域番号を取得
    allocate(internal_node_id(N), source = 0)
    do i = 1, N
      internal_node_id(i) = nid(i)
    enddo
    call monolis_qsort_int(internal_node_id, 1, N)

    allocate(outer_dom_id_all(n_outer), source = 0)
    outer_dom_id_all(:) = commsize + 1

    aa:do i = 1, commsize
      !> 自領域であればスキップ
      if(i == myrank + 1) cycle
      !> 他領域の外点と自領域の内点が重複するか判定
      jS = displs(i) + 1
      jE = displs(i + 1)
      do j = jS, jE
        id = outer_node_id_all(j)
        call monolis_bsearch_int(internal_node_id, 1, N, id, idx)
        if(idx /= -1)then
          outer_dom_id_all(j) = myrank
        endif
      enddo
    enddo aa

    !> reduce に修正して効率化可能
    call monolis_allreduce_I(n_outer, outer_dom_id_all, monolis_min, monolis%COM%comm)

    write(100+monolis_global_myrank(),*)"outer_dom_id_all"
    write(100+monolis_global_myrank(),*)outer_dom_id_all

    !> 隣接領域の取得
    allocate(is_neib(commsize), source = 0)

    do i = 1, n_outer
      j = outer_dom_id_all(i)
      is_neib(j + 1) = 1
    enddo
    is_neib(myrank + 1) = 0

    n_neib_recv = 0
    do i = 1, commsize
      if(is_neib(i) == 1) n_neib_recv = n_neib_recv + 1
    enddo

    allocate(neib_id(n_neib_recv), source = 0)

    j = 0
    do i = 1, commsize
      if(is_neib(i) == 1)then
        j = j + 1
        neib_id(j) = i - 1
      endif
    enddo

    write(100+monolis_global_myrank(),*)"n_neib_recv"
    write(100+monolis_global_myrank(),*)n_neib_recv
    write(100+monolis_global_myrank(),*)"neib_id"
    write(100+monolis_global_myrank(),*)neib_id

    allocate(recv_list(n_neib_recv))

    !> recv の作成
    allocate(local_nid(NP), source = 0)
    allocate(temp(NP), source = 0)

    temp(:) = nid(:)
    do i = 1, NP
      local_nid(i) = i
    enddo

    call monolis_qsort_int_with_perm(temp, 1, NP, local_nid)

    do i = 1, n_neib_recv
      recv_rank = neib_id(i)
      in = myrank + 1
      jS = displs(in) + 1
      jE = displs(in + 1)

      n_data = 0
      do j = jS, jE
        id = outer_dom_id_all(j)
        if(recv_rank == id)then
          n_data = n_data + 1
        endif
      enddo

      recv_list(i)%nnode = n_data
      recv_list(i)%domid = recv_rank
      allocate(recv_list(i)%local_nid(n_data), source = 0)

      n_data = 0
      do j = jS, jE
        id = outer_dom_id_all(j)
        if(recv_rank == id)then
          n_data = n_data + 1
          global_id = outer_node_id_all(j)
          call monolis_bsearch_int(temp, 1, NP, global_id, idx)
          recv_list(i)%local_nid(n_data) = local_nid(idx)
        endif
      enddo

      write(100+monolis_global_myrank(),*)"recv_list(i)%nnode"
      write(100+monolis_global_myrank(),*)recv_list(i)%nnode
      write(100+monolis_global_myrank(),*)"recv_list(i)%domid"
      write(100+monolis_global_myrank(),*)recv_list(i)%domid
      write(100+monolis_global_myrank(),*)"recv_list(i)%local_nid"
      write(100+monolis_global_myrank(),*)recv_list(i)%local_nid
    enddo

    !> send の作成
    !> slave から master に個数を送信
    allocate(send_n_list(commsize), source = 0)

    do i = 1, n_neib_recv
      id = recv_list(i)%domid
      in = recv_list(i)%nnode
      send_n_list(id + 1) = in
    enddo

    call mpi_alltoall(send_n_list, 1, MPI_INTEGER, &
      send_n_list, 1, MPI_INTEGER, monolis%COM%comm, ierr)

    write(100+monolis_global_myrank(),*)"send_n_list"
    write(100+monolis_global_myrank(),*)send_n_list

    !> send 個数の確保
    n_neib_send = 0
    do i = 1, commsize
      if(send_n_list(i) > 0) n_neib_send = n_neib_send + 1
    enddo

    write(100+monolis_global_myrank(),*)"n_neib_send"
    write(100+monolis_global_myrank(),*)n_neib_send

    allocate(send_list(n_neib_send))

    n_neib_send = 0
    do i = 1, commsize
      if(send_n_list(i) > 0)then
        n_neib_send = n_neib_send + 1
        send_list(n_neib_send)%domid = i - 1
        n_data = send_n_list(i)
        send_list(n_neib_send)%nnode = n_data
        allocate(send_list(n_neib_send)%local_nid(n_data), source = 0)
        !write(100+monolis_global_myrank(),*)"send_list(i)%nnode"
        !write(100+monolis_global_myrank(),*)send_list(n_neib_send)%nnode
        !write(100+monolis_global_myrank(),*)"send_list(i)%domid"
        !write(100+monolis_global_myrank(),*)send_list(n_neib_send)%domid
      endif
    enddo

    !> monolis com の構築
    !> recv
    monolis%COM%recv_n_neib = n_neib_recv
    allocate(monolis%COM%recv_neib_pe(n_neib_recv), source = 0)
    do i = 1, n_neib_recv
      monolis%COM%recv_neib_pe(i) = recv_list(i)%domid
    enddo
    allocate(monolis%COM%recv_index(0:n_neib_recv), source = 0)
    do i = 1, n_neib_recv
      monolis%COM%recv_index(i) = monolis%COM%recv_index(i-1) + recv_list(i)%nnode
    enddo
    in = monolis%COM%recv_index(n_neib_recv)
    allocate(monolis%COM%recv_item(in), source = 0)
    in = 0
    do i = 1, n_neib_recv
      jE = recv_list(i)%nnode
      do j = 1, jE
        in = in + 1
        idx = recv_list(i)%local_nid(j)
        monolis%COM%recv_item(in) = idx
      enddo
    enddo

    !> send
    monolis%COM%send_n_neib = n_neib_send
    allocate(monolis%COM%send_neib_pe(n_neib_send), source = 0)
    do i = 1, n_neib_send
      monolis%COM%send_neib_pe(i) = send_list(i)%domid
    enddo
    allocate(monolis%COM%send_index(0:n_neib_send), source = 0)
    do i = 1, n_neib_send
      monolis%COM%send_index(i) = monolis%COM%send_index(i-1) + send_list(i)%nnode
    enddo
    in = monolis%COM%send_index(n_neib_send)
    allocate(monolis%COM%send_item(in), source = 0)

    !> slave から master に global_nid を送信
    !allocate(ws(ns))
    !allocate(wr(nr))
    allocate(sta1(monolis_status_size,monolis%COM%recv_n_neib))
    allocate(sta2(monolis_status_size,monolis%COM%send_n_neib))
    allocate(req1(monolis%COM%recv_n_neib))
    allocate(req2(monolis%COM%send_n_neib))

    in = monolis%COM%recv_index(n_neib_recv)
    allocate(ws(in), source = 0)

    do i = 1, monolis%COM%recv_n_neib
      id = recv_list(i)%domid
      in = recv_list(i)%nnode
      jS = monolis%COM%recv_index(i-1) + 1
      jE = monolis%COM%recv_index(i)
      do j = jS, jE
        idx = monolis%COM%recv_item(j)
        ws(j) = nid(idx)
      enddo
      call MPI_Isend(ws(jS:jE), in, MPI_INTEGER, id, 0, monolis%COM%comm, req1(i), ierr)
    enddo

    in = monolis%COM%send_index(n_neib_send)
    allocate(wr(in), source = 0)
    do i = 1, monolis%COM%send_n_neib
      id = send_list(i)%domid
      in = send_list(i)%nnode
      jS = monolis%COM%send_index(i-1) + 1
      jE = monolis%COM%send_index(i)
      call MPI_Irecv(wr(jS:jE), in, MPI_INTEGER, id, 0, monolis%COM%comm, req2(i), ierr)
    enddo

    call MPI_waitall(monolis%COM%recv_n_neib, req2, sta2, ierr)

    write(100+monolis_global_myrank(),*)"wr"
    write(100+monolis_global_myrank(),*)wr

    call MPI_waitall(monolis%COM%send_n_neib, req1, sta1, ierr)

    !> local_nid に変換
    in = monolis%COM%send_index(n_neib_send)
    do i = 1, in
      call monolis_bsearch_int(temp, 1, NP, wr(i), id)
      monolis%COM%send_item(i) = local_nid(id)
    enddo

    write(100+monolis_global_myrank(),*)"temp"
    write(100+monolis_global_myrank(),*)temp
    write(100+monolis_global_myrank(),*)"local_nid"
    write(100+monolis_global_myrank(),*)local_nid
    write(100+monolis_global_myrank(),*)"wr"
    write(100+monolis_global_myrank(),*)wr

    write(100+monolis_global_myrank(),*)"monolis%COM%send_item"
    write(100+monolis_global_myrank(),*)monolis%COM%send_item
  end subroutine monolis_com_get_comm_table

end module mod_monolis_util_com
