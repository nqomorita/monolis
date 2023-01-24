module mod_monolis_com
  use mod_monolis_prm
  use iso_c_binding
  implicit none
#ifdef WITH_MPI
  include 'mpif.h'
#endif

  type monolis_com
    integer(kint)          :: myrank
    integer(kint)          :: comm
    integer(kint)          :: commsize
    integer(kint)          :: internal_nnode
    integer(kint)          :: internal_nelem
    logical :: is_overlap = .true.

    integer(kint)          :: recv_n_neib
    integer(kint), pointer :: recv_neib_pe(:) => null()
    integer(kint), pointer :: recv_index(:)   => null()
    integer(kint), pointer :: recv_item(:)    => null()

    integer(kint)          :: send_n_neib
    integer(kint), pointer :: send_neib_pe(:) => null()
    integer(kint), pointer :: send_index(:)   => null()
    integer(kint), pointer :: send_item(:)    => null()

    integer(kint), pointer :: global_node_id(:) => null()
    integer(kint), pointer :: global_elem_id(:) => null()
  end type monolis_com

  integer(kint), parameter :: monolis_sum = 1
  integer(kint), parameter :: monolis_max = 2
  integer(kint), parameter :: monolis_min = 3
#ifdef WITH_MPI
  integer(kint), parameter :: monolis_status_size = MPI_STATUS_SIZE
#else
  integer(kint), parameter :: monolis_status_size = 1
#endif

contains

  subroutine monolis_mpi_initialize()
    implicit none
    integer(kint) :: ierr
#ifdef WITH_MPI
    call MPI_init(ierr)
#endif
  end subroutine monolis_mpi_initialize

  subroutine monolis_mpi_finalize()
    implicit none
    integer(kint) :: ierr
#ifdef WITH_MPI
    call MPI_finalize(ierr)
#endif
  end subroutine monolis_mpi_finalize

  subroutine monolis_com_initialize(monoCOM, is_entire_in)
    implicit none
    type(monolis_com) :: monoCOM
    integer(kint) :: ierr, commsize, myrank
    logical :: is_entire
    logical, optional :: is_entire_in

    is_entire = .false.
    if(present(is_entire_in)) is_entire = is_entire_in

    monoCOM%myrank = 0
    monoCOM%comm = 0
    monoCOM%commsize = 1
    monoCOM%internal_nnode = 0
    monoCOM%internal_nelem = 0

    monoCOM%recv_n_neib = 0
    monoCOM%recv_neib_pe => null()
    monoCOM%recv_index => null()
    monoCOM%recv_item => null()

    monoCOM%send_n_neib = 0
    monoCOM%send_neib_pe => null()
    monoCOM%send_index => null()
    monoCOM%send_item => null()

#ifdef WITH_MPI
    call MPI_comm_size(MPI_COMM_WORLD, commsize, ierr)
    call MPI_comm_rank(MPI_COMM_WORLD, myrank,   ierr)
    monoCOM%comm = MPI_COMM_WORLD
    monoCOM%commsize = commsize
    monoCOM%myrank = myrank

    if(is_entire)then
      monoCOM%comm = MPI_COMM_SELF
      monoCOM%commsize = 1
      monoCOM%myrank = 0
    endif
#endif
  end subroutine monolis_com_initialize

  subroutine monolis_com_finalize(monoCOM)
    implicit none
    type(monolis_com) :: monoCOM

    if(associated(monoCOM%recv_neib_pe)) deallocate(monoCOM%recv_neib_pe)
    if(associated(monoCOM%recv_index)) deallocate(monoCOM%recv_index)
    if(associated(monoCOM%recv_item)) deallocate(monoCOM%recv_item)

    if(associated(monoCOM%send_neib_pe)) deallocate(monoCOM%send_neib_pe)
    if(associated(monoCOM%send_index)) deallocate(monoCOM%send_index)
    if(associated(monoCOM%send_item)) deallocate(monoCOM%send_item)

    if(associated(monoCOM%global_node_id)) deallocate(monoCOM%global_node_id)
    if(associated(monoCOM%global_elem_id)) deallocate(monoCOM%global_elem_id)

    monoCOM%recv_neib_pe => null()
    monoCOM%recv_index => null()
    monoCOM%recv_item => null()

    monoCOM%send_neib_pe => null()
    monoCOM%send_index => null()
    monoCOM%send_item => null()

    monoCOM%global_node_id => null()
    monoCOM%global_elem_id => null()
  end subroutine monolis_com_finalize

  subroutine monolis_com_input_comm_table(monoCOM, fname_in)
    implicit none
    type(monolis_com) :: monoCOM
    integer(kint) :: i, j, nitem
    character :: cnum*5, header*128
    character(*) :: fname_in

    if(monoCOM%commsize <= 1)then
      monoCOM%commsize = 1
      return
    endif

    write(cnum,"(i0)")monoCOM%myrank
    header = trim(fname_in)//"/parted.0/"

    open(10, file=trim(header)//"monolis.send."//trim(cnum), status='old')
      !> for overlap
      read(10,*) monoCOM%send_n_neib, nitem
      allocate(monoCOM%send_neib_pe(monoCOM%send_n_neib))
      do i = 1, monoCOM%send_n_neib
        read(10,*) monoCOM%send_neib_pe(i)
      enddo
      allocate(monoCOM%send_index(0:monoCOM%send_n_neib), source = 0)
      allocate(monoCOM%send_item(nitem), source = 0)

      if(monoCOM%send_n_neib == 0) return

      do i = 0, monoCOM%send_n_neib
        read(10,*) monoCOM%send_index(i)
      enddo
      do i = 1, nitem
        read(10,*) monoCOM%send_item(i)
      enddo
    close(10)

    open(10, file=trim(header)//"monolis.recv."//trim(cnum), status='old')
      read(10,*)monoCOM%recv_n_neib, nitem
      if(monoCOM%send_n_neib /= monoCOM%recv_n_neib)then
        stop "** error: monolis_com_input_comm_table"
      endif

      allocate(monoCOM%recv_neib_pe(monoCOM%recv_n_neib))
      do i = 1, monoCOM%recv_n_neib
        read(10,*) monoCOM%recv_neib_pe(i)
      enddo
      allocate(monoCOM%recv_index(0:monoCOM%recv_n_neib), source = 0)
      allocate(monoCOM%recv_item(nitem), source = 0)

      do i = 0, monoCOM%recv_n_neib
        read(10,*) monoCOM%recv_index(i)
      enddo
      do i = 1, nitem
        read(10,*) monoCOM%recv_item(i)
      enddo
    close(10)

    open(10, file=trim(header)//"node.id."//trim(cnum), status='old')
      read(10,*)nitem
      allocate(monoCOM%global_node_id(nitem), source = 0)
      do i = 1, nitem
        read(10,*) j, j, monoCOM%global_node_id(i)
      enddo
    close(10)

    open(10, file=trim(header)//"connectivity.id."//trim(cnum), status='old')
      read(10,*)nitem
      allocate(monoCOM%global_elem_id(nitem), source = 0)
      do i = 1, nitem
        read(10,*) j, j, monoCOM%global_elem_id(i)
      enddo
    close(10)

    open(10, file=trim(header)//"node.n_internal."//trim(cnum), status='old')
      read(10,*)monoCOM%internal_nnode
    close(10)

    open(10, file=trim(header)//"connectivity.n_internal."//trim(cnum), status='old')
      read(10,*)monoCOM%internal_nelem
    close(10)
  end subroutine monolis_com_input_comm_table

  subroutine monolis_com_copy(monoCOM, monoCOM_reorder)
    implicit none
    type(monolis_com) :: monoCOM
    type(monolis_com) :: monoCOM_reorder

    monoCOM_reorder%myrank = monoCOM%myrank
    monoCOM_reorder%comm = monoCOM%comm
    monoCOM_reorder%commsize = monoCOM%commsize
    monoCOM_reorder%internal_nnode = monoCOM%internal_nnode

    monoCOM_reorder%recv_n_neib = monoCOM%recv_n_neib
    monoCOM_reorder%recv_neib_pe => monoCOM%recv_neib_pe
    monoCOM_reorder%recv_index => monoCOM%recv_index
    monoCOM_reorder%recv_item => monoCOM%recv_item

    monoCOM_reorder%send_n_neib = monoCOM%send_n_neib
    monoCOM_reorder%send_neib_pe => monoCOM%send_neib_pe
    monoCOM_reorder%send_index => monoCOM%send_index
    monoCOM_reorder%send_item => monoCOM%send_item

    monoCOM_reorder%global_node_id => monoCOM%global_node_id
    monoCOM_reorder%global_elem_id => monoCOM%global_elem_id
  end subroutine monolis_com_copy

  function monolis_global_comm()
    implicit none
    integer(kint) :: monolis_global_comm

#ifdef WITH_MPI
    monolis_global_comm = MPI_COMM_WORLD
#else
    monolis_global_comm = 0
#endif
  end function monolis_global_comm

  function monolis_global_commsize()
    implicit none
    integer(kint) :: monolis_global_commsize, ierr

#ifdef WITH_MPI
    call MPI_comm_size(MPI_COMM_WORLD, monolis_global_commsize, ierr)
#else
    monolis_global_commsize = 1
#endif
  end function monolis_global_commsize

  function monolis_global_myrank()
    implicit none
    integer(kint) :: monolis_global_myrank, ierr

#ifdef WITH_MPI
    call MPI_comm_rank(MPI_COMM_WORLD, monolis_global_myrank, ierr)
#else
    monolis_global_myrank = 0
#endif
  end function monolis_global_myrank

  subroutine monolis_barrier_(comm)
    implicit none
    integer(kint) :: comm
    integer(kint) :: ierr
#ifdef WITH_MPI
    call MPI_barrier(comm, ierr)
#endif
  end subroutine monolis_barrier_

end module mod_monolis_com