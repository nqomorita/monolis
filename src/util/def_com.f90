module mod_monolis_com
  use mod_monolis_prm
  use iso_c_binding
  implicit none
#ifdef WITH_MPI
  include 'mpif.h'
#endif

  type monolis_com
    integer(kind=kint)          :: myrank
    integer(kind=kint)          :: comm
    integer(kind=kint)          :: commsize
    integer(kind=kint)          :: n_neib
    integer(kind=kint), pointer :: neib_pe(:)    => null()
    integer(kind=kint), pointer :: recv_index(:) => null()
    integer(kind=kint), pointer :: recv_item(:)  => null()
    integer(kind=kint), pointer :: send_index(:) => null()
    integer(kind=kint), pointer :: send_item(:)  => null()
  end type monolis_com

  type, bind(c) :: monolis_com_c
    integer(c_int) :: myrank
    integer(c_int) :: comm
    integer(c_int) :: commsize
    integer(c_int) :: n_neib
    type(c_ptr) :: neib_pe
    type(c_ptr) :: recv_index
    type(c_ptr) :: recv_item
    type(c_ptr) :: send_index
    type(c_ptr) :: send_item
  end type monolis_com_c

  integer(kind=kint), parameter :: monolis_sum = 1
  integer(kind=kint), parameter :: monolis_max = 2
  integer(kind=kint), parameter :: monolis_min = 3
#ifdef WITH_MPI
  integer(kind=kint), parameter :: monolis_status_size = MPI_STATUS_SIZE
#else
  integer(kind=kint), parameter :: monolis_status_size = 1
#endif

contains

  subroutine monolis_com_initialize(monoCOM)
    implicit none
    type(monolis_com) :: monoCOM
    integer(kind=kint) :: ierr, commsize, myrank

    monoCOM%myrank = 0
    monoCOM%comm = 0
    monoCOM%commsize = 0
    monoCOM%n_neib = 0
    monoCOM%neib_pe => null()
    monoCOM%recv_index => null()
    monoCOM%recv_item => null()
    monoCOM%send_index => null()
    monoCOM%send_item => null()

#ifdef WITH_MPI
    call MPI_init(ierr)
    call MPI_comm_size(MPI_COMM_WORLD, commsize, ierr)
    call MPI_comm_rank(MPI_COMM_WORLD, myrank,   ierr)
    monoCOM%comm = MPI_COMM_WORLD
    monoCOM%commsize = commsize
    monoCOM%myrank = myrank
#endif
  end subroutine monolis_com_initialize

  subroutine monolis_com_initialize_c(monoCOM_c) bind(c, name="monolis_com_initialize")
    use iso_c_binding
    implicit none
    type(monolis_com_c) :: monoCOM_c
    integer(kind=kint) :: ierr, commsize, myrank

    monoCOM_c%myrank = 0
    monoCOM_c%comm = 0
    monoCOM_c%commsize = 0
    monoCOM_c%n_neib = 0

#ifdef WITH_MPI
    call MPI_init(ierr)
    call MPI_comm_size(MPI_COMM_WORLD, commsize, ierr)
    call MPI_comm_rank(MPI_COMM_WORLD, myrank,   ierr)
    monoCOM_c%comm = MPI_COMM_WORLD
    monoCOM_c%commsize = commsize
    monoCOM_c%myrank = myrank
#endif
  end subroutine monolis_com_initialize_c

  subroutine monolis_com_finalize(monoCOM)
    implicit none
    type(monolis_com) :: monoCOM
    integer(kind=kint) :: ierr

    if(associated(monoCOM%neib_pe)) deallocate(monoCOM%neib_pe)
    if(associated(monoCOM%recv_index)) deallocate(monoCOM%recv_index)
    if(associated(monoCOM%recv_item)) deallocate(monoCOM%recv_item)
    if(associated(monoCOM%send_index)) deallocate(monoCOM%send_index)
    if(associated(monoCOM%send_item)) deallocate(monoCOM%send_item)
    monoCOM%neib_pe => null()
    monoCOM%recv_index => null()
    monoCOM%recv_item => null()
    monoCOM%send_index => null()
    monoCOM%send_item => null()

#ifdef WITH_MPI
    call MPI_finalize(ierr)
#endif
  end subroutine monolis_com_finalize

  subroutine monolis_com_finalize_c(monoCOM_c) bind(c, name="monolis_com_finalize")
    use iso_c_binding
    implicit none
    type(monolis_com_c) :: monoCOM_c
    integer(kind=kint) :: ierr

#ifdef WITH_MPI
    call MPI_finalize(ierr)
#endif
  end subroutine monolis_com_finalize_c

  subroutine monolis_com_copy(monoCOM, monoCOM_reorder)
    implicit none
    type(monolis_com) :: monoCOM
    type(monolis_com) :: monoCOM_reorder

    monoCOM_reorder%myrank = monoCOM%myrank
    monoCOM_reorder%comm = monoCOM%comm
    monoCOM_reorder%commsize = monoCOM%commsize
    monoCOM_reorder%n_neib = monoCOM%n_neib
    monoCOM_reorder%neib_pe => monoCOM%neib_pe
    monoCOM_reorder%recv_index => monoCOM%recv_index
    monoCOM_reorder%recv_item => monoCOM%recv_item
    monoCOM_reorder%send_index => monoCOM%send_index
    monoCOM_reorder%send_item => monoCOM%send_item
  end subroutine monolis_com_copy

  subroutine monolis_com_size(size, comm)
    implicit none
    integer(kind=kint) :: size, comm
    integer(kind=kint) :: ierr

    size = 1
#ifdef WITH_MPI
    call MPI_comm_size(comm, size, ierr)
#endif
  end subroutine monolis_com_size

  subroutine monolis_barrier(comm)
    implicit none
    integer(kind=kint) :: comm
    integer(kind=kint) :: ierr
#ifdef WITH_MPI
    call MPI_barrier(comm, ierr)
#endif
  end subroutine monolis_barrier

  function monolis_wtime()
    implicit none
#ifdef WITH_MPI
    real(kind=kdouble) :: monolis_wtime
    monolis_wtime = MPI_wtime()
#else
    integer(kind=kint) :: monolis_wtime
    call system_clock(monolis_wtime)
#endif
  end function monolis_wtime

end module mod_monolis_com