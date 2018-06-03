module mod_monolis_com
  use mod_monolis_prm
  implicit none
#ifdef WITHMPI
  include 'mpif.h'
#endif

  type monolis_com
    integer(kind=kint)          :: myrank
    integer(kind=kint)          :: comm
    integer(kind=kint)          :: commsize
    integer(kind=kint)          :: n_neib
    integer(kind=kint), pointer :: neib_pe(:)
    integer(kind=kint), pointer :: recv_index(:)
    integer(kind=kint), pointer :: recv_item(:)
    integer(kind=kint), pointer :: send_index(:)
    integer(kind=kint), pointer :: send_item(:)
  end type monolis_com

  integer(kind=kint), parameter :: monolis_sum = 1
  integer(kind=kint), parameter :: monolis_max = 2
  integer(kind=kint), parameter :: monolis_min = 3

  integer(kind=kint), parameter :: monolis_status_size = MPI_STATUS_SIZE

contains

  subroutine monolis_com_initialize(monoCOM)
    implicit none
    type(monolis_com) :: monoCOM

  end subroutine monolis_com_initialize

  subroutine monolis_com_finalize(monoCOM)
    implicit none
    type(monolis_com) :: monoCOM

  end subroutine monolis_com_finalize

  subroutine monolis_com_size(size, comm)
    implicit none
    integer(kind=kint) :: size, comm
    integer(kind=kint) :: ierr

    size = 1
#ifdef WITHMPI
    call MPI_comm_size(comm, size, ierr)
#endif
  end subroutine monolis_com_size

  subroutine monolis_barrier(comm)
    implicit none
    integer(kind=kint) :: comm
    integer(kind=kint) :: ierr
#ifdef WITHMPI
    call MPI_barrier(comm, ierr)
#endif
  end subroutine monolis_barrier

  function monolis_wtime()
    implicit none
    real(kind=kdouble) :: monolis_wtime
#ifdef WITHMPI
    monolis_wtime = MPI_wtime()
#else
    call system_clock(monolis_wtime)
#endif
  end function monolis_wtime

end module mod_monolis_com