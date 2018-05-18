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
    integer(kind=kint), pointer :: NeibPE(:)
    integer(kind=kint), pointer :: RecvIndex(:)
    integer(kind=kint), pointer :: RecvItem(:)
    integer(kind=kint), pointer :: SendIndex(:)
    integer(kind=kint), pointer :: SendItem(:)
  end type monolis_com

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
    call MPI_Comm_size(comm, size, ierr)
#endif
  end subroutine monolis_com_size

  subroutine monolis_barrier(comm)
    implicit none
    integer(kind=kint) :: comm
    integer(kind=kint) :: ierr
#ifdef WITHMPI
    call MPI_BARRIER(comm, ierr)
#endif
  end subroutine monolis_barrier

end module mod_monolis_com