module mod_monolis_comm
  use mod_monolis_util
  implicit none

  type monolis_comm
    integer(kind=kint)          :: myrank
    integer(kind=kint)          :: comm
    integer(kind=kint)          :: commsize
    integer(kind=kint), pointer :: NeibPE(:)
    integer(kind=kint), pointer :: RecvIndex(:)
    integer(kind=kint), pointer :: RecvItem(:)
    integer(kind=kint), pointer :: SendIndex(:)
    integer(kind=kint), pointer :: SendItem(:)
  end type monolis_comm

  type(monolis_comm), save :: monoCOM

contains

  subroutine monolis_comm_initialize(monoCOM)
    implicit none
    type(monolis_comm) :: monoCOM

  end subroutine monolis_comm_initialize

end module mod_monolis_comm