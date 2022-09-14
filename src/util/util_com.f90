module mod_monolis_util_com
  use mod_monolis_util
  use mod_monolis_prm
  use mod_monolis_com
  implicit none

contains

  !> set parameter section
  subroutine monolis_com_set_communicator(monolis, comm)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: comm
    monolis%COM%comm = comm
  end subroutine monolis_com_set_communicator

end module mod_monolis_util_com
