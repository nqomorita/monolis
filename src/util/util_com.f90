module mod_monolis_util_com
  use mod_monolis_util
  use mod_monolis_prm
  use mod_monolis_com
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

end module mod_monolis_util_com
