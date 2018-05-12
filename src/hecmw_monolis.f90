module mod_monolis_hecmw
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  implicit none

contains

  subroutine solve_hecmw_monolis()
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    write(*,*)"** monolis hello!"

  end subroutine solve_hecmw_monolis

end module mod_monolis_hecmw