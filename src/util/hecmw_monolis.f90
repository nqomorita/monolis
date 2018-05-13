module mod_monolis_hecmw
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  implicit none

contains

  subroutine monolis_solve_hecmw_inner()
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    write(*,*)"** monolis hello!"

  end subroutine monolis_solve_hecmw_inner

end module mod_monolis_hecmw