module mod_monolis_direct
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  implicit none

contains

  subroutine monolis_solver_LU(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

  end subroutine monolis_solver_LU
end module mod_monolis_direct
