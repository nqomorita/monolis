module mod_monolis_solve
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_iterative
  implicit none

contains

  subroutine monolis_solve(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    call monolis_iterative(monoPRM, monoCOM, monoMAT)

  end subroutine monolis_solve

end module mod_monolis_solve