module mod_monolis_iterative
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_solver_CG
  implicit none

contains

  subroutine monolis_iterative(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    select case(monoPRM%method)
      case (1)
        call monolis_solver_CG(monoPRM, monoCOM, monoMAT)
    end select

  end subroutine monolis_iterative

end module mod_monolis_iterative