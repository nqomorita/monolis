module mod_monolis_solve
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_iterative
  use mod_monolis_scaling
  use mod_monolis_precond
  implicit none

contains

  subroutine monolis_solve(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    call monolis_scaling_fw(monoPRM, monoCOM, monoMAT)

    call monolis_precond_setup(monoPRM, monoCOM, monoMAT)

    call monolis_iterative(monoPRM, monoCOM, monoMAT)

    call monolis_precond_clear(monoPRM, monoCOM, monoMAT)

    call monolis_scaling_bk(monoPRM, monoCOM, monoMAT)

  end subroutine monolis_solve

end module mod_monolis_solve