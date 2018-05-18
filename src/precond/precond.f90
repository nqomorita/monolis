module mod_monolis_precond
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  implicit none

contains

  subroutine monolis_precond(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    call monolis_prm_initialize(monoPRM)
    call monolis_com_initialize(monoCOM)
    call monolis_mat_initialize(monoMAT)
  end subroutine monolis_precond

end module mod_monolis_precond