module mod_monolis_precond_mumps
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat

  implicit none

contains

  subroutine monolis_precond_mumps_setup(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

  end subroutine monolis_precond_mumps_setup

  subroutine monolis_precond_mumps_apply(monoPRM, monoCOM, monoMAT, X, Y)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    real(kind=kdouble) :: X(:), Y(:)

  end subroutine monolis_precond_mumps_apply

  subroutine monolis_precond_mumps_clear(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

  end subroutine monolis_precond_mumps_clear
end module mod_monolis_precond_mumps