module mod_monolis_precond_ilu
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_precond_ilu_33
  use mod_monolis_precond_ilu_nn
  implicit none

contains

  subroutine monolis_precond_ilu_setup(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    if(monoMAT%NDOF == 3)then
      call monolis_precond_ilu_33_setup(monoMAT)
    else
      call monolis_precond_ilu_nn_setup(monoMAT)
    endif
  end subroutine monolis_precond_ilu_setup

  subroutine monolis_precond_ilu_apply(monoPRM, monoCOM, monoMAT, X, Y)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    real(kind=kdouble) :: X(:), Y(:)

    if(monoMAT%NDOF == 3)then
      call monolis_precond_ilu_33_apply(monoMAT, X, Y)
    else
      call monolis_precond_ilu_nn_apply(monoMAT, X, Y)
    endif
  end subroutine monolis_precond_ilu_apply

  subroutine monolis_precond_ilu_clear(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    if(monoMAT%NDOF == 3)then
      call monolis_precond_ilu_33_clear()
    else
      call monolis_precond_ilu_nn_clear()
    endif
  end subroutine monolis_precond_ilu_clear

end module mod_monolis_precond_ilu