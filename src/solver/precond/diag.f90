module mod_monolis_precond_diag
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_precond_diag_33
  use mod_monolis_precond_diag_nn

  implicit none

contains

  subroutine monolis_precond_diag_setup(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    if(monoMAT%NDOF == 3)then
      call monolis_precond_diag_33_setup(monoMAT)
    else
      call monolis_precond_diag_nn_setup(monoMAT)
    endif
  end subroutine monolis_precond_diag_setup

  subroutine monolis_precond_diag_apply(monoPRM, monoCOM, monoMAT, X, Y)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    real(kind=kdouble) :: X(:), Y(:)

    if(monoMAT%NDOF == 3)then
      call monolis_precond_diag_33_apply(monoMAT, X, Y)
    else
      call monolis_precond_diag_nn_apply(monoMAT, X, Y)
    endif
  end subroutine monolis_precond_diag_apply

  subroutine monolis_precond_diag_clear(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    if(monoMAT%NDOF == 3)then
      call monolis_precond_diag_33_clear(monoMAT)
    else
      call monolis_precond_diag_nn_clear(monoMAT)
    endif
  end subroutine monolis_precond_diag_clear
end module mod_monolis_precond_diag