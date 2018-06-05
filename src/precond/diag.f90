module mod_monolis_precond_diag
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_precond_diag_33
  implicit none

contains

  subroutine monolis_precond_diag_setup(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    if(monoMAT%NDOF == 3)then
      call monolis_precond_diag_33_setup(monoPRM, monoCOM, monoMAT)
    endif
  end subroutine monolis_precond_diag_setup

  subroutine monolis_precond_diag_apply(monoPRM, monoCOM, monoMAT, X, Y)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: i
    real(kind=kdouble) :: X(:), Y(:)

    if(monoMAT%NDOF == 3)then
      call monolis_precond_diag_33_apply(monoPRM, monoCOM, monoMAT, X, Y)
    endif
  end subroutine monolis_precond_diag_apply

  subroutine monolis_precond_diag_clear(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    if(monoMAT%NDOF == 3)then
      call monolis_precond_diag_33_clear(monoPRM, monoCOM, monoMAT)
    endif
  end subroutine monolis_precond_diag_clear

end module mod_monolis_precond_diag