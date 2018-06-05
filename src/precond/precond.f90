module mod_monolis_precond
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_precond_diag
  implicit none

contains

  subroutine monolis_precond_setup(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    if(monoPRM%precond == monolis_prec_DIAG)then
      call monolis_precond_diag_setup(monoPRM, monoCOM, monoMAT)
    endif
  end subroutine monolis_precond_setup

  subroutine monolis_precond_apply(monoPRM, monoCOM, monoMAT, X, Y)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: i
    real(kind=kdouble) :: X(:), Y(:)

    if(monoPRM%precond == monolis_prec_DIAG)then
      call monolis_precond_diag_apply(monoPRM, monoCOM, monoMAT, X, Y)
    else
      do i=1, monoMAT%N*monoMAT%NDOF
        Y(i) = X(i)
      enddo
    endif
  end subroutine monolis_precond_apply

  subroutine monolis_precond_clear(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    if(monoPRM%precond == monolis_prec_DIAG)then
      call monolis_precond_diag_clear(monoPRM, monoCOM, monoMAT)
    endif
  end subroutine monolis_precond_clear

end module mod_monolis_precond