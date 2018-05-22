module mod_monolis_prm
  implicit none

  integer(4), parameter :: kint    = 4
  integer(4), parameter :: kdouble = 8

  type monolis_prm
    integer(kind=kint) :: method
    integer(kind=kint) :: precond
    integer(kind=kint) :: maxiter
  end type monolis_prm

contains

  subroutine monolis_prm_initialize(monoPRM)
    implicit none
    type(monolis_prm) :: monoPRM

  end subroutine monolis_prm_initialize

  subroutine monolis_prm_finalize(monoPRM)
    implicit none
    type(monolis_prm) :: monoPRM

  end subroutine monolis_prm_finalize

end module mod_monolis_prm