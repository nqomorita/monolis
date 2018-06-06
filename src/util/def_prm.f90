module mod_monolis_prm
  implicit none

  integer(4), parameter :: kint    = 4
  integer(4), parameter :: kdouble = 8

  integer(4), parameter :: monolis_iter_CG       = 1
  integer(4), parameter :: monolis_iter_BiCGSTAB = 2
  integer(4), parameter :: monolis_iter_BiCGSTAB_noprec = 3
  integer(4), parameter :: monolis_iter_GropCG   = 14
  integer(4), parameter :: monolis_iter_PipeCG   = 4

  integer(4), parameter :: monolis_prec_DIAG   = 1
  integer(4), parameter :: monolis_prec_ILU    = 2
  integer(4), parameter :: monolis_prec_SOR    = 3
  integer(4), parameter :: monolis_prec_JACOBI = 4
  integer(4), parameter :: monolis_prec_SAINV  = 5
  integer(4), parameter :: monolis_prec_RIF    = 6
  integer(4), parameter :: monolis_prec_SPIKE  = 7

  type monolis_prm
    integer(kind=kint) :: method
    integer(kind=kint) :: precond
    integer(kind=kint) :: maxiter
    real(kind=kdouble) :: tol
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