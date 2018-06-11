module mod_monolis_util
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  implicit none

  private
  public :: monolis_initialize
  public :: monolis_finalize
  public :: monolis_timer_initialize
  public :: monolis_timer_finalize

  !> tsol = tspmv + tprec + tcomm + others
  real(kind=kdouble) :: tsol  = 0.0d0
  real(kind=kdouble) :: tspmv = 0.0d0
  real(kind=kdouble) :: tprec = 0.0d0
  real(kind=kdouble) :: tcomm = 0.0d0

contains

  subroutine monolis_initialize(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    call monolis_prm_initialize(monoPRM)
    call monolis_com_initialize(monoCOM)
    call monolis_mat_initialize(monoMAT)
  end subroutine monolis_initialize

  subroutine monolis_finalize(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    call monolis_prm_finalize(monoPRM)
    call monolis_com_finalize(monoCOM)
    call monolis_mat_finalize(monoMAT)
  end subroutine monolis_finalize

  subroutine monolis_timer_initialize()
    implicit none

    tsol  = 0.0d0
    tspmv = 0.0d0
    tprec = 0.0d0
    tcomm = 0.0d0
  end subroutine monolis_timer_initialize

  subroutine monolis_timer_finalize(monoCOM)
    implicit none
    type(monolis_com) :: monoCOM

    if(monoCOM%myrank == 0) write(*,"(a,i8,1p4e12.5)")" ** monolis solved:", 0, tsol, tspmv, tprec, tcomm
  end subroutine monolis_timer_finalize
end module mod_monolis_util