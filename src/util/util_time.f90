module mod_monolis_util_time
  use mod_monolis_util
  use mod_monolis_prm
  use mod_monolis_com
  implicit none

contains

  subroutine monolis_timer_initialize(monoPRM, monoCOM)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM

    call monolis_debug_header("monolis_timer_initialize")

    call monolis_barrier_(monoCOM%comm)
    monoPRM%tsol  = monolis_get_time()
    monoPRM%tprep = 0.0d0
    monoPRM%tspmv = 0.0d0
    monoPRM%tdotp = 0.0d0
    monoPRM%tprec = 0.0d0
    monoPRM%tcomm_dotp = 0.0d0
    monoPRM%tcomm_spmv = 0.0d0
  end subroutine monolis_timer_initialize

end module mod_monolis_util_time
