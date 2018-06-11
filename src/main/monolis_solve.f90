module mod_monolis_solve
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_iterative
  use mod_monolis_scaling
  use mod_monolis_precond
  use mod_monolis_util
  implicit none

contains

  subroutine monolis_solve(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    call monolis_timer_initialize()
    call monolis_scaling_fw(monoPRM, monoCOM, monoMAT)
    call monolis_precond_setup(monoPRM, monoCOM, monoMAT)
    call monolis_iterative(monoPRM, monoCOM, monoMAT)
    call monolis_precond_clear(monoPRM, monoCOM, monoMAT)
    call monolis_scaling_bk(monoPRM, monoCOM, monoMAT)
    call monolis_timer_finalize(monoCOM)
  end subroutine monolis_solve

  subroutine monolis_solve_test(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: i, j

    do i = 1, 9
      monoPRM%method = i
      do j = 1, 1
        monoPRM%precond = j
        call monolis_timer_initialize()
        call monolis_scaling_fw(monoPRM, monoCOM, monoMAT)
        call monolis_precond_setup(monoPRM, monoCOM, monoMAT)
        call monolis_iterative(monoPRM, monoCOM, monoMAT)
        call monolis_precond_clear(monoPRM, monoCOM, monoMAT)
        call monolis_scaling_bk(monoPRM, monoCOM, monoMAT)
        call monolis_timer_finalize(monoCOM)
      enddo
    enddo
  end subroutine monolis_solve_test

end module mod_monolis_solve