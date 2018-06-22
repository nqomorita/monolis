module mod_monolis_precond_ilu
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat

  implicit none

contains

  subroutine  monolis_precond_ilu_setup(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    !call monolis_solver_direct_init(monoPRM, monoCOM, monoMAT)
    !call monolis_solver_direct_fact(monoPRM, monoCOM, monoMAT)
  end subroutine monolis_precond_ilu_setup

  subroutine monolis_precond_ilu_apply(monoPRM, monoCOM, monoMAT, X, Y)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    real(kind=kdouble) :: X(:), Y(:)

    !call monolis_solver_direct_solv(monoPRM, monoCOM, monoMAT)
  end subroutine monolis_precond_ilu_apply

  subroutine monolis_precond_ilu_clear(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    !call monolis_solver_direct_clear(monoPRM, monoCOM, monoMAT)
  end subroutine monolis_precond_ilu_clear

!  subroutine monolis_solver_direct_init(monoPRM, monoCOM, monoMAT)
!    implicit none
!    type(monolis_prm) :: monoPRM
!    type(monolis_com) :: monoCOM
!    type(monolis_mat) :: monoMAT
!
!    if(monoCOM%commsize == 0) isEntire = .true.
!    call monolis_matrix_get_fillin(hecMESH, hecT, idxU, itemU, NPU)
!    call hecmw_matrix_copy_with_fillin(hecMESH, hecT, idxU, itemU, AU, NPU)
!  end subroutine monolis_solver_direct_init
end module mod_monolis_precond_ilu