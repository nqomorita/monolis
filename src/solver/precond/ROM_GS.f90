module mod_monolis_precond_rom_gs
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_precond_rom_gs_33
  !use mod_monolis_precond_rom_gs_nn

  implicit none

contains

  subroutine monolis_precond_rom_gs_setup(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    if(monoMAT%NDOF == 3)then
      call monolis_precond_rom_gs_33_setup(monoMAT)
    else
      stop
      !call monolis_precond_rom_gs_nn_setup(monoMAT)
    endif
  end subroutine monolis_precond_rom_gs_setup

  subroutine monolis_precond_rom_gs_apply(monoPRM, monoCOM, monoMAT, X, Y, iter)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: iter
    real(kdouble) :: X(:), Y(:)

    if(monoMAT%NDOF == 3)then
      call monolis_precond_rom_gs_33_apply(monoMAT, X, Y, iter)
    else
      stop
      !call monolis_precond_rom_gs_nn_apply(monoMAT, X, Y)
    endif
  end subroutine monolis_precond_rom_gs_apply

  subroutine monolis_precond_rom_gs_clear(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    if(monoMAT%NDOF == 3)then
      call monolis_precond_rom_gs_33_clear(monoMAT)
      stop
    else
      !call monolis_precond_rom_gs_nn_clear(monoMAT)
    endif
  end subroutine monolis_precond_rom_gs_clear
end module mod_monolis_precond_rom_gs