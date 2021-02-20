module mod_monolis_fact_MF_11
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat

  implicit none
  private
  public :: monolis_init_MF_inner_11
  public :: monolis_fact_MF_inner_11
  public :: monolis_solv_MF_inner_11
  public :: monolis_clear_MF_inner_11

contains

  subroutine monolis_init_MF_inner_11(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat_LDU) :: monoMAT

  end subroutine monolis_init_MF_inner_11

  subroutine monolis_clear_MF_inner_11(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat_LDU) :: monoMAT

  end subroutine monolis_clear_MF_inner_11

  subroutine monolis_fact_MF_inner_11(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat_LDU) :: monoMAT

  end subroutine monolis_fact_MF_inner_11

  subroutine monolis_solv_MF_inner_11(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat_LDU) :: monoMAT

  end subroutine monolis_solv_MF_inner_11
end module mod_monolis_fact_MF_11
