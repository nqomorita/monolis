module mod_monolis_fact_LU
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_fact_LU_33

  implicit none

contains

  subroutine monolis_init_LU_inner(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    if(monoMAT%NDOF == 3)then
      call monolis_init_LU_inner_33(monoPRM, monoCOM, monoMAT)
    endif
  end subroutine monolis_init_LU_inner

  subroutine monolis_fact_LU_inner(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    if(monoMAT%NDOF == 3)then
      call monolis_fact_LU_inner_33(monoPRM, monoCOM, monoMAT)
    endif
  end subroutine monolis_fact_LU_inner

  subroutine monolis_solv_LU_inner(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    if(monoMAT%NDOF == 3)then
      call monolis_solv_LU_inner_33(monoPRM, monoCOM, monoMAT)
    endif
  end subroutine monolis_solv_LU_inner

  subroutine monolis_clear_LU_inner(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    if(monoMAT%NDOF == 3)then
      call monolis_clear_LU_inner_33(monoPRM, monoCOM, monoMAT)
    endif
  end subroutine monolis_clear_LU_inner
end module mod_monolis_fact_LU
