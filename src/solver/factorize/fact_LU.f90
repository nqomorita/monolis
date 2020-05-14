module mod_monolis_fact_LU
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_fact_LU_11
  use mod_monolis_fact_LU_33
  use mod_monolis_fact_LU_nn

  implicit none

contains

  subroutine monolis_init_LU_inner(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat_LDU) :: monoMAT

    if(monoMAT%NDOF == 3)then
      !call monolis_init_LU_inner_33(monoPRM, monoCOM, monoMAT)
      call monolis_init_LU_inner_nn(monoPRM, monoCOM, monoMAT)
    elseif(monoMAT%NDOF == 1)then
      call monolis_init_LU_inner_11(monoPRM, monoCOM, monoMAT)
    else
      call monolis_init_LU_inner_nn(monoPRM, monoCOM, monoMAT)
    endif
  end subroutine monolis_init_LU_inner

  subroutine monolis_fact_LU_inner(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat_LDU) :: monoMAT

    if(monoMAT%NDOF == 3)then
      !call monolis_fact_LU_inner_33(monoPRM, monoCOM, monoMAT)
      call monolis_fact_LU_inner_nn(monoPRM, monoCOM, monoMAT)
    elseif(monoMAT%NDOF == 1)then
      call monolis_fact_LU_inner_11(monoPRM, monoCOM, monoMAT)
    else
      call monolis_fact_LU_inner_nn(monoPRM, monoCOM, monoMAT)
    endif
  end subroutine monolis_fact_LU_inner

  subroutine monolis_solv_LU_inner(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat_LDU) :: monoMAT

    if(monoMAT%NDOF == 3)then
      !call monolis_solv_LU_inner_33(monoPRM, monoCOM, monoMAT)
      call monolis_solv_LU_inner_nn(monoPRM, monoCOM, monoMAT)
    elseif(monoMAT%NDOF == 1)then
      call monolis_solv_LU_inner_11(monoPRM, monoCOM, monoMAT)
    else
      call monolis_solv_LU_inner_nn(monoPRM, monoCOM, monoMAT)
    endif
  end subroutine monolis_solv_LU_inner

  subroutine monolis_clear_LU_inner(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat_LDU) :: monoMAT

    if(monoMAT%NDOF == 3)then
      !call monolis_clear_LU_inner_33(monoPRM, monoCOM, monoMAT)
      call monolis_clear_LU_inner_nn(monoPRM, monoCOM, monoMAT)
    elseif(monoMAT%NDOF == 1)then
      call monolis_clear_LU_inner_11(monoPRM, monoCOM, monoMAT)
    else
      call monolis_clear_LU_inner_nn(monoPRM, monoCOM, monoMAT)
    endif
  end subroutine monolis_clear_LU_inner
end module mod_monolis_fact_LU
