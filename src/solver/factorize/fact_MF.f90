module mod_monolis_fact_MF
  use mod_monolis_prm
  use mod_monolis_mat
  use mod_monolis_fact_MF_11

  implicit none

contains

  subroutine monolis_init_MF_inner(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat_LDU) :: monoMAT

!    if(monoMAT%NDOF == 3)then
!      call monolis_init_MF_inner_33(monoPRM, monoCOM, monoMAT)
!    elseif(monoMAT%NDOF == 1)then
      call monolis_init_MF_inner_11(monoPRM, monoCOM, monoMAT)
!    else
!      call monolis_init_MF_inner_nn(monoPRM, monoCOM, monoMAT)
!    endif
  end subroutine monolis_init_MF_inner

  subroutine monolis_fact_MF_inner(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat_LDU) :: monoMAT

!    if(monoMAT%NDOF == 3)then
!      call monolis_fact_MF_inner_33(monoPRM, monoCOM, monoMAT)
!    elseif(monoMAT%NDOF == 1)then
      call monolis_fact_MF_inner_11(monoPRM, monoCOM, monoMAT)
!    else
!      call monolis_fact_MF_inner_nn(monoPRM, monoCOM, monoMAT)
!    endif
  end subroutine monolis_fact_MF_inner

  subroutine monolis_solv_MF_inner(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat_LDU) :: monoMAT

!    if(monoMAT%NDOF == 3)then
!      call monolis_solv_MF_inner_33(monoPRM, monoCOM, monoMAT)
!    elseif(monoMAT%NDOF == 1)then
      call monolis_solv_MF_inner_11(monoPRM, monoCOM, monoMAT)
!    else
!      call monolis_solv_MF_inner_nn(monoPRM, monoCOM, monoMAT)
!    endif
  end subroutine monolis_solv_MF_inner

  subroutine monolis_clear_MF_inner(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat_LDU) :: monoMAT

!    if(monoMAT%NDOF == 3)then
!      call monolis_clear_MF_inner_33(monoPRM, monoCOM, monoMAT)
!    elseif(monoMAT%NDOF == 1)then
      call monolis_clear_MF_inner_11(monoPRM, monoCOM, monoMAT)
!    else
!      call monolis_clear_MF_inner_nn(monoPRM, monoCOM, monoMAT)
!    endif
  end subroutine monolis_clear_MF_inner
end module mod_monolis_fact_MF
