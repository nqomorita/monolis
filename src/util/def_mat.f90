module mod_monolis_mat
  use mod_monolis_prm
  implicit none

  type monolis_mat
    integer(kind=kint) :: N, NP, NPU, NPL, NDOF
    integer(kind=kint), pointer :: indexU(:)
    integer(kind=kint), pointer :: indexL(:)
    integer(kind=kint), pointer :: itemU(:)
    integer(kind=kint), pointer :: itemL(:)
    real(kind=kdouble), pointer :: D(:)
    real(kind=kdouble), pointer :: AU(:)
    real(kind=kdouble), pointer :: AL(:)
    real(kind=kdouble), pointer :: X(:)
    real(kind=kdouble), pointer :: B(:)
  end type monolis_mat

contains

  subroutine monolis_mat_initialize(monoMAT)
    implicit none
    type(monolis_mat) :: monoMAT

  end subroutine monolis_mat_initialize

  subroutine monolis_mat_finalize(monoMAT)
    implicit none
    type(monolis_mat) :: monoMAT

  end subroutine monolis_mat_finalize
end module mod_monolis_mat