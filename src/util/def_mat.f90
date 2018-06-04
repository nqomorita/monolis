module mod_monolis_mat
  use mod_monolis_prm
  implicit none

  type monolis_mat
    integer(kind=kint) :: N, NP, NPU, NPL, NDOF
    integer(kind=kint), pointer :: indexU(:) => null()
    integer(kind=kint), pointer :: indexL(:) => null()
    integer(kind=kint), pointer :: itemU(:) => null()
    integer(kind=kint), pointer :: itemL(:) => null()
    real(kind=kdouble), pointer :: D(:) => null()
    real(kind=kdouble), pointer :: AU(:) => null()
    real(kind=kdouble), pointer :: AL(:) => null()
    real(kind=kdouble), pointer :: X(:) => null()
    real(kind=kdouble), pointer :: B(:) => null()
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