module mod_monolis_util
  implicit none

  integer(4), parameter :: kint    = 4
  integer(4), parameter :: kdouble = 8

  type monolis_matrix
    integer(kind=kint) :: N, NP, NPU, NPL, NDOF
    integer(kind=kint), pointer :: indexU(:)
    integer(kind=kint), pointer :: indexL(:)
    integer(kind=kint), pointer :: itemU(:)
    integer(kind=kint), pointer :: itemL(:)
    real(kind=kdouble), pointer :: D(:)
    real(kind=kdouble), pointer :: AU(:)
    real(kind=kdouble), pointer :: AL(:)m
    real(kind=kdouble), pointer :: X(:)
    real(kind=kdouble), pointer :: B(:)
  end type monolis_matrix

  type(monolis_matrix), save :: monoMAT

contains

  subroutine monolis_mat_initialize(monoMAT)
    implicit none
    type(monolis_matrix) :: monoMAT

  end subroutine monolis_mat_initialize

  subroutine monolis_mat_finalize(monoMAT)
    implicit none
    type(monolis_matrix) :: monoMAT

  end subroutine monolis_mat_finalize
end module mod_monolis_util