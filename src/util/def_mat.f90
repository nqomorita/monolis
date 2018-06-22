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

  subroutine monolis_mat_nullify(monoMAT)
    implicit none
    type(monolis_mat) :: monoMAT

    monoMAT%N = 0
    monoMAT%NP = 0
    monoMAT%NPU = 0
    monoMAT%NPL = 0
    monoMAT%NDOF = 0
    monoMAT%indexU => null()
    monoMAT%indexL => null()
    monoMAT%itemU => null()
    monoMAT%itemL => null()
    monoMAT%D => null()
    monoMAT%AU => null()
    monoMAT%AL => null()
    monoMAT%X => null()
    monoMAT%B => null()
  end subroutine monolis_mat_nullify

  subroutine monolis_mat_initialize(monoMAT)
    implicit none
    type(monolis_mat) :: monoMAT

    monoMAT%N = 0
    monoMAT%NP = 0
    monoMAT%NPU = 0
    monoMAT%NPL = 0
    monoMAT%NDOF = 0
    monoMAT%indexU => null()
    monoMAT%indexL => null()
    monoMAT%itemU => null()
    monoMAT%itemL => null()
    monoMAT%D => null()
    monoMAT%AU => null()
    monoMAT%AL => null()
    monoMAT%X => null()
    monoMAT%B => null()
  end subroutine monolis_mat_initialize

  subroutine monolis_mat_copy(monoA, monoB)
    implicit none
    type(monolis_mat) :: monoA
    type(monolis_mat) :: monoB

    monoB%N = monoA%N
    monoB%NP = monoA%NP
    monoB%NPU = monoA%NPU
    monoB%NPL = monoA%NPL
    monoB%NDOF = monoA%NDOF
    monoB%indexU => monoA%indexU
    monoB%indexL => monoA%indexL
    monoB%itemU => monoA%itemU
    monoB%itemL => monoA%itemL
    monoB%D => monoA%D
    monoB%AU => monoA%AU
    monoB%AL => monoA%AL
    monoB%X => monoA%X
    monoB%B => monoA%B
  end subroutine monolis_mat_copy

  subroutine monolis_mat_finalize(monoMAT)
    implicit none
    type(monolis_mat) :: monoMAT

    if(associated(monoMAT%indexU)) deallocate(monoMAT%indexU)
    if(associated(monoMAT%indexL)) deallocate(monoMAT%indexL)
    if(associated(monoMAT%itemU)) deallocate(monoMAT%itemU)
    if(associated(monoMAT%itemL)) deallocate(monoMAT%itemL)
    if(associated(monoMAT%D)) deallocate(monoMAT%D)
    if(associated(monoMAT%AU)) deallocate(monoMAT%AU)
    if(associated(monoMAT%AL)) deallocate(monoMAT%AL)
    if(associated(monoMAT%X)) deallocate(monoMAT%X)
    if(associated(monoMAT%B)) deallocate(monoMAT%B)
    monoMAT%indexU => null()
    monoMAT%indexL => null()
    monoMAT%itemU => null()
    monoMAT%itemL => null()
    monoMAT%D => null()
    monoMAT%AU => null()
    monoMAT%AL => null()
    monoMAT%X => null()
    monoMAT%B => null()
  end subroutine monolis_mat_finalize
end module mod_monolis_mat