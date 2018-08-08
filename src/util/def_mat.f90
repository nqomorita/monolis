module mod_monolis_mat
  use mod_monolis_prm
  implicit none

  type monolis_mat
    integer(kind=kint) :: N, NP, NZ, NDOF
    integer(kind=kint), pointer :: index(:) => null()
    integer(kind=kint), pointer :: item(:) => null()
    integer(kind=kint), pointer :: perm(:) => null()
    integer(kind=kint), pointer :: iperm(:) => null()
    real(kind=kdouble), pointer :: A(:) => null()
    real(kind=kdouble), pointer :: X(:) => null()
    real(kind=kdouble), pointer :: B(:) => null()
    real(kind=kdouble), pointer :: diag(:) => null()
  end type monolis_mat

  type monolis_mat_LDU
    integer(kind=kint) :: N, NP, NPU, NPL, NDOF
    integer(kind=kint), pointer :: indexU(:) => null()
    integer(kind=kint), pointer :: itemU(:) => null()
    integer(kind=kint), pointer :: indexL(:) => null()
    integer(kind=kint), pointer :: itemL(:) => null()
    real(kind=kdouble), pointer :: U(:) => null()
    real(kind=kdouble), pointer :: D(:) => null()
    real(kind=kdouble), pointer :: L(:) => null()
    real(kind=kdouble), pointer :: X(:) => null()
    real(kind=kdouble), pointer :: B(:) => null()
  end type monolis_mat_LDU

contains

  subroutine monolis_mat_nullify(monoMAT)
    implicit none
    type(monolis_mat) :: monoMAT

    monoMAT%N = 0
    monoMAT%NP = 0
    monoMAT%NZ = 0
    monoMAT%NDOF = 0
    monoMAT%index => null()
    monoMAT%item => null()
    monoMAT%A => null()
    monoMAT%X => null()
    monoMAT%B => null()
  end subroutine monolis_mat_nullify

  subroutine monolis_mat_initialize(monoMAT)
    implicit none
    type(monolis_mat) :: monoMAT

    monoMAT%N = 0
    monoMAT%NP = 0
    monoMAT%NZ = 0
    monoMAT%NDOF = 0
    monoMAT%index => null()
    monoMAT%item => null()
    monoMAT%A => null()
    monoMAT%X => null()
    monoMAT%B => null()
  end subroutine monolis_mat_initialize

  subroutine monolis_mat_copy(monoA, monoB)
    implicit none
    type(monolis_mat) :: monoA
    type(monolis_mat) :: monoB

    monoB%N = monoA%N
    monoB%NP = monoA%NP
    monoB%NZ = monoA%NZ
    monoB%NDOF = monoA%NDOF
    monoB%index => monoA%index
    monoB%item => monoA%item
    monoB%A => monoA%A
    monoB%X => monoA%X
    monoB%B => monoA%B
  end subroutine monolis_mat_copy

  subroutine monolis_mat_finalize(monoMAT)
    implicit none
    type(monolis_mat) :: monoMAT

    if(associated(monoMAT%index)) deallocate(monoMAT%index)
    if(associated(monoMAT%item)) deallocate(monoMAT%item)
    if(associated(monoMAT%A)) deallocate(monoMAT%A)
    if(associated(monoMAT%X)) deallocate(monoMAT%X)
    if(associated(monoMAT%B)) deallocate(monoMAT%B)
    monoMAT%index => null()
    monoMAT%item => null()
    monoMAT%A => null()
    monoMAT%X => null()
    monoMAT%B => null()
  end subroutine monolis_mat_finalize
end module mod_monolis_mat