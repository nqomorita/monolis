module mod_monolis_mat
  use mod_monolis_prm
  implicit none

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

  type monolis_mat
    type(monolis_mat_LDU) :: monoTREE
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

  subroutine monolis_mat_copy_by_pointer(min, mout)
    implicit none
    type(monolis_mat) :: min
    type(monolis_mat) :: mout

    mout%N = min%N
    mout%NP = min%NP
    mout%NZ = min%NZ
    mout%NDOF = min%NDOF
    mout%index => min%index
    mout%item => min%item
    mout%A => min%A
    mout%X => min%X
    mout%B => min%B
  end subroutine monolis_mat_copy_by_pointer

  subroutine monolis_mat_copy_all(min, mout)
    implicit none
    type(monolis_mat) :: min
    type(monolis_mat) :: mout
    integer(kind=kint) :: i, NZ

    mout%N = min%N
    mout%NP = min%NP
    mout%NZ = min%NZ
    mout%NDOF = min%NDOF

    NZ = min%index(min%NP)
    allocate(mout%index(0:min%NP))
    allocate(mout%item(NZ))
    allocate(mout%A(min%NDOF*min%NDOF*NZ))
    allocate(mout%X(min%NDOF*min%NP))
    allocate(mout%B(min%NDOF*min%NP))

    mout%index(0:min%NP) = min%index(0:min%NP)
    mout%item = min%item
    mout%A = min%A
    mout%X = min%X
    mout%B = min%B
  end subroutine monolis_mat_copy_all

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