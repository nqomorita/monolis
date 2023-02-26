!> 行列構造体の定義
module mod_monolis_def_mat
  use mod_monolis_utils
  implicit none

  !> 行列構造体（セパレート CSR 構造）
  type monolis_mat_LDU
    integer(kint) :: N, NP, NPU, NPL, NDOF
    integer(kint), pointer :: indexU(:) => null()
    integer(kint), pointer :: itemU(:) => null()
    integer(kint), pointer :: indexL(:) => null()
    integer(kint), pointer :: itemL(:) => null()
    real(kdouble), pointer :: U(:) => null()
    real(kdouble), pointer :: D(:) => null()
    real(kdouble), pointer :: L(:) => null()
    real(kdouble), pointer :: X(:) => null()
    real(kdouble), pointer :: B(:) => null()
  end type monolis_mat_LDU

  !> 行列構造体（CSR 構造）
  type monolis_mat
    type(monolis_mat_LDU) :: monoTREE
    integer(kint) :: N, NP, NZ, NDOF
    integer(kint), pointer :: index(:) => null()
    integer(kint), pointer :: item(:) => null()
    integer(kint), pointer :: perm(:) => null()
    integer(kint), pointer :: iperm(:) => null()

    !> for CRR format
    integer(kint), pointer :: indexR(:) => null()
    integer(kint), pointer :: itemR(:) => null()
    integer(kint), pointer :: permR(:) => null()

    real(kdouble), pointer :: A(:) => null()
    real(kdouble), pointer :: X(:) => null()
    real(kdouble), pointer :: B(:) => null()
    real(kdouble), pointer :: diag(:) => null()
  end type monolis_mat

contains

  subroutine monolis_mat_initialize(monoMAT)
    implicit none
    type(monolis_mat) :: monoMAT

    monoMAT%N = 0
    monoMAT%NP = 0
    monoMAT%NZ = 0
    monoMAT%NDOF = 0
    monoMAT%index => null()
    monoMAT%item => null()
    monoMAT%perm => null()
    monoMAT%iperm => null()
    monoMAT%indexR => null()
    monoMAT%itemR => null()
    monoMAT%permR => null()
    monoMAT%A => null()
    monoMAT%X => null()
    monoMAT%B => null()
  end subroutine monolis_mat_initialize

  subroutine monolis_mat_finalize(monoMAT)
    implicit none
    type(monolis_mat) :: monoMAT

    call monolis_mat_tree_finalize(monoMAT%monoTREE)

    if(associated(monoMAT%index)) deallocate(monoMAT%index)
    if(associated(monoMAT%item)) deallocate(monoMAT%item)
    if(associated(monoMAT%A)) deallocate(monoMAT%A)
    if(associated(monoMAT%X)) deallocate(monoMAT%X)
    if(associated(monoMAT%B)) deallocate(monoMAT%B)
    if(associated(monoMAT%diag)) deallocate(monoMAT%diag)
    if(associated(monoMAT%perm)) deallocate(monoMAT%perm)
    if(associated(monoMAT%iperm)) deallocate(monoMAT%iperm)
    if(associated(monoMAT%indexR)) deallocate(monoMAT%indexR)
    if(associated(monoMAT%itemR)) deallocate(monoMAT%itemR)
    if(associated(monoMAT%permR)) deallocate(monoMAT%permR)
    monoMAT%index => null()
    monoMAT%item => null()
    monoMAT%A => null()
    monoMAT%X => null()
    monoMAT%B => null()
    monoMAT%perm => null()
    monoMAT%iperm => null()
    monoMAT%diag => null()
    monoMAT%indexR => null()
    monoMAT%itemR => null()
    monoMAT%permR => null()
  end subroutine monolis_mat_finalize

  subroutine monolis_mat_tree_finalize(monoMAT)
    implicit none
    type(monolis_mat_LDU) :: monoMAT

    if(associated(monoMAT%indexU)) deallocate(monoMAT%indexU)
    if(associated(monoMAT%indexL)) deallocate(monoMAT%indexL)
    if(associated(monoMAT%itemU)) deallocate(monoMAT%itemU)
    if(associated(monoMAT%itemL)) deallocate(monoMAT%itemL)
    if(associated(monoMAT%L)) deallocate(monoMAT%L)
    if(associated(monoMAT%D)) deallocate(monoMAT%D)
    if(associated(monoMAT%U)) deallocate(monoMAT%U)
    if(associated(monoMAT%X)) deallocate(monoMAT%X)
    if(associated(monoMAT%B)) deallocate(monoMAT%B)
    monoMAT%indexU => null()
    monoMAT%indexL => null()
    monoMAT%itemU => null()
    monoMAT%itemL => null()
    monoMAT%L => null()
    monoMAT%D => null()
    monoMAT%U => null()
    monoMAT%X => null()
    monoMAT%B => null()
  end subroutine monolis_mat_tree_finalize

end module mod_monolis_def_mat