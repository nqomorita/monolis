module mod_monolis_alloc_matrix
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat

  implicit none

contains

  subroutine monolis_convert_alloc_matrix(NP, NDOF, NPU, NPL, D, AU, AL, indexU, indexL, itemU, itemL, X, B)
    implicit none
    real(kind=kdouble), pointer :: X(:), B(:)
    real(kind=kdouble), pointer :: D(:), AU(:), AL(:)
    integer(kind=kint), pointer :: indexU(:)
    integer(kind=kint), pointer :: indexL(:)
    integer(kind=kint), pointer :: itemU(:)
    integer(kind=kint), pointer :: itemL(:)
    integer(kind=kint) :: NP, NDOF, NPU, NPL

    allocate(X(NP*NDOF))
    allocate(B(NP*NDOF))
    allocate(D(NP*NDOF*NDOF))
    allocate(AU(NPU*NDOF*NDOF))
    allocate(AL(NPL*NDOF*NDOF))
    allocate(indexU(0:NP))
    allocate(indexL(0:NP))
    allocate(itemU(NPU))
    allocate(itemL(NPL))
    X = 0.0d0
    B = 0.0d0
    D = 0.0d0
    AU = 0.0d0
    AL = 0.0d0
    indexU = 0
    indexL = 0
    itemU = 0
    itemL = 0
  end subroutine monolis_convert_alloc_matrix
end module mod_monolis_alloc_matrix

module mod_monolis_convert
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_convert_full
  use mod_monolis_convert_coo
  use mod_monolis_convert_csr
  use mod_monolis_alloc_matrix
end module mod_monolis_convert
