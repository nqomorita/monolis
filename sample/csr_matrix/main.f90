program main
  use mod_monolis
  implicit none
  type(monolis_prm) :: monoPRM
  type(monolis_com) :: monoCOM
  type(monolis_mat) :: monoMAT
  integer(kind=kint) :: i, j, N, NZ, NDOF, NPU, NPL
  integer(kind=kint) :: method, precond, maxiter, is_scaling
  integer(kind=kint), pointer :: index(:) => NULL()
  integer(kind=kint), pointer :: item(:) => NULL()
  integer(kind=kint), pointer :: indexU(:) => NULL()
  integer(kind=kint), pointer :: indexL(:) => NULL()
  integer(kind=kint), pointer :: itemU(:) => NULL()
  integer(kind=kint), pointer :: itemL(:) => NULL()
  real(kind=kdouble), pointer :: A(:) => NULL()
  real(kind=kdouble), pointer :: D(:) => NULL()
  real(kind=kdouble), pointer :: AU(:) => NULL()
  real(kind=kdouble), pointer :: AL(:) => NULL()
  real(kind=kdouble), pointer :: X(:) => NULL()
  real(kind=kdouble), pointer :: B(:) => NULL()
  real(kind=kdouble) :: tol
  character :: ctemp*10

  write(*,"(a)")"* monolis csr_matrix test"

  !> reference: matrix market
  !> https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/bcsstruc2/bcsstk14.html
  open(20, file="./test.mtx", status="old")
    read(20,*)ctemp
    read(20,*)N, NZ
    allocate(index(0:N))
    allocate(item(NZ))
    allocate(A(NZ))
    read(20,*)(index(i), i = 0, N)
    read(20,*)(item(i), i = 1, NZ)
    read(20,*)(A(i), i = 1, NZ)
  close(20)

  NDOF = 1

  call monolis_convert_csr_get_size(N, NZ, index, item, NPU, NPL)
  call monolis_convert_alloc_matrix(N, NDOF, NPU, NPL, D, AU, AL, indexU, indexL, itemU, itemL, X, B)
  call monolis_convert_csr_get_index(N, NZ, index, item, NPU, NPL, indexU, itemU, indexL, itemL)
  call monolis_convert_csr_update_matrix_entry(N, NZ, NDOF, A, index, item, NPU, NPL, &
       & D, AU, AL, indexU, itemU, indexL, itemL)

  B = 1.0d0
  method = 1
  precond = 1
  maxiter = 1000
  tol = 1.0d-8
  is_scaling = 1

  call monolis_serial(N, NDOF, NPU, NPL, D, AU, AL, X, B, &
  & indexU, itemU, indexL, itemL, &
  & method, precond, maxiter, tol, is_scaling)

  call monolis_convert_dealloc_matrix(D, AU, AL, indexU, indexL, itemU, itemL, X, B)
end program main
