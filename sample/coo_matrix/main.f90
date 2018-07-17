program main
  use mod_monolis
  implicit none
  type(monolis_prm) :: monoPRM
  type(monolis_com) :: monoCOM
  type(monolis_mat) :: monoMAT
  integer(kind=kint) :: i, j, N, NDOF, NPU, NPL, NZ
  integer(kind=kint) :: method, precond, maxiter
  integer(kind=kint) :: is_scaling, is_reordering, is_init_x, show_iteration
  integer(kind=kint), pointer :: indexI(:) => NULL()
  integer(kind=kint), pointer :: indexJ(:) => NULL()
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

  write(*,"(a)")"* monolis coo_matrix test"

  !> reference: matrix market
  !> https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/bcsstruc2/bcsstk14.html
  open(20, file="./test.mtx", status="old")
    read(20,*)ctemp
    read(20,*)N, N, NZ
    allocate(indexI(NZ))
    allocate(indexJ(NZ))
    allocate(A(NZ))
    do i = 1, NZ
      read(20,*)indexI(i), indexJ(i), A(i)
    enddo
  close(20)

  NDOF = 1

  call monolis_convert_coo_get_size(N, NZ, indexI, indexJ, NPU, NPL)
  call monolis_convert_alloc_matrix(N, NDOF, NPU, NPL, D, AU, AL, indexU, indexL, itemU, itemL, X, B)
  call monolis_convert_coo_get_index(N, NZ, indexI, indexJ, NPU, NPL, indexU, itemU, indexL, itemL)
  call monolis_convert_coo_update_matrix_entry(N, NZ, NDOF, A, indexI, indexJ, NPU, NPL, &
       & D, AU, AL, indexU, itemU, indexL, itemL)

  B = 1.0d0
  method = 1
  precond = 1
  maxiter = 1000
  tol = 1.0d-8
  is_scaling = 1
  is_reordering = 1
  is_init_x = 1
  show_iteration = 1

  call monolis_serial(N, NDOF, NPU, NPL, D, AU, AL, X, B, &
  & indexU, itemU, indexL, itemL, &
  & method, precond, maxiter, tol, is_scaling, is_reordering, is_init_x, show_iteration)

  write(*,"(a)")"* monolis result"
  write(*,"(1p3e12.5)")X

  B = 0.0d0
  call monolis_matvec_serial(N, NDOF, NPU, NPL, D, AU, AL, indexU, itemU, indexL, itemL, X, B)

  write(*,"(a)")"* monolis b = Ax"
  write(*,"(1p3e12.5)")B

  call  monolis_convert_dealloc_matrix(D, AU, AL, indexU, indexL, itemU, itemL, X, B)
end program main
