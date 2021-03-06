program main
  use mod_monolis
  implicit none
  type(monolis_com) :: monoCOM
  integer(kind=kint) :: i, ni, N, NZ, NDOF
  integer(kind=kint) :: method, precond, maxiter
  logical :: is_scaling, is_reordering, is_init_x
  logical :: show_iterlog, show_time, show_summary
  integer(kind=kint), pointer :: indexI(:) => NULL()
  integer(kind=kint), pointer :: indexJ(:) => NULL()
  integer(kind=kint), pointer :: index(:) => NULL()
  integer(kind=kint), pointer :: item(:) => NULL()
  real(kind=kdouble), pointer :: A(:) => NULL()
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

  call monolis_com_initialize(monoCOM)

  allocate(index(0:N))
  allocate(item(NZ))
  index = 0
  item  = 0

  do i = 1, NZ
    ni = indexI(i)
    index(ni) = index(ni) + 1
    item(i) = indexJ(i)
  enddo

  do i = 1, N
    index(i) = index(i) + index(i-1)
  enddo

  allocate(X(N))
  allocate(B(N))

  X = 0.0d0
  B = 1.0d0
  method = 1
  precond = 1
  maxiter = 1000
  tol = 1.0d-8
  is_scaling    = .false.
  is_reordering = .false.
  is_init_x     = .false.
  show_iterlog  = .true.
  show_time     = .true.
  show_summary  = .true.

  write(*,"(a)")"* call monolis"
  call monolis(monoCOM, N, N, NZ, NDOF, A, X, B, index, item, &
    & method, precond, maxiter, tol, &
    & is_scaling, is_reordering, is_init_x, show_iterlog, show_time, show_summary)

  write(*,"(a)")"* monolis result"
  write(*,"(1p3e12.5)")X

  B = 0.0d0
  call monolis_matvec_wrapper(monoCOM, N, N, NZ, NDOF, A, index, item, X, B)

  write(*,"(a)")"* monolis b = Ax"
  write(*,"(1p3e12.5)")B

  deallocate(X)
  deallocate(B)

  call monolis_com_finalize(monoCOM)
end program main
