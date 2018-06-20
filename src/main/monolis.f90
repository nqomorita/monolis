module mod_monolis
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_solve
end module mod_monolis

subroutine monolis(N, NP, NDOF, NPU, NPL, D, AU, AL, X, B, &
  & indexU, itemU, indexL, itemL, myrank, comm, commsize, n_neib, &
  & neib_pe, recv_index, send_index, recv_item, send_item, &
  & method, precond, maxiter, tol, is_scaling)
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_solve
  implicit none
  type(monolis_prm), save :: monoPRM
  type(monolis_com), save :: monoCOM
  type(monolis_mat), save :: monoMAT
  !> for monoMAT
  integer(kind=kint), intent(in) :: N, NP, NDOF, NPU, NPL
  integer(kind=kint), intent(in), pointer :: indexU(:)
  integer(kind=kint), intent(in), pointer :: indexL(:)
  integer(kind=kint), intent(in), pointer :: itemU(:)
  integer(kind=kint), intent(in), pointer :: itemL(:)
  real(kind=kdouble), intent(in), pointer :: D(:)
  real(kind=kdouble), intent(in), pointer :: AU(:)
  real(kind=kdouble), intent(in), pointer :: AL(:)
  real(kind=kdouble), intent(in), pointer :: B(:)
  real(kind=kdouble), intent(out),pointer :: X(:)
  !> for monoCOM
  integer(kind=kint), intent(in)          :: myrank
  integer(kind=kint), intent(in)          :: comm
  integer(kind=kint), intent(in)          :: commsize
  integer(kind=kint), intent(in)          :: n_neib
  integer(kind=kint), intent(in), pointer :: neib_pe(:)
  integer(kind=kint), intent(in), pointer :: recv_index(:)
  integer(kind=kint), intent(in), pointer :: recv_item(:)
  integer(kind=kint), intent(in), pointer :: send_index(:)
  integer(kind=kint), intent(in), pointer :: send_item(:)
  !> for monoPRM
  integer(kind=kint), intent(in) :: method
  integer(kind=kint), intent(in) :: precond
  integer(kind=kint), intent(in) :: maxiter
  real(kind=kdouble), intent(in) :: tol
  logical, intent(in) :: is_scaling

  !> for monoMAT
  monoMAT%N = N
  monoMAT%NP = NP
  monoMAT%NPU = NPU
  monoMAT%NPL = NPL
  monoMAT%NDOF = NDOF
  monoMAT%D  => D
  monoMAT%AU => AU
  monoMAT%AL => AL
  monoMAT%X  => X
  monoMAT%B  => B
  monoMAT%indexU => indexU
  monoMAT%indexL => indexL
  monoMAT%itemU => itemU
  monoMAT%itemL => itemL
  !> for monoCOM
  monoCOM%myrank = myrank
  monoCOM%comm = comm
  monoCOM%commsize = commsize
  monoCOM%n_neib = n_neib
  if(n_neib /= 0)then
    monoCOM%neib_pe => neib_pe
    monoCOM%recv_index => recv_index
    monoCOM%recv_item  => recv_item
    monoCOM%send_index => send_index
    monoCOM%send_item  => send_item
  endif
  !> for monoPRM
  monoPRM%method = method
  monoPRM%precond = precond
  monoPRM%maxiter = maxiter
  monoPRM%tol = tol
  monoPRM%is_scaling = is_scaling

  call monolis_solve(monoPRM, monoCOM, monoMAT)
  !call monolis_solve_test(monoPRM, monoCOM, monoMAT)
end subroutine monolis

subroutine monolis_serial(N, NDOF, NPU, NPL, D, AU, AL, X, B, &
  & indexU, itemU, indexL, itemL, &
  & method, precond, maxiter, tol, is_scaling)
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_solve
  implicit none
  type(monolis_prm), save :: monoPRM
  type(monolis_com), save :: monoCOM
  type(monolis_mat), save :: monoMAT
  !> for monoMAT
  integer(kind=kint), intent(in) :: N, NDOF, NPU, NPL
  integer(kind=kint), intent(in), pointer :: indexU(:)
  integer(kind=kint), intent(in), pointer :: indexL(:)
  integer(kind=kint), intent(in), pointer :: itemU(:)
  integer(kind=kint), intent(in), pointer :: itemL(:)
  real(kind=kdouble), intent(in), pointer :: D(:)
  real(kind=kdouble), intent(in), pointer :: AU(:)
  real(kind=kdouble), intent(in), pointer :: AL(:)
  real(kind=kdouble), intent(in), pointer :: B(:)
  real(kind=kdouble), intent(out),pointer :: X(:)
  !> for monoPRM
  integer(kind=kint), intent(in) :: method
  integer(kind=kint), intent(in) :: precond
  integer(kind=kint), intent(in) :: maxiter
  real(kind=kdouble), intent(in) :: tol
  logical, intent(in) :: is_scaling

  !> for monoMAT
  monoMAT%N = N
  monoMAT%NP = N
  monoMAT%NPU = NPU
  monoMAT%NPL = NPL
  monoMAT%NDOF = NDOF
  monoMAT%D  => D
  monoMAT%AU => AU
  monoMAT%AL => AL
  monoMAT%X  => X
  monoMAT%B  => B
  monoMAT%indexU => indexU
  monoMAT%indexL => indexL
  monoMAT%itemU => itemU
  monoMAT%itemL => itemL
  !> for monoCOM
  monoCOM%myrank = 0
  monoCOM%comm = 0
  monoCOM%commsize = 0
  monoCOM%n_neib = 0
  !> for monoPRM
  monoPRM%method = method
  monoPRM%precond = precond
  monoPRM%maxiter = maxiter
  monoPRM%tol = tol
  monoPRM%is_scaling = is_scaling

  call monolis_solve(monoPRM, monoCOM, monoMAT)
  !call monolis_solve_test(monoPRM, monoCOM, monoMAT)
end subroutine monolis_serial

subroutine monolis_wrapper(N, D, AU, AL, X, B, indexU, itemU, indexL, itemL, &
  & icomm, neib_pe, recv_index, send_index, recv_item, send_item, &
  & iparam, rparam, lparam)
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_solve
  implicit none
  type(monolis_prm), save :: monoPRM
  type(monolis_com), save :: monoCOM
  type(monolis_mat), save :: monoMAT
  !> for monoMAT
  integer(kind=kint), intent(in) :: N(5)
  integer(kind=kint), intent(in), pointer :: indexU(:)
  integer(kind=kint), intent(in), pointer :: indexL(:)
  integer(kind=kint), intent(in), pointer :: itemU(:)
  integer(kind=kint), intent(in), pointer :: itemL(:)
  real(kind=kdouble), intent(in), pointer :: D(:)
  real(kind=kdouble), intent(in), pointer :: AU(:)
  real(kind=kdouble), intent(in), pointer :: AL(:)
  real(kind=kdouble), intent(in), pointer :: B(:)
  real(kind=kdouble), intent(out),pointer :: X(:)
  !> for monoCOM
  integer(kind=kint), intent(in)          :: icomm(4)
  integer(kind=kint), intent(in), pointer :: neib_pe(:)
  integer(kind=kint), intent(in), pointer :: recv_index(:)
  integer(kind=kint), intent(in), pointer :: recv_item(:)
  integer(kind=kint), intent(in), pointer :: send_index(:)
  integer(kind=kint), intent(in), pointer :: send_item(:)
  !> for monoPRM
  integer(kind=kint), intent(in) :: iparam(3)
  real(kind=kdouble), intent(in) :: rparam(1)
  logical :: lparam(1)

  !> for monoMAT
  monoMAT%N = N(1)
  monoMAT%NP = N(2)
  monoMAT%NPU = N(3)
  monoMAT%NPL = N(4)
  monoMAT%NDOF = N(5)
  monoMAT%D  => D
  monoMAT%AU => AU
  monoMAT%AL => AL
  monoMAT%X  => X
  monoMAT%B  => B
  monoMAT%indexU => indexU
  monoMAT%indexL => indexL
  monoMAT%itemU => itemU
  monoMAT%itemL => itemL
  !> for monoCOM
  monoCOM%myrank = icomm(1)
  monoCOM%comm = icomm(2)
  monoCOM%commsize = icomm(3)
  monoCOM%n_neib = icomm(4)
  if(icomm(4) /= 0)then
    monoCOM%neib_pe => neib_pe
    monoCOM%recv_index => recv_index
    monoCOM%recv_item  => recv_item
    monoCOM%send_index => send_index
    monoCOM%send_item  => send_item
  endif
  !> for monoPRM
  monoPRM%method = iparam(1)
  monoPRM%precond = iparam(2)
  monoPRM%maxiter = iparam(3)
  monoPRM%tol = rparam(1)
  monoPRM%is_scaling = lparam(1)

  call monolis_solve(monoPRM, monoCOM, monoMAT)
  !call monolis_solve_test(monoPRM, monoCOM, monoMAT)
end subroutine monolis_wrapper
