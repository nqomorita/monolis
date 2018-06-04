module mod_monolis
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_solve
end module mod_monolis

subroutine monolis(N, NP, NPU, NPL, NDOF, indexU, indexL, itemU, itemL, D, AU, AL, X, B, &
  & myrank, comm, commsize, n_neib, neib_pe, recv_index, send_index, recv_item, send_item, &
  & method, precond, maxiter, tol)
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_solve
  implicit none
  type(monolis_prm) :: monoPRM
  type(monolis_com) :: monoCOM
  type(monolis_mat) :: monoMAT
  !> for monoMAT
  integer(kind=kint) :: N, NP, NPU, NPL, NDOF
  integer(kind=kint), pointer :: indexU(:)
  integer(kind=kint), pointer :: indexL(:)
  integer(kind=kint), pointer :: itemU(:)
  integer(kind=kint), pointer :: itemL(:)
  real(kind=kdouble), pointer :: D(:)
  real(kind=kdouble), pointer :: AU(:)
  real(kind=kdouble), pointer :: AL(:)
  real(kind=kdouble), pointer :: X(:)
  real(kind=kdouble), pointer :: B(:)
  !> for monoCOM
  integer(kind=kint)          :: myrank
  integer(kind=kint)          :: comm
  integer(kind=kint)          :: commsize
  integer(kind=kint)          :: n_neib
  integer(kind=kint), pointer :: neib_pe(:)
  integer(kind=kint), pointer :: recv_index(:)
  integer(kind=kint), pointer :: recv_item(:)
  integer(kind=kint), pointer :: send_index(:)
  integer(kind=kint), pointer :: send_item(:)
  !> for monoPRM
  integer(kind=kint) :: method
  integer(kind=kint) :: precond
  integer(kind=kint) :: maxiter
  real(kind=kdouble) :: tol

  !> for monoMAT
  monoMAT%N = N
  monoMAT%NP = NP
  monoMAT%NPU = NPU
  monoMAT%NPL = NPL
  monoMAT%NDOF = NDOF
  monoMAT%indexU => indexU
  monoMAT%indexL => indexL
  monoMAT%itemU  => itemU
  monoMAT%itemL  => itemL
  monoMAT%D  => D
  monoMAT%AU => AU
  monoMAT%AL => AL
  monoMAT%X  => X
  monoMAT%B  => B

  !> for monoCOM
  monoCOM%myrank = myrank
  monoCOM%comm = comm
  monoCOM%commsize = commsize
  monoCOM%n_neib = n_neib
  if(1 < n_neib)then
    monoCOM%neib_pe => neib_pe
    monoCOM%recv_index => recv_index
    monoCOM%send_index => send_index
    monoCOM%recv_item  => recv_item
    monoCOM%send_item  => send_item
  endif

  !> for monoPRM
  monoPRM%method = method
  monoPRM%precond = precond
  monoPRM%maxiter = maxiter
  monoPRM%tol = tol

  call monolis_solve(monoPRM, monoCOM, monoMAT)

end subroutine monolis
