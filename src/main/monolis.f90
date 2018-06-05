module mod_monolis
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_solve
end module mod_monolis

subroutine monolis(N, NP, NDOF, NPU, NPL, D, AU, AL, X, B, &
  & indexU, itemU, indexL, itemL, myrank, comm, commsize, n_neib, &
  & neib_pe, recv_index, send_index, recv_item, send_item, &
  & method, precond, maxiter, tol)
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
  real(kind=kdouble), intent(out), pointer :: X(:)
  !> for monoCOM
  integer(kind=kint), intent(in)         :: myrank
  integer(kind=kint), intent(in)         :: comm
  integer(kind=kint), intent(in)         :: commsize
  integer(kind=kint), intent(in)         :: n_neib
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

  call monolis_solve(monoPRM, monoCOM, monoMAT)

end subroutine monolis
