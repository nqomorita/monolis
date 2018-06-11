interface
  subroutine monolis(N, NP, NDOF, NPU, NPL, D, AU, AL, X, B, &
    & indexU, itemU, indexL, itemL, myrank, comm, commsize, n_neib, &
    & neib_pe, recv_index, send_index, recv_item, send_item, &
    & method, precond, maxiter, tol, is_scaling)
    implicit none
    !> for monoMAT
    integer, parameter :: kint = 4
    integer, parameter :: kdouble = 8
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
  end subroutine monolis
end interface

interface
  subroutine monolis_wrapper(N, D, AU, AL, X, B, indexU, itemU, indexL, itemL, &
  & icomm, neib_pe, recv_index, send_index, recv_item, send_item, &
  & iparam, rparam, lparam)
    implicit none
    !> for monoMAT
    integer, parameter :: kint = 4
    integer, parameter :: kdouble = 8
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
  end subroutine monolis_wrapper
end interface
