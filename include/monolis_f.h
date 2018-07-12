!interface
!  subroutine monolis(N, NP, NDOF, NPU, NPL, D, AU, AL, X, B, &
!    & indexU, itemU, indexL, itemL, myrank, comm, commsize, n_neib, &
!    & neib_pe, recv_index, send_index, recv_item, send_item, &
!    & method, precond, maxiter, tol, is_scaling)
!    implicit none
!    !> for monoMAT
!    integer, parameter :: kint = 4
!    integer, parameter :: kdouble = 8
!    integer(kind=kint), intent(in) :: N, NP, NDOF, NPU, NPL
!    integer(kind=kint), intent(in), pointer :: indexU(:)
!    integer(kind=kint), intent(in), pointer :: indexL(:)
!    integer(kind=kint), intent(in), pointer :: itemU(:)
!    integer(kind=kint), intent(in), pointer :: itemL(:)
!    real(kind=kdouble), intent(in), pointer :: D(:)
!    real(kind=kdouble), intent(in), pointer :: AU(:)
!    real(kind=kdouble), intent(in), pointer :: AL(:)
!    real(kind=kdouble), intent(in), pointer :: B(:)
!    real(kind=kdouble), intent(out),pointer :: X(:)
!    !> for monoCOM
!    integer(kind=kint), intent(in)          :: myrank
!    integer(kind=kint), intent(in)          :: comm
!    integer(kind=kint), intent(in)          :: commsize
!    integer(kind=kint), intent(in)          :: n_neib
!    integer(kind=kint), intent(in), pointer :: neib_pe(:)
!    integer(kind=kint), intent(in), pointer :: recv_index(:)
!    integer(kind=kint), intent(in), pointer :: recv_item(:)
!    integer(kind=kint), intent(in), pointer :: send_index(:)
!    integer(kind=kint), intent(in), pointer :: send_item(:)
!    !> for monoPRM
!    integer(kind=kint), intent(in) :: method
!    integer(kind=kint), intent(in) :: precond
!    integer(kind=kint), intent(in) :: maxiter
!    real(kind=kdouble), intent(in) :: tol
!    logical, intent(in) :: is_scaling
!  end subroutine monolis
!end interface
!
!interface
!  subroutine monolis_serial(N, NDOF, NPU, NPL, D, AU, AL, X, B, &
!    & indexU, itemU, indexL, itemL, &
!    & method, precond, maxiter, tol, is_scaling)
!    use mod_monolis_prm
!    use mod_monolis_com
!    use mod_monolis_mat
!    use mod_monolis_solve
!    implicit none
!    type(monolis_prm), save :: monoPRM
!    type(monolis_com), save :: monoCOM
!    type(monolis_mat), save :: monoMAT
!    !> for monoMAT
!    integer(kind=kint), intent(in) :: N, NDOF, NPU, NPL
!    integer(kind=kint), intent(in), pointer :: indexU(:)
!    integer(kind=kint), intent(in), pointer :: indexL(:)
!    integer(kind=kint), intent(in), pointer :: itemU(:)
!    integer(kind=kint), intent(in), pointer :: itemL(:)
!    real(kind=kdouble), intent(in), pointer :: D(:)
!    real(kind=kdouble), intent(in), pointer :: AU(:)
!    real(kind=kdouble), intent(in), pointer :: AL(:)
!    real(kind=kdouble), intent(in), pointer :: B(:)
!    real(kind=kdouble), intent(out),pointer :: X(:)
!    !> for monoPRM
!    integer(kind=kint), intent(in) :: method
!    integer(kind=kint), intent(in) :: precond
!    integer(kind=kint), intent(in) :: maxiter
!    real(kind=kdouble), intent(in) :: tol
!    logical, intent(in) :: is_scaling
!  end subroutine monolis_serial
!end interface
