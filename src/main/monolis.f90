module mod_monolis
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_util
  use mod_monolis_solve
  use mod_monolis_convert
  use mod_monolis_hash

contains

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
    integer(kind=kint), intent(in) :: is_scaling
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
    if(is_scaling == 1) monoPRM%is_scaling = .true.
    monoPRM%is_reordering = .true.

#ifdef DTEST_ALL
    call monolis_solve_test(monoPRM, monoCOM, monoMAT)
#else
    call monolis_solve(monoPRM, monoCOM, monoMAT)
#endif
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
    integer(kind=kint), intent(in) :: is_scaling
    real(kind=kdouble), intent(in) :: tol

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
    monoCOM%neib_pe => NULL()
    monoCOM%recv_index => NULL()
    monoCOM%recv_item  => NULL()
    monoCOM%send_index => NULL()
    monoCOM%send_item  => NULL()
    !> for monoPRM
    monoPRM%method = method
    monoPRM%precond = precond
    monoPRM%maxiter = maxiter
    monoPRM%tol = tol
    if(is_scaling == 1) monoPRM%is_scaling = .true.
    monoPRM%is_reordering = .true.

    call monolis_com_initialize(monoCOM)

#ifdef DTEST_ALL
    call monolis_solve_test(monoPRM, monoCOM, monoMAT)
#else
    call monolis_solve(monoPRM, monoCOM, monoMAT)
#endif
  end subroutine monolis_serial

  subroutine monolis_serial_c(N, NDOF, NPU, NPL, D, AU, AL, X, B, &
    & indexU, itemU, indexL, itemL, &
    & method, precond, maxiter, tol, is_scaling) bind(c, name="monolis_serial")
    use iso_c_binding
    use mod_monolis_prm
    use mod_monolis_com
    use mod_monolis_mat
    use mod_monolis_solve
    implicit none
    type(monolis_prm), save :: monoPRM
    type(monolis_com), save :: monoCOM
    type(monolis_mat), save :: monoMAT
    !> for monoMAT
    integer(c_int), value   :: N, NDOF, NPU, NPL
    integer(c_int), intent(in), target :: indexU(0:N)
    integer(c_int), intent(in), target :: indexL(0:N)
    integer(c_int), intent(in), target :: itemU(NPU)
    integer(c_int), intent(in), target :: itemL(NPL)
    real(c_double), intent(in), target :: D(N*NDOF*NDOF)
    real(c_double), intent(in), target :: AU(NPU*NDOF*NDOF)
    real(c_double), intent(in), target :: AL(NPL*NDOF*NDOF)
    real(c_double), intent(in), target :: B(N*NDOF)
    real(c_double), intent(out),target :: X(N*NDOF)
    !> for monoPRM
    integer(c_int), value :: method, precond, maxiter, is_scaling
    real(c_double), value :: tol

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
    monoCOM%neib_pe => NULL()
    monoCOM%recv_index => NULL()
    monoCOM%recv_item  => NULL()
    monoCOM%send_index => NULL()
    monoCOM%send_item  => NULL()
    !> for monoPRM
    monoPRM%method = method
    monoPRM%precond = precond
    monoPRM%maxiter = maxiter
    monoPRM%tol = tol
    if(is_scaling == 1) monoPRM%is_scaling = .true.
    monoPRM%is_reordering = .true.

    call monolis_com_initialize(monoCOM)

#ifdef DTEST_ALL
    call monolis_solve_test(monoPRM, monoCOM, monoMAT)
#else
    call monolis_solve(monoPRM, monoCOM, monoMAT)
#endif
  end subroutine monolis_serial_c
end module mod_monolis
