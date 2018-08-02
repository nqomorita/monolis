module mod_monolis
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_util
  use mod_monolis_solve
  use mod_monolis_convert
  use mod_monolis_hash

contains

  subroutine monolis(N, NP, NZ, NDOF, A, X, B, index, item, &
    & myrank, comm, commsize, n_neib, neib_pe, recv_index, send_index, recv_item, send_item, &
    & method, precond, maxiter, tol, &
    & is_scaling, is_reordering, is_init_x, show_iterlog, show_time, show_summary)
    use mod_monolis_prm
    use mod_monolis_com
    use mod_monolis_mat
    use mod_monolis_solve
    implicit none
    type(monolis_prm), save :: monoPRM
    type(monolis_com), save :: monoCOM
    type(monolis_mat), save :: monoMAT
    !> for monoMAT
    integer(kind=kint), intent(in) :: N, NP, NZ, NDOF
    integer(kind=kint), intent(in), pointer :: index(:)
    integer(kind=kint), intent(in), pointer :: item(:)
    real(kind=kdouble), intent(in), pointer :: A(:)
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
    integer(kind=kint), intent(in) :: is_reordering
    integer(kind=kint), intent(in) :: is_init_x
    integer(kind=kint), intent(in) :: show_iterlog
    integer(kind=kint), intent(in) :: show_time
    integer(kind=kint), intent(in) :: show_summary
    real(kind=kdouble), intent(in) :: tol

    !> for monoMAT
    monoMAT%N = N
    monoMAT%NP = NP
    monoMAT%NZ = NZ
    monoMAT%NDOF = NDOF
    monoMAT%A => A
    monoMAT%X => X
    monoMAT%B => B
    monoMAT%index => index
    monoMAT%item => item
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
    monoPRM%is_scaling     = .false.
    monoPRM%is_reordering  = .false.
    monoPRM%is_init_x      = .false.
    monoPRM%show_iterlog = .false.
    monoPRM%show_time      = .false.
    monoPRM%show_summary   = .false.
    if(is_scaling     == 1) monoPRM%is_scaling     = .true.
    if(is_reordering  == 1) monoPRM%is_reordering  = .true.
    if(is_init_x      == 1) monoPRM%is_init_x      = .true.
    if(show_iterlog   == 1) monoPRM%show_iterlog   = .true.
    if(show_time      == 1) monoPRM%show_time      = .true.
    if(show_summary   == 1) monoPRM%show_summary   = .true.

#ifdef DTEST_ALL
    call monolis_solve_test(monoPRM, monoCOM, monoMAT)
#else
    call monolis_solve(monoPRM, monoCOM, monoMAT)
#endif
  end subroutine monolis

  subroutine monolis_c(N, NP, NZ, NDOF, A, X, B, index, item, &
    & myrank, comm, commsize, n_neib, neib_pe, recv_index, send_index, recv_item, send_item, &
    & method, precond, maxiter, tol, &
    & is_scaling, is_reordering, is_init_x, show_iterlog, show_time, show_summary) &
    & bind(c, name="monolis_serial")
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
    integer(c_int), value   :: N, NP, NZ, NDOF
    integer(c_int), intent(in), target :: index(0:N)
    integer(c_int), intent(in), target :: item(NZ)
    real(c_double), intent(in), target :: A(NZ*NDOF*NDOF)
    real(c_double), intent(in), target :: B(N*NDOF)
    real(c_double), intent(out),target :: X(N*NDOF)
    !> for monoCOM
    integer(c_int), intent(in), value  :: myrank, comm, commsize, n_neib
    integer(c_int), intent(in), target :: neib_pe(n_neib)
    integer(c_int), intent(in), target :: recv_index(0:n_neib)
    integer(c_int), intent(in), target :: recv_item(:)
    integer(c_int), intent(in), target :: send_index(0:n_neib)
    integer(c_int), intent(in), target :: send_item(:)
    !> for monoPRM
    integer(c_int), value :: method, precond, maxiter
    integer(c_int), value :: is_scaling, is_reordering, is_init_x
    integer(c_int), value :: show_iterlog, show_time, show_summary
    real(c_double), value :: tol

    !> for monoMAT
    monoMAT%N = N
    monoMAT%NP = NP
    monoMAT%NZ = NZ
    monoMAT%NDOF = NDOF
    monoMAT%A => A
    monoMAT%X => X
    monoMAT%B => B
    monoMAT%index => index
    monoMAT%item => item
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
    monoPRM%is_scaling = .false.
    monoPRM%is_reordering = .false.
    monoPRM%is_init_x = .false.
    monoPRM%show_iterlog = .false.
    monoPRM%show_time = .false.
    monoPRM%show_summary = .false.
    if(is_scaling     == 1) monoPRM%is_scaling    = .true.
    if(is_reordering  == 1) monoPRM%is_reordering = .true.
    if(is_init_x      == 1) monoPRM%is_init_x     = .true.
    if(show_iterlog   == 1) monoPRM%show_iterlog  = .true.
    if(show_time      == 1) monoPRM%show_time     = .true.
    if(show_summary   == 1) monoPRM%show_summary  = .true.

#ifdef DTEST_ALL
    call monolis_solve_test(monoPRM, monoCOM, monoMAT)
#else
    call monolis_solve(monoPRM, monoCOM, monoMAT)
#endif
  end subroutine monolis_c
end module mod_monolis
