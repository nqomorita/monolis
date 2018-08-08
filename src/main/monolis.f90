module mod_monolis
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_util
  use mod_monolis_solve
  !use mod_monolis_convert
  use mod_monolis_hash

contains

  subroutine monolis(monoCOM, N, NP, NZ, NDOF, A, X, B, index, item, &
    & method, precond, maxiter, tol, &
    & is_scaling, is_reordering, is_init_x, show_iterlog, show_time, show_summary)
    use mod_monolis_prm
    use mod_monolis_com
    use mod_monolis_mat
    use mod_monolis_solve
    implicit none
    type(monolis_com) :: monoCOM
    type(monolis_prm), save :: monoPRM
    type(monolis_mat), save :: monoMAT
    !> for monoMAT
    integer(kind=kint), intent(in) :: N, NP, NZ, NDOF
    integer(kind=kint), intent(in), pointer :: index(:)
    integer(kind=kint), intent(in), pointer :: item(:)
    real(kind=kdouble), intent(in), pointer :: A(:)
    real(kind=kdouble), intent(in), pointer :: B(:)
    real(kind=kdouble), intent(out),pointer :: X(:)
    !> for monoPRM
    integer(kind=kint), intent(in) :: method
    integer(kind=kint), intent(in) :: precond
    integer(kind=kint), intent(in) :: maxiter
    real(kind=kdouble), intent(in) :: tol
    logical, intent(in) :: is_scaling
    logical, intent(in) :: is_reordering
    logical, intent(in) :: is_init_x
    logical, intent(in) :: show_iterlog
    logical, intent(in) :: show_time
    logical, intent(in) :: show_summary

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
    !> for monoPRM
    monoPRM%method = method
    monoPRM%precond = precond
    monoPRM%maxiter = maxiter
    monoPRM%tol = tol
    monoPRM%is_scaling    = is_scaling
    monoPRM%is_reordering = is_reordering
    monoPRM%is_init_x     = is_init_x
    monoPRM%show_iterlog  = show_iterlog
    monoPRM%show_time     = show_time
    monoPRM%show_summary  = show_summary

#ifdef DTEST_ALL
    call monolis_solve_test(monoPRM, monoCOM, monoMAT)
#else
    call monolis_solve(monoPRM, monoCOM, monoMAT)
#endif
  end subroutine monolis

  subroutine monolis_c(monoCOM_c, N, NP, NZ, NDOF, A, X, B, index, item, &
    & method, precond, maxiter, tol, &
    & is_scaling, is_reordering, is_init_x, show_iterlog, show_time, show_summary) &
    & bind(c, name="monolis")
    use iso_c_binding
    use mod_monolis_prm
    use mod_monolis_com
    use mod_monolis_mat
    use mod_monolis_solve
    implicit none
    type(monolis_com_c) :: monoCOM_c
    type(monolis_com), save :: monoCOM
    type(monolis_prm), save :: monoPRM
    type(monolis_mat), save :: monoMAT
    !> for monoMAT
    integer(c_int), value :: N, NP, NZ, NDOF
    integer(c_int), intent(in), target :: index(0:NP)
    integer(c_int), intent(in), target :: item(NZ)
    real(c_double), intent(in), target :: A(NZ*NDOF*NDOF)
    real(c_double), intent(in), target :: B(NP*NDOF)
    real(c_double), intent(out),target :: X(NP*NDOF)
    !> for monoPRM
    integer(c_int), value :: method, precond, maxiter
    real(c_double), value :: tol
    integer(c_int), value :: is_scaling, is_reordering, is_init_x
    integer(c_int), value :: show_iterlog, show_time, show_summary

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
    monoMAT%item = monoMAT%item + 1
    !> for monoCOM
    monoCOM%myrank = monoCOM_c%myrank
    monoCOM%comm = monoCOM_c%comm
    monoCOM%commsize = monoCOM_c%commsize
    monoCOM%n_neib = monoCOM_c%n_neib
    !> for monoPRM
    monoPRM%method = method
    monoPRM%precond = precond
    monoPRM%maxiter = maxiter
    monoPRM%tol = tol
    monoPRM%is_scaling    = .false.
    monoPRM%is_reordering = .false.
    monoPRM%is_init_x     = .false.
    monoPRM%show_iterlog  = .false.
    monoPRM%show_time     = .false.
    monoPRM%show_summary  = .false.
    if(is_scaling     == 1) monoPRM%is_scaling    = .true.
    if(is_reordering  == 1) monoPRM%is_reordering = .true.
    if(is_init_x      == 1) monoPRM%is_init_x     = .true.
    if(show_iterlog   == 1) monoPRM%show_iterlog  = .true.
    if(show_time      == 1) monoPRM%show_time     = .true.
    if(show_summary   == 1) monoPRM%show_summary  = .true.

#ifdef TEST_ALL
    call monolis_solve_test(monoPRM, monoCOM, monoMAT)
#else
    call monolis_solve(monoPRM, monoCOM, monoMAT)
#endif
    monoMAT%item = monoMAT%item - 1
  end subroutine monolis_c
end module mod_monolis
