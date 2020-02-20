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
    real(kind=kdouble), pointer :: X(:)
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

#ifdef TEST_ALL
    call monolis_solve_test(monoPRM, monoCOM, monoMAT)
#else
    call monolis_solve(monoPRM, monoCOM, monoMAT)
#endif
  end subroutine monolis
end module mod_monolis
