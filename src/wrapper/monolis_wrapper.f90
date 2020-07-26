module mod_monolis_wrapper
  use iso_c_binding
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_util
  use mod_monolis_stdlib
  implicit none

contains

  !> initialize section
  subroutine monolis_global_initialize_c() &
    & bind(c, name = "monolis_global_initialize")
    implicit none
    call monolis_global_initialize()
  end subroutine monolis_global_initialize_c

  subroutine monolis_global_finalize_c() &
    & bind(c, name = "monolis_global_finalize")
    implicit none
    call monolis_global_finalize()
  end subroutine monolis_global_finalize_c

  !> mat
  subroutine monolis_add_sparse_matrix_c() &
    & bind(c, name = "monolis_add_sparse_matrix_c_main")
    implicit none
    type(monolis_structure) :: monolis


  end subroutine monolis_add_sparse_matrix_c

  subroutine monolis_set_Dirichlet_bc_c() &
    & bind(c, name = "monolis_set_Dirichlet_bc_c_main")
    implicit none
    type(monolis_structure) :: monolis


  end subroutine monolis_set_Dirichlet_bc_c

  !> solve
  subroutine monolis_solve_c() &
    & bind(c, name = "monolis_solve_c_main")
    implicit none
    type(monolis_structure) :: monolis


  end subroutine monolis_solve_c

  !> std lib
  subroutine monolis_qsort_int_c(array, iS, iE) &
    & bind(c, name = "monolis_qsort_int")
    implicit none
    integer(c_int), target :: array(1:iE-iS+1)
    integer(c_int), value :: iS, iE

    iS = iS
    iE = iE
    call monolis_qsort_int(array, iS, iE)
  end subroutine monolis_qsort_int_c

!  subroutine monolis_c(monoCOM_c, N, NP, NZ, NDOF, A, X, B, index, item, &
!    & method, precond, maxiter, tol, &
!    & is_scaling, is_reordering, is_init_x, show_iterlog, show_time, show_summary) &
!    & bind(c, name="monolis")
!    use iso_c_binding
!    use mod_monolis_prm
!    use mod_monolis_com
!    use mod_monolis_mat
!    use mod_monolis_solve
!    implicit none
!    type(monolis_com_c) :: monoCOM_c
!    type(monolis_com) :: monoCOM
!    type(monolis_prm) :: monoPRM
!    type(monolis_mat) :: monoMAT
!    !> for monoMAT
!    integer(c_int), value :: N, NP, NZ, NDOF
!    integer(c_int), intent(in), target :: index(0:NP)
!    integer(c_int), intent(in), target :: item(NZ)
!    real(c_double), intent(in), target :: A(NZ*NDOF*NDOF)
!    real(c_double), intent(in), target :: B(NP*NDOF)
!    real(c_double), intent(out),target :: X(NP*NDOF)
!    !> for monoPRM
!    integer(c_int), value :: method, precond, maxiter
!    real(c_double), value :: tol
!    integer(c_int), value :: is_scaling, is_reordering, is_init_x
!    integer(c_int), value :: show_iterlog, show_time, show_summary
!
!    !> for monoMAT
!    monoMAT%N = N
!    monoMAT%NP = NP
!    monoMAT%NZ = NZ
!    monoMAT%NDOF = NDOF
!    monoMAT%A => A
!    monoMAT%X => X
!    monoMAT%B => B
!    monoMAT%index => index
!    monoMAT%item => item
!    monoMAT%item = monoMAT%item + 1
!    !> for monoCOM
!    monoCOM%myrank = monoCOM_c%myrank
!    monoCOM%comm = monoCOM_c%comm
!    monoCOM%commsize = monoCOM_c%commsize
!    monoCOM%n_neib = monoCOM_c%n_neib
!    !> for monoPRM
!    monoPRM%method = method
!    monoPRM%precond = precond
!    monoPRM%maxiter = maxiter
!    monoPRM%tol = tol
!    monoPRM%is_scaling    = .false.
!    monoPRM%is_reordering = .false.
!    monoPRM%is_init_x     = .false.
!    monoPRM%show_iterlog  = .false.
!    monoPRM%show_time     = .false.
!    monoPRM%show_summary  = .false.
!    if(is_scaling     == 1) monoPRM%is_scaling    = .true.
!    if(is_reordering  == 1) monoPRM%is_reordering = .true.
!    if(is_init_x      == 1) monoPRM%is_init_x     = .true.
!    if(show_iterlog   == 1) monoPRM%show_iterlog  = .true.
!    if(show_time      == 1) monoPRM%show_time     = .true.
!    if(show_summary   == 1) monoPRM%show_summary  = .true.
!
!#ifdef TEST_ALL
!    call monolis_solve_test(monoPRM, monoCOM, monoMAT)
!#else
!    call monolis_solve(monoPRM, monoCOM, monoMAT)
!#endif
!    monoMAT%item = monoMAT%item - 1
!  end subroutine monolis_c
!
!  subroutine monolis_matvec_wrapper(monoCOM, N, NP, NZ, NDOF, A, index, item, X, Y, tcomm)
!    implicit none
!    type(monolis_com) :: monoCOM
!    type(monolis_mat) :: monoMAT
!    integer(kind=kint) :: N, NP, NZ, NDOF
!    integer(kind=kint), pointer :: index(:), item(:)
!    real(kind=kdouble), pointer :: A(:)
!    real(kind=kdouble) :: X(:), Y(:)
!    real(kind=kdouble), optional :: tcomm
!
!    !> for monoMAT
!    monoMAT%N = N
!    monoMAT%NP = NP
!    monoMAT%NZ = NZ
!    monoMAT%NDOF = NDOF
!    monoMAT%A => A
!    monoMAT%index => index
!    monoMAT%item => item
!
!    call monolis_matvec(monoCOM, monoMAT, X, Y, tcomm)
!  end subroutine monolis_matvec_wrapper
!
!  subroutine monolis_matvec_wrapper_c(monoCOM_c, N, NP, NZ, NDOF, A, index, item, X_c, Y_c) &
!    & bind(c, name="monolis_matvec_wrapper")
!    use iso_c_binding
!    implicit none
!    type(monolis_com_c) :: monoCOM_c
!    type(monolis_com) :: monoCOM
!    type(monolis_mat) :: monoMAT
!    !> for monoMAT
!    integer(c_int), value :: N, NP, NZ, NDOF
!    integer(c_int), intent(in), target :: index(0:NP)
!    integer(c_int), intent(in), target :: item(NZ)
!    real(c_double), intent(in), target :: A(NZ*NDOF*NDOF)
!    real(c_double), intent(in), target :: X_c(NP*NDOF)
!    real(c_double), intent(out),target :: Y_c(NP*NDOF)
!    real(kind=kdouble), pointer :: X(:), Y(:)
!    !> for monoMAT
!
!    monoMAT%N = N
!    monoMAT%NP = NP
!    monoMAT%NZ = NZ
!    monoMAT%NDOF = NDOF
!    monoMAT%A => A
!    monoMAT%index => index
!    monoMAT%item => item
!    monoMAT%item = monoMAT%item + 1
!    !> for monoCOM
!    monoCOM%myrank = monoCOM_c%myrank
!    monoCOM%comm = monoCOM_c%comm
!    monoCOM%commsize = monoCOM_c%commsize
!    monoCOM%n_neib = monoCOM_c%n_neib
!
!    X => X_c
!    Y => Y_c
!    call monolis_matvec(monoCOM, monoMAT, X, Y)
!    monoMAT%item = monoMAT%item - 1
!  end subroutine monolis_matvec_wrapper_c

end module mod_monolis_wrapper
