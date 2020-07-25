module mod_monolis_wrapper
  use iso_c_binding
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_util

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

  subroutine monolis_initialize_c( &
      prm_method, prm_precond, prm_curiter, prm_maxiter, prm_ierr, &
      prm_tol, prm_curresid, &
      !prm_is_scaling, prm_is_reordering, prm_is_init_x, prm_is_debug, &
      !prm_show_iterlog, prm_show_time, prm_show_summary
      N, NP, NZ, NDOF, &
      myrank, comm, commsize, recv_n_neib, send_n_neib) &
      bind(c, name = "monolis_initialize_c_main")
    implicit none
    type(monolis_structure) :: monolis
    !> prm
    integer(c_int), value :: prm_method
    integer(c_int), value :: prm_precond
    integer(c_int), value :: prm_curiter
    integer(c_int), value :: prm_maxiter
    integer(c_int), value :: prm_ierr
    real(c_double), value :: prm_tol
    real(c_double), value :: prm_curresid
    logical :: prm_is_scaling
    logical :: prm_is_reordering
    logical :: prm_is_init_x
    logical :: prm_is_debug
    logical :: prm_show_iterlog
    logical :: prm_show_time
    logical :: prm_show_summary
    !> mat
    integer(c_int), value :: N
    integer(c_int), value :: NP
    integer(c_int), value :: NZ
    integer(c_int), value :: NDOF
    !> comm
    integer(c_int), value :: myrank
    integer(c_int), value :: comm
    integer(c_int), value :: commsize
    integer(c_int), value :: recv_n_neib
    integer(c_int), value :: send_n_neib

    call monolis_initialize(monolis)

    !> prm
    prm_method = monolis%PRM%method
    prm_precond = monolis%PRM%precond
    prm_curiter = monolis%PRM%curiter
    prm_maxiter = monolis%PRM%maxiter
    prm_ierr = monolis%PRM%ierr
    prm_tol = monolis%PRM%tol
    prm_curresid = monolis%PRM%curresid
    !prm_is_scaling = monolis%PRM%is_scaling
    !prm_is_reordering = monolis%PRM%is_reordering
    !prm_is_init_x = monolis%PRM%is_init_x
    !prm_is_debug = monolis%PRM%is_debug
    !prm_show_iterlog = monolis%PRM%show_iterlog
    !prm_show_time = monolis%PRM%show_time
    !prm_show_summary = monolis%PRM%show_summary

    !> mat
    N = monolis%MAT%N
    NP = monolis%MAT%NP
    NZ = monolis%MAT%NZ
    NDOF = monolis%MAT%NDOF

    !> comm
    myrank = monolis%COM%myrank
    comm = monolis%COM%comm
    commsize = monolis%COM%commsize
    recv_n_neib = monolis%COM%recv_n_neib
    send_n_neib = monolis%COM%send_n_neib
  end subroutine monolis_initialize_c

  subroutine monolis_finalize_c() &
    & bind(c, name = "monolis_finalize_c_main")
    implicit none
    type(monolis_structure) :: monolis

  end subroutine monolis_finalize_c

  !> mat
  subroutine monolis_get_nonzero_pattern_c( &
      nnode_c, nbase_func_c, ndof_c, nelem_c, &
      elem_c, index_c, item_c, A_c, B_c, X_c) &
      bind(c, name = "monolis_get_nonzero_pattern_c_main")
    use mod_monolis_sparse_util
    implicit none
    type(monolis_structure) :: monolis
    integer(c_int), value :: nnode_c, nbase_func_c, ndof_c, nelem_c
    !integer(kint) :: nnode, nbase_func, ndof, nelem
    !integer(kint) :: test(4,4)
    type(c_ptr) :: elem_c
    type(c_ptr) :: index_c, item_c
    type(c_ptr) :: A_c, B_c, X_c
    integer(c_int), pointer :: elem(:,:)
    integer(c_int), pointer :: index(:), item(:)
    real(c_double), pointer :: A(:), B(:), X(:)

    !write(*,*)nnode_c, nbase_func_c, ndof_c, nelem_c
    !call c_f_pointer(elem_c, elem, [nbase_func_c,nnode_c])
    !write(*,*)elem

    !call monolis_get_nonzero_pattern(monolis, nnode_c, nbase_func_c, ndof_c, nelem_c, test)
  end subroutine monolis_get_nonzero_pattern_c

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
