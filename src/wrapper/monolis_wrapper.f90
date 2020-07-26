module mod_monolis_wrapper
  use iso_c_binding
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_util
  use mod_monolis_sparse_util
  use mod_monolis_stdlib
  use mod_monolis_solve
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
  subroutine monolis_get_CRR_format_c(N, NZ, index, item, indexR, itemR, permR) &
    & bind(c, name = "monolis_get_CRR_format")
    implicit none
    integer(c_int), intent(in), value :: N, NZ
    integer(c_int), intent(in), target :: index(0:N)
    integer(c_int), intent(in), target :: item(NZ)
    integer(c_int), target :: indexR(0:N)
    integer(c_int), target :: itemR(NZ)
    integer(c_int), target :: permR(NZ)
    integer(kint), pointer :: indexRt(:)
    integer(kint), pointer :: itemRt(:)
    integer(kint), pointer :: permRt(:)

    indexRt => indexR
    itemRt => itemR
    permRt => permR
    call monolis_get_CRR_format(N, NZ, index, item, indexRt, itemRt, permRt)
  end subroutine monolis_get_CRR_format_c

  subroutine monolis_add_sparse_matrix_c(N, NZ, NDOF, NBF, index, item, A, conn, mat) &
    & bind(c, name = "monolis_add_sparse_matrix_c_main")
    implicit none
    integer(c_int), intent(in), value :: N, NZ, NDOF, NBF
    integer(c_int), intent(in), target :: index(0:N)
    integer(c_int), intent(in), target :: item(NZ)
    integer(c_int), intent(in), target :: conn(NBF)
    real(c_double), target :: A(NDOF*NDOF*NZ)
    real(c_double), target :: mat(NDOF*NDOF*NBF*NBF)
    integer(kint) :: conn_t(NBF)

    conn_t = conn + 1
    call monolis_sparse_matrix_assemble(index, item, A, NBF, NDOF, conn_t, conn_t, mat)
  end subroutine monolis_add_sparse_matrix_c

  subroutine monolis_set_Dirichlet_bc_c(N, NZ, NDOF, index, item, indexR, itemR, permR, &
    & A, B, nid, ndof_bc, val) &
    & bind(c, name = "monolis_set_Dirichlet_bc_c_main")
    implicit none
    integer(c_int), intent(in), value :: N, NZ, NDOF, nid, ndof_bc
    integer(c_int), intent(in), target :: index(0:N)
    integer(c_int), intent(in), target :: item(NZ)
    integer(c_int), intent(in), target :: indexR(0:N)
    integer(c_int), intent(in), target :: itemR(NZ)
    integer(c_int), intent(in), target :: permR(NZ)
    real(c_double), intent(in), value :: val
    real(c_double), target :: A(NDOF*NDOF*NZ)
    real(c_double), target :: B(NDOF*N)

    call monolis_sparse_matrix_add_bc(index, item, A, B, indexR, itemR, permR, &
      & ndof, nid, ndof_bc, val)
  end subroutine monolis_set_Dirichlet_bc_c

  !> solve
  subroutine monolis_solve_c(N, NP, NZ, NDOF, A, X, B, index, item, &
    myrank, comm, commsize, recv_n_neib, send_n_neib, &
    method, precond, maxiter, tol) &
    & bind(c, name = "monolis_solve_c_main")
    implicit none
    type(monolis_structure) :: monolis
    integer(c_int), intent(in), value :: N, NP, NZ, NDOF
    integer(c_int), intent(in), value :: myrank, comm, commsize, recv_n_neib, send_n_neib
    integer(c_int), intent(in), value :: method, precond, maxiter
    integer(c_int), intent(in), target :: index(0:N)
    integer(c_int), intent(in), target :: item(NZ)
    real(c_double), intent(in), value :: tol
    real(c_double), intent(in), target :: A(NDOF*NDOF*NZ)
    real(c_double), intent(in), target :: X(NDOF*N)
    real(c_double), intent(in), target :: B(NDOF*N)

    !> for monoMAT
    monolis%MAT%N = N
    monolis%MAT%NP = NP
    monolis%MAT%NZ = NZ
    monolis%MAT%NDOF = NDOF
    monolis%MAT%A => A
    monolis%MAT%X => X
    monolis%MAT%B => B
    monolis%MAT%index => index
    monolis%MAT%item => item

    !> for monoCOM
    monoliS%COM%myrank = myrank
    monoliS%COM%comm = comm
    monoliS%COM%commsize = commsize
    monoliS%COM%recv_n_neib = recv_n_neib
    monoliS%COM%send_n_neib = send_n_neib

    !> for monoPRM
    monolis%PRM%method = method
    monolis%PRM%precond = precond
    monolis%PRM%maxiter = maxiter
    monolis%PRM%tol = tol
    monolis%PRM%is_scaling    = .false.
    monolis%PRM%is_reordering = .false.
    monolis%PRM%is_init_x     = .true.
    monolis%PRM%is_debug      = .true.
    monolis%PRM%show_iterlog  = .true.
    monolis%PRM%show_time     = .true.
    monolis%PRM%show_summary  = .true.

    call monolis_solve(monolis, B, X)
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

end module mod_monolis_wrapper
