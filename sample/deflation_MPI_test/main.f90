program main
  use mod_monolis
  use mod_gedatsu
  implicit none

  call monolis_global_initialize()

  call monolis_std_log_string("monolis_solver_parallel_test")

  call monolis_solver_parallel_R_test()

  call monolis_global_finalize()

  contains

  subroutine monolis_solver_parallel_R_test()
    implicit none
    type(monolis_structure) :: mat !> 疎行列変数
    type(monolis_com) :: com
    integer(kint) :: n_node, n_elem, n_base, n_id
    integer(kint) :: n_coef, eid(2)
    integer(kint) :: i, j, iter_conv
    integer(kint) :: n_get_eigen
    integer(kint) :: n_deflation_mode
    real(kdouble) :: val, res_conv
    character(monolis_charlen) :: fname
    integer(kint), allocatable :: elem(:,:), global_eid(:)
    real(kdouble), allocatable :: coef(:), node(:,:)
    real(kdouble), allocatable :: a(:), b(:), c(:)
    real(kdouble), allocatable :: eig_val(:), deflation_mode(:,:)
    logical, allocatable :: is_bc(:)

    call monolis_std_log_string("monolis_solver_parallel_test linear")

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.dat")
    call monolis_input_node(fname, n_node, node)

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "elem.dat")
    call monolis_input_elem(fname, n_elem, n_base, elem)

    if(monolis_mpi_get_global_comm_size() > 1)then
      fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "elem.dat.id")
      call monolis_input_global_id(fname, n_id, global_eid)
    else
      call monolis_alloc_I_1d(global_eid, n_elem)
      do i = 1, n_elem
        global_eid(i) = i
      enddo
    endif

    call monolis_initialize(mat)

    call monolis_com_initialize_by_parted_files(com, &
      monolis_mpi_get_global_comm(), &
      MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.dat")

    call monolis_get_nonzero_pattern_by_simple_mesh_R(mat, n_node, 2, 1, n_elem, elem)

    open(20, file = "coef.dat", status = "old")
      read(20,*) n_coef
      call monolis_alloc_R_1d(coef, n_coef)
      do i = 1, n_coef
        read(20,*) coef(i)
      enddo
    close(20)

    do i = 1, n_elem
      eid = elem(:,i)
      val = coef(global_eid(i))
      if(eid(1) == eid(2))then
        call monolis_add_scalar_to_sparse_matrix_R(mat, eid(1), eid(2), 1, 1, val)
      else
        call monolis_add_scalar_to_sparse_matrix_R(mat, eid(1), eid(2), 1, 1, val)
        call monolis_add_scalar_to_sparse_matrix_R(mat, eid(2), eid(1), 1, 1, val)
      endif
    enddo

    call monolis_alloc_R_1d(a, n_node)
    call monolis_alloc_R_1d(b, n_node)
    call monolis_alloc_R_1d(c, n_node)

    a = 1.0d0
    call monolis_matvec_product_R(mat, com, a, c)

    !> get global deflation mode
    n_deflation_mode = 2
    allocate(deflation_mode(n_node,n_deflation_mode), source = 0.0d0)
    allocate(eig_val(n_deflation_mode), source = 0.0d0)
    allocate(is_bc(n_node), source = .false.)

    call monolis_set_method(mat, monolis_iter_CG)
    call monolis_set_precond(mat, monolis_prec_SOR)
    call monolis_set_tolerance(mat, 1.0d-8)
    call monolis_show_iterlog(mat, .false.)
    call monolis_show_timelog(mat, .false.)
    call monolis_show_summary(mat, .false.)

    is_bc = .false.
    call monolis_eigen_inverted_standard_lanczos_R &
      & (mat, com, n_deflation_mode, 1.0d-10, 1000, eig_val, deflation_mode, is_bc)

    call monolis_set_deflation_mode(mat, n_deflation_mode, deflation_mode)

    !> solve deflation CG
    !call monolis_set_method(mat, monolis_iter_CG)
    call monolis_set_method(mat, monolis_iter_DeflatedCG1)
    !call monolis_set_precond(mat, monolis_prec_None)
    call monolis_set_precond(mat, monolis_prec_Diag)
    !call monolis_set_precond(mat, monolis_prec_SOR)
    call monolis_set_maxiter(mat, 5000)
    call monolis_set_tolerance(mat, 1.0d-10)
    call monolis_show_timelog(mat, .false.)
    call monolis_show_iterlog(mat, .true.)
    call monolis_show_summary(mat, .true.)

    a = 0.0d0
    call monolis_solve_R(mat, com, c, a)

    call monolis_mpi_global_barrier()

    b = 1.0d0

    call monolis_test_check_eq_R("monolis_solver_parallel_R_test", a, b)

    call monolis_get_converge_iter(mat, iter_conv)
    if(iter_conv <= 1)then
      call monolis_test_assert_fail("monolis_solver_parallel_R_test", "conv iter is less than 1")
    endif

    call monolis_get_converge_residual(mat, res_conv)

    if(res_conv > 1.0d-10)then
      call monolis_test_assert_fail("monolis_solver_parallel_R_test", "residual is greater than ths")
    endif

    call monolis_mpi_global_barrier()

    call monolis_finalize(mat)
  end subroutine monolis_solver_parallel_R_test
end program main
