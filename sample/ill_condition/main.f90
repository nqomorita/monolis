program main
  use mod_monolis
  use mod_gedatsu
  implicit none

  call monolis_global_initialize()

  call monolis_std_log_string("monolis_ill_condition_test")

  call monolis_solver_parallel_R_test()

  call monolis_global_finalize()

  contains

  subroutine monolis_solver_parallel_R_test()
    implicit none
    type(monolis_structure) :: mat !> 疎行列変数
    type(monolis_com) :: com
    integer(kint) :: n_node, n_elem, n_base, n_id, eid(2), nid(2)
    integer(kint) :: i, j, in, jS, jE, n_coef, idx
    real(kdouble) :: val, ans
    character(monolis_charlen) :: fname
    integer(kint), allocatable :: elem(:,:), global_nid(:), global_eid(:), coef_dof(:,:), perm(:)
    real(kdouble), allocatable :: a(:), b(:), c(:)
    real(kdouble), allocatable :: node(:,:), coef(:)

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.f.dat")
    call monolis_input_node(fname, n_node, node)

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "elem.f.dat")
    call monolis_input_elem(fname, n_elem, n_base, elem)

    if(monolis_mpi_get_global_comm_size() > 1)then
      fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.f.dat.id")
      call monolis_input_global_id(fname, n_id, global_nid)

      fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "elem.f.dat.id")
      call monolis_input_global_id(fname, n_id, global_eid)
    else
      call monolis_alloc_I_1d(global_nid, n_node)
      do i = 1, n_node
        global_nid(i) = i
      enddo

      call monolis_alloc_I_1d(global_eid, n_elem)
      do i = 1, n_elem
        global_eid(i) = i
      enddo
    endif

    call monolis_initialize(mat)

    if(monolis_mpi_get_global_comm_size() > 1)then
      call monolis_com_initialize_by_parted_files(com, &
        monolis_mpi_get_global_comm(), &
        MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.f.dat")
    else
      call monolis_com_initialize_by_self(com)
    endif

    call monolis_get_nonzero_pattern_by_simple_mesh_R(mat, n_node, 2, 1, n_elem, elem)

    open(20, file = "coef.f.dat", status = "old")
      read(20,*) n_coef
      call monolis_alloc_R_1d(coef, n_coef)
      call monolis_alloc_I_2d(coef_dof, 2, n_coef)
      do i = 1, n_coef
        read(20,*) coef_dof(1,i), coef_dof(2,i), coef(i)
      enddo
    close(20)

    call monolis_alloc_I_1d(perm, n_node)
    do i = 1, n_node
      perm(i) = i
    enddo
    call monolis_qsort_I_2d(global_nid, perm, 1, n_node)

    do i = 1, n_elem
      in = global_eid(i)
      call monolis_bsearch_I(global_nid, 1, n_node, coef_dof(1,in), idx)
      eid(1) = perm(idx)

      call monolis_bsearch_I(global_nid, 1, n_node, coef_dof(2,in), idx)
      eid(2) = perm(idx)

      val = coef(in)
!write(*,*)eid(1), eid(2), val
      call monolis_add_scalar_to_sparse_matrix_R(mat, eid(1), eid(2), 1, 1, val)
    enddo

    call monolis_alloc_R_1d(a, n_node)
    call monolis_alloc_R_1d(b, n_node)
    call monolis_alloc_R_1d(c, n_node)

    a = 1.0d0

    call monolis_matvec_product_R(mat, com, a, b)

    call monolis_set_maxiter(mat, 1000)
    call monolis_set_tolerance(mat, 1.0d-10)
    call monolis_show_timelog(mat, .true.)
    call monolis_show_iterlog(mat, .true.)
    call monolis_show_summary(mat, .true.)

    call monolis_set_method(mat, monolis_iter_BiCGSTAB)
    call monolis_set_precond(mat, monolis_prec_DIAG)

    call monolis_solve_R(mat, com, b, c)

    b = 1.0d0
    call monolis_test_check_eq_R("monolis_ill_condition_test", c, b)

    call monolis_finalize(mat)
  end subroutine monolis_solver_parallel_R_test

end program main
