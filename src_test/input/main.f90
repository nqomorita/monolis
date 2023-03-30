program main
  use mod_monolis
  use mod_gedatsu
  implicit none
  integer(kint) :: iter

  call monolis_global_initialize()

  call monolis_std_log_string("monolis_solver_parallel_test")

  call monolis_solver_parallel_test()

  call monolis_global_finalize()

  contains

  subroutine monolis_solver_parallel_test()
    implicit none
    type(monolis_structure) :: mat !> 疎行列変数
    integer(kint) :: n_node, n_elem, n_base, n_id
    integer(kint) :: n_coef, eid(2)
    integer(kint) :: i
    real(kdouble) :: val
    character(monolis_charlen) :: fname
    integer(kint), allocatable :: elem(:,:), global_eid(:)
    real(kdouble), allocatable :: coef(:), node(:,:)
    real(kdouble), allocatable :: a(:), b(:)

    fname = monolis_get_global_input_file_name("parted.0", "node.dat")
    call monolis_input_node(fname, n_node, node)

    fname = monolis_get_global_input_file_name("parted.0", "elem.dat")
    call monolis_input_elem(fname, n_elem, n_base, elem)

    if(monolis_mpi_get_global_comm_size() > 1)then
      fname = monolis_get_global_input_file_name("parted.0", "elem.dat.id")
      call monolis_input_global_id(fname, n_id, global_eid)
    else
      call monolis_alloc_I_1d(global_eid, n_elem)
      do i = 1, n_elem
        global_eid(i) = i
      enddo
    endif

    call monolis_initialize(mat)

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

    a = 1.0d0

    call monolis_matvec_product_R(mat, a, b)

    call monolis_set_method(mat, monolis_iter_CG)
    call monolis_set_precond(mat, monolis_prec_NONE)
    call monolis_set_maxiter(mat, 1000)
    call monolis_set_tolerance(mat, 1.0d-8)
    call monolis_show_timelog(mat, .true.)
    call monolis_show_iterlog(mat, .true.)
    call monolis_show_summary(mat, .true.)

    a = 0.0d0

    call monolis_solve_R(mat, b, a)

    b = 1.0d0

    call monolis_test_check_eq_R("monolis_solver_parallel_test", a, b)

    call monolis_finalize(mat)
  end subroutine monolis_solver_parallel_test
end program main
