program main
  use mod_monolis
  use mod_gedatsu
  implicit none
  integer(kint) :: iter, prec
  real(kdouble) :: condition_number_lanczos

  call monolis_global_initialize()

  call monolis_std_log_string("monolis_solver_parallel_arbit_dof_test")

  call monolis_solver_parallel_R_test()

  call monolis_solver_parallel_C_test()

  call monolis_condition_number_R_test()

  call monolis_global_finalize()

  contains

  subroutine monolis_condition_number_R_test()
    implicit none
    type(monolis_structure) :: mat !> 疎行列変数
    type(monolis_com) :: com
    integer(kint) :: n_node, n_elem, n_base, n_id
    integer(kint) :: n_coef, eid(2), i, j, shift, id1, id2, ndof
    real(kdouble) :: val, condition_number, rmax, rmin
    character(monolis_charlen) :: fname
    integer(kint), allocatable :: elem(:,:), global_eid(:), global_nid(:), vtxdist(:)
    real(kdouble), allocatable :: dense(:,:)
    real(kdouble), allocatable :: coef(:), node(:,:)

    call monolis_std_log_string("monolis_get_condition_number_R")

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.dat")
    call monolis_input_node(fname, n_node, node)

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "elem.dat")
    call monolis_input_elem(fname, n_elem, n_base, elem)

    if(monolis_mpi_get_global_comm_size() > 1)then
      fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.dat.id")
      call monolis_input_global_id(fname, n_id, global_nid)

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

    ndof = 1
    call monolis_get_nonzero_pattern_by_simple_mesh_R(mat, n_node, 2, ndof, n_elem, elem)

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

    !> monolis_convert_sparse_matrix_to_dense_matrix_R
    call monolis_convert_sparse_matrix_to_dense_matrix_R(mat%MAT, com, dense)

    if(monolis_mpi_get_global_comm_size() == 1)then
      do i = 1, 10
        call monolis_test_check_eq_R1("monolis_convert_sparse_matrix_to_dense_matrix_R test", dense(i,i), 4.0d0)
      enddo
      do i = 1, 9
        call monolis_test_check_eq_R1("monolis_convert_sparse_matrix_to_dense_matrix_R test", dense(i,i+1), 1.0d0)
        call monolis_test_check_eq_R1("monolis_convert_sparse_matrix_to_dense_matrix_R test", dense(i+1,i), 1.0d0)
      enddo

    elseif(monolis_mpi_get_global_comm_size() == 3)then
      call monolis_com_n_vertex_list(com%n_internal_vertex, com%comm, vtxdist)
      do i = 1, com%n_internal_vertex
        shift = vtxdist(monolis_mpi_get_global_my_rank() + 1)
        call monolis_test_check_eq_R1("monolis_convert_sparse_matrix_to_dense_matrix_R test", &
          & dense(i,shift + i), 4.0d0)

        id1 = global_nid(i)
        do j = 1, com%n_internal_vertex
          id2 = global_nid(j)
          if(id1 == id2 - 1) call monolis_test_check_eq_R1("monolis_convert_sparse_matrix_to_dense_matrix_R test", &
            & dense(i,shift + j), 1.0d0)
          if(id1 == id2 + 1) call monolis_test_check_eq_R1("monolis_convert_sparse_matrix_to_dense_matrix_R test", &
            & dense(i,shift + j), 1.0d0)
        enddo
      enddo
    endif

    !> monolis_get_condition_number_R
    call monolis_get_condition_number_R(mat, com, rmax, rmin)
    condition_number = rmax/rmin

    call monolis_test_check_eq_R1("monolis_get_condition_number_R test", &
      & condition_number, condition_number_lanczos)

    call monolis_finalize(mat)
  end subroutine monolis_condition_number_R_test

  subroutine monolis_solver_parallel_R_test()
    implicit none
    type(monolis_structure) :: mat !> 疎行列変数
    type(monolis_com) :: com
    integer(kint) :: n_node, n_elem, n_base, n_id
    integer(kint) :: n_coef, eid(2)
    integer(kint) :: i, j, iter_conv, n
    integer(kint) :: n_get_eigen
    real(kdouble) :: val, res_conv
    character(monolis_charlen) :: fname
    integer(kint), allocatable :: global_nid(:)
    integer(kint), allocatable :: elem(:,:), global_eid(:), n_dof_list(:)
    real(kdouble), allocatable :: coef(:), node(:,:)
    real(kdouble), allocatable :: a(:), b(:), c(:)
    real(kdouble), allocatable :: eig_val1(:), eig_mode1(:,:)
    real(kdouble), allocatable :: eig_val2(:), eig_mode2(:,:)
    logical, allocatable :: is_bc(:)

    call monolis_std_log_string("monolis_solver_parallel_R_test linear")

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.dat")
    call monolis_input_node(fname, n_node, node)

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "elem.dat")
    call monolis_input_elem(fname, n_elem, n_base, elem)

    if(monolis_mpi_get_global_comm_size() > 1)then
      fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.dat.id")
      call monolis_input_global_id(fname, n_id, global_nid)
    else
      call monolis_alloc_I_1d(global_nid, n_node)
      do i = 1, n_node
        global_nid(i) = i
      enddo
    endif

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

    call monolis_alloc_I_1d(n_dof_list, n_node)

    do i = 1, n_node
      j = global_nid(i)
      if(mod(j,2) /= 0)then
        n_dof_list(i) = 1
      else
        !n_dof_list(i) = 1
        n_dof_list(i) = 2
      endif
    enddo

    call monolis_get_nonzero_pattern_by_simple_mesh_R( &
      mat, n_node, 2, 1, n_elem, elem)

    !call monolis_get_nonzero_pattern_by_simple_mesh_V_R( &
    !  mat, n_node, 2, n_dof_list, n_elem, elem)
  
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
        if(mod(eid(1),2) /= 0)then
          call monolis_add_scalar_to_sparse_matrix_R(mat, eid(1), eid(2), 1, 1, val)
        else
          call monolis_add_scalar_to_sparse_matrix_R(mat, eid(1), eid(2), 1, 1, val)
          !call monolis_add_scalar_to_sparse_matrix_R(mat, eid(1), eid(2), 2, 2, val)
          !call monolis_add_scalar_to_sparse_matrix_R(mat, eid(1), eid(2), 1, 2, 0.25d0*val)
          !scall monolis_add_scalar_to_sparse_matrix_R(mat, eid(1), eid(2), 2, 1, 0.25d0*val)
        endif
      else
        call monolis_add_scalar_to_sparse_matrix_R(mat, eid(1), eid(2), 1, 1, val)
        call monolis_add_scalar_to_sparse_matrix_R(mat, eid(2), eid(1), 1, 1, val)
      endif
    enddo

    n = mat%MAT%n_dof_index(n_node + 1)

    call monolis_alloc_R_1d(a, n)
    call monolis_alloc_R_1d(b, n)
    call monolis_alloc_R_1d(c, n)

    a = 1.0d0

    call monolis_matvec_product_R(mat, com, a, c)

    call monolis_set_maxiter(mat, 1000)
    call monolis_set_tolerance(mat, 1.0d-10)
    call monolis_show_timelog(mat, .true.)
    call monolis_show_iterlog(mat, .true.)
    call monolis_show_summary(mat, .true.)

    do iter = monolis_iter_CG, monolis_iter_IDRS
    do prec = monolis_prec_NONE, monolis_prec_SOR
      a = 0.0d0
      b = c

      call monolis_set_method(mat, iter)
      call monolis_set_precond(mat, prec)

      call monolis_solve_R(mat, com, b, a)

      call monolis_mpi_global_barrier();

      b = 1.0d0

      if(monolis_mpi_get_global_my_rank() == 0)then
        write(*,*)"iter", iter, ", prec", prec
      endif
      call monolis_test_check_eq_R("monolis_solver_parallel_R_test", a, b)

      call monolis_get_converge_iter(mat, iter_conv)
      if(iter_conv <= 1)then
        call monolis_test_assert_fail("monolis_solver_parallel_R_test", "conv iter is less than 1")
      endif

      call monolis_get_converge_residual(mat, res_conv)

      if(res_conv > 1.0d-10)then
        call monolis_test_assert_fail("monolis_solver_parallel_R_test", "residual is greater than ths")
      endif

      call monolis_mpi_global_barrier();
    enddo
    enddo


    !> eigen test
    call monolis_std_log_string("monolis_solver_parallel_test eigen")

    n_get_eigen = 10

    call monolis_alloc_R_1d(eig_val1, n_get_eigen)
    call monolis_alloc_R_1d(eig_val2, n_get_eigen)
    call monolis_alloc_R_2d(eig_mode1, n, n_get_eigen)
    call monolis_alloc_R_2d(eig_mode2, n, n_get_eigen)
    call monolis_alloc_L_1d(is_bc, n)

    call monolis_set_method(mat, monolis_iter_CG)
    call monolis_set_precond(mat, monolis_prec_SOR)
    call monolis_show_timelog(mat, .false.)
    call monolis_show_iterlog(mat, .false.)
    call monolis_show_summary(mat, .false.)

    call monolis_eigen_standard_lanczos_R &
      & (mat, com, n_get_eigen, 1.0d-6, 100, eig_val1, eig_mode1, is_bc)

    call monolis_eigen_inverted_standard_lanczos_R &
      & (mat, com, n_get_eigen, 1.0d-6, 100, eig_val2, eig_mode2, is_bc)

    condition_number_lanczos = eig_val1(1)/eig_val1(10)

    do i = 1, n_get_eigen
      j = n_get_eigen - i + 1

      call monolis_test_check_eq_R1("monolis_solver_parallel_R_test eigen value", eig_val1(i), eig_val2(j))

      call monolis_test_check_eq_R ("monolis_solver_parallel_R_test eigen mode", &
        & dabs(eig_mode1(:,i)), dabs(eig_mode2(:,j)))
    enddo

    call monolis_finalize(mat)
  end subroutine monolis_solver_parallel_R_test

  subroutine monolis_solver_parallel_C_test()
    implicit none
    type(monolis_structure) :: mat !> 疎行列変数
    type(monolis_com) :: com
    integer(kint) :: n_node, n_elem, n_base, n_id, n
    integer(kint) :: n_coef, eid(2)
    integer(kint) :: i, j, iter_conv
    real(kdouble) :: r, res_conv
    complex(kdouble) :: val
    character(monolis_charlen) :: fname
    integer(kint), allocatable :: global_nid(:)
    integer(kint), allocatable :: elem(:,:), global_eid(:), n_dof_list(:)
    real(kdouble), allocatable :: coef(:), node(:,:)
    complex(kdouble), allocatable :: a(:), b(:), c(:)

    call monolis_std_log_string("monolis_solver_parallel_C_test linear")

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.dat")
    call monolis_input_node(fname, n_node, node)

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "elem.dat")
    call monolis_input_elem(fname, n_elem, n_base, elem)

    if(monolis_mpi_get_global_comm_size() > 1)then
      fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.dat.id")
      call monolis_input_global_id(fname, n_id, global_nid)
    else
      call monolis_alloc_I_1d(global_nid, n_node)
      do i = 1, n_node
        global_nid(i) = i
      enddo
    endif

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

    call monolis_alloc_I_1d(n_dof_list, n_node)

    do i = 1, n_node
      j = global_nid(i)
      if(mod(j,2) /= 0)then
        n_dof_list(i) = 1
      else
        n_dof_list(i) = 1
        !n_dof_list(i) = 2
      endif
    enddo

    call monolis_get_nonzero_pattern_by_simple_mesh_C( &
      & mat, n_node, 2, 1, n_elem, elem)

    !call monolis_get_nonzero_pattern_by_simple_mesh_V_C( &
    !  mat, n_node, 2, n_dof_list, n_elem, elem)

    open(20, file = "coef.dat", status = "old")
      read(20,*) n_coef
      call monolis_alloc_R_1d(coef, n_coef)
      do i = 1, n_coef
        read(20,*) coef(i)
      enddo
    close(20)

    do i = 1, n_elem
      eid = elem(:,i)
      r = coef(global_eid(i))
      val = complex(r, r)
      if(eid(1) == eid(2))then
        if(mod(eid(1),2) /= 0)then
          call monolis_add_scalar_to_sparse_matrix_C(mat, eid(1), eid(2), 1, 1, val)
        else
          call monolis_add_scalar_to_sparse_matrix_C(mat, eid(1), eid(2), 1, 1, val)
          !call monolis_add_scalar_to_sparse_matrix_C(mat, eid(1), eid(2), 2, 2, val)
          !call monolis_add_scalar_to_sparse_matrix_C(mat, eid(1), eid(2), 1, 2, 0.25d0*val)
          !call monolis_add_scalar_to_sparse_matrix_C(mat, eid(1), eid(2), 2, 1, 0.25d0*val)
        endif
      else
        call monolis_add_scalar_to_sparse_matrix_C(mat, eid(1), eid(2), 1, 1, val)
        call monolis_add_scalar_to_sparse_matrix_C(mat, eid(2), eid(1), 1, 1, val)
      endif
    enddo

    n = mat%MAT%n_dof_index(n_node + 1)

    call monolis_alloc_C_1d(a, n)
    call monolis_alloc_C_1d(b, n)
    call monolis_alloc_C_1d(c, n)

    a = (1.0d0, 1.0d0)

    call monolis_matvec_product_C(mat, com, a, c)

    call monolis_set_maxiter(mat, 1000)
    call monolis_set_tolerance(mat, 1.0d-10)
    call monolis_show_timelog(mat, .true.)
    call monolis_show_iterlog(mat, .true.)
    call monolis_show_summary(mat, .true.)

    do iter = monolis_iter_COCG, monolis_iter_COCG
    do prec = monolis_prec_NONE, monolis_prec_SOR
      a = (0.0d0, 0.0d0)
      b = c

      call monolis_set_method(mat, iter)
      call monolis_set_precond(mat, prec)

      call monolis_solve_C(mat, com, b, a)

      call monolis_mpi_global_barrier();

      b = (1.0d0, 1.0d0)

      call monolis_test_check_eq_C("monolis_solver_parallel_C_test", a, b)

      call monolis_get_converge_iter(mat, iter_conv)
      if(iter_conv <= 1)then
        call monolis_test_assert_fail("monolis_solver_parallel_R_test", "conv iter is less than 1")
      endif

      call monolis_get_converge_residual(mat, res_conv)

      if(res_conv > 1.0d-10)then
        call monolis_test_assert_fail("monolis_solver_parallel_R_test", "residual is greater than ths")
      endif

      call monolis_mpi_global_barrier();
    enddo
    enddo

    call monolis_finalize(mat)
  end subroutine monolis_solver_parallel_C_test
end program main
