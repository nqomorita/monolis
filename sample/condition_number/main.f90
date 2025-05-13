program main
  use mod_monolis
  use mod_gedatsu
  implicit none
  integer(kint) :: iter, prec
  real(kdouble) :: condition_number_lanczos

  call monolis_global_initialize()

  call monolis_std_log_string("monolis_condition_number_R_test")

  !call monolis_solver_parallel_R_test()

  call monolis_condition_number_R_test()

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
    real(kdouble) :: val, res_conv
    character(monolis_charlen) :: fname
    integer(kint), allocatable :: elem(:,:), global_eid(:)
    real(kdouble), allocatable :: coef(:), node(:,:)
    real(kdouble), allocatable :: a(:), b(:), c(:)
    real(kdouble), allocatable :: eig_val1(:), eig_mode1(:,:)
    real(kdouble), allocatable :: eig_val2(:), eig_mode2(:,:)
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
        !call monolis_add_scalar_to_sparse_matrix_R(mat, eid(2), eid(1), 1, 1, val)
      endif
    enddo

    !> eigen test
    call monolis_std_log_string("monolis_solver_parallel_test eigen")

    n_get_eigen = 10

    call monolis_alloc_R_1d(eig_val1, n_get_eigen)
    call monolis_alloc_R_1d(eig_val2, n_get_eigen)
    call monolis_alloc_R_2d(eig_mode1, n_node, n_get_eigen)
    call monolis_alloc_R_2d(eig_mode2, n_node, n_get_eigen)
    call monolis_alloc_L_1d(is_bc, n_node)

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

    if(monolis_mpi_get_global_my_rank() == 0) write(*,*)"condition_number_lanczos", condition_number_lanczos

    call monolis_finalize(mat)
  end subroutine monolis_solver_parallel_R_test

  subroutine monolis_condition_number_R_test()
    implicit none
    type(monolis_structure) :: mat !> 疎行列変数
    type(monolis_com) :: com
    integer(kint) :: n_node, n_elem, n_base, n_id
    integer(kint) :: n_coef, eid(2), i, j, shift, id1, id2, ndof
    real(kdouble) :: val, condition_number
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
        !call monolis_add_scalar_to_sparse_matrix_R(mat, eid(2), eid(1), 1, 1, val)
      endif
    enddo

    !> monolis_convert_sparse_matrix_to_dense_matrix_R
    call monolis_convert_sparse_matrix_to_dense_matrix_R(mat%MAT, com, dense)

    !> monolis_get_condition_number_R
    call monolis_get_condition_number_R(mat, com, condition_number)

    if(monolis_mpi_get_global_my_rank() == 0) write(*,*)"condition_number        ", condition_number

    !call monolis_test_check_eq_R1("monolis_get_condition_number_R test", &
    !  & condition_number, condition_number_lanczos)

    call monolis_finalize(mat)
  end subroutine monolis_condition_number_R_test
end program main