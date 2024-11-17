program main
  use mod_monolis
  use mod_gedatsu
  implicit none

  call monolis_global_initialize()

  call monolis_std_log_string("monolis_conv_symmetric_test")

  call monolis_solver_parallel_R_test()

  call monolis_global_finalize()

  contains

  subroutine monolis_solver_parallel_R_test()
    implicit none
    type(monolis_structure) :: mat !> 疎行列変数
    type(monolis_com) :: com
    integer(kint) :: n_node, n_elem, n_base, n_id, eid(2), nid(2)
    integer(kint) :: i, j, in, jS, jE
    real(kdouble) :: val, ans
    character(monolis_charlen) :: fname
    integer(kint), allocatable :: elem(:,:), global_nid(:)
    real(kdouble), allocatable :: node(:,:)

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.dat")
    call monolis_input_node(fname, n_node, node)

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "elem.dat")
    call monolis_input_elem(fname, n_elem, n_base, elem)

    if(monolis_mpi_get_global_comm_size() > 1)then
      fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.dat.id")
      call monolis_input_global_id(fname, n_id, global_nid)
    else
      call monolis_alloc_I_1d(global_nid, n_elem)
      do i = 1, n_elem
        global_nid(i) = i
      enddo
    endif

    call monolis_initialize(mat)

    call monolis_com_initialize_by_parted_files(com, &
      monolis_mpi_get_global_comm(), &
      MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.dat")

    call monolis_get_nonzero_pattern_by_simple_mesh_R(mat, n_node, 2, 2, n_elem, elem)

    do i = 1, n_elem
      eid = elem(:,i)
      nid(1) = global_nid(elem(1,i))
      nid(2) = global_nid(elem(2,i))
      if(eid(1) == eid(2))then
        val = 10000.0d0*nid(1) + 100.0d0*nid(2)
        call monolis_add_scalar_to_sparse_matrix_R(mat, eid(1), eid(2), 1, 1, val + 1.0d0)
        call monolis_add_scalar_to_sparse_matrix_R(mat, eid(1), eid(2), 1, 2, val + 2.0d0)
        call monolis_add_scalar_to_sparse_matrix_R(mat, eid(1), eid(2), 2, 1, val + 2.0d0)
        call monolis_add_scalar_to_sparse_matrix_R(mat, eid(1), eid(2), 2, 2, val + 3.0d0)
      elseif(eid(1) < eid(2))then
        val = 10000.0d0*nid(1) + 100.0d0*nid(2)
        call monolis_add_scalar_to_sparse_matrix_R(mat, eid(1), eid(2), 1, 1, val + 1.0d0)
        call monolis_add_scalar_to_sparse_matrix_R(mat, eid(1), eid(2), 1, 2, val + 2.0d0)
        call monolis_add_scalar_to_sparse_matrix_R(mat, eid(1), eid(2), 2, 1, val + 3.0d0)
        call monolis_add_scalar_to_sparse_matrix_R(mat, eid(1), eid(2), 2, 2, val + 4.0d0)
        val = 10000.0d0*nid(2) + 100.0d0*nid(1)
        call monolis_add_scalar_to_sparse_matrix_R(mat, eid(2), eid(1), 1, 1, val + 1.0d0)
        call monolis_add_scalar_to_sparse_matrix_R(mat, eid(2), eid(1), 1, 2, val + 3.0d0)
        call monolis_add_scalar_to_sparse_matrix_R(mat, eid(2), eid(1), 2, 1, val + 2.0d0)
        call monolis_add_scalar_to_sparse_matrix_R(mat, eid(2), eid(1), 2, 2, val + 4.0d0)
      endif
    enddo

    call monolis_std_log_string("monolis_matrix_convert_to_symmetric_inner_R")
    call monolis_matrix_convert_to_symmetric_inner_R(mat%MAT)

    call monolis_std_log_string("monolis_matrix_convert_to_symmetric_outer_R")
    call monolis_matrix_convert_to_symmetric_outer_R(mat%MAT, com)

    do i = 1, mat%MAT%N
      jS = mat%MAT%CSR%index(i) + 1
      jE = mat%MAT%CSR%index(i + 1)
      do j = jS, jE
        in = mat%MAT%CSR%item(j)
        nid(1) = global_nid(i)
        nid(2) = global_nid(in)
        if(mat%MAT%N < in)then
          !> no update region
            ans = 10000.0d0*nid(1) + 100.0d0*nid(2) + 1.0d0
            if(ans /= mat%MAT%R%A(4*j-3)) stop "A1"
            ans = 10000.0d0*nid(1) + 100.0d0*nid(2) + 2.0d0
            if(ans /= mat%MAT%R%A(4*j-2)) stop "A2"
            ans = 10000.0d0*nid(1) + 100.0d0*nid(2) + 3.0d0
            if(ans /= mat%MAT%R%A(4*j-1)) stop "A3"
            ans = 10000.0d0*nid(1) + 100.0d0*nid(2) + 4.0d0
            if(ans /= mat%MAT%R%A(4*j  )) stop "A4"
        else
          !> update region
          if(i == in)then
            ans = 5050.0d0*nid(1) + 5050.0d0*nid(2) + 1.0d0
            if(ans /= mat%MAT%R%A(4*j-3)) stop "B1"
            ans = 5050.0d0*nid(1) + 5050.0d0*nid(2) + 2.0d0
            if(ans /= mat%MAT%R%A(4*j-2)) stop "B2"
            if(ans /= mat%MAT%R%A(4*j-1)) stop "B3"
            ans = 5050.0d0*nid(1) + 5050.0d0*nid(2) + 3.0d0
            if(ans /= mat%MAT%R%A(4*j  )) stop "B4"
          elseif(i < in)then
            ans = 5050.0d0*nid(1) + 5050.0d0*nid(2) + 1.0d0
            if(ans /= mat%MAT%R%A(4*j-3)) stop "C1"
            ans = 5050.0d0*nid(1) + 5050.0d0*nid(2) + 2.0d0
            if(ans /= mat%MAT%R%A(4*j-2)) stop "C2"
            ans = 5050.0d0*nid(1) + 5050.0d0*nid(2) + 3.0d0
            if(ans /= mat%MAT%R%A(4*j-1)) stop "C3"
            ans = 5050.0d0*nid(1) + 5050.0d0*nid(2) + 4.0d0
            if(ans /= mat%MAT%R%A(4*j  )) stop "C4"
          elseif(in < i)then
            ans = 5050.0d0*nid(1) + 5050.0d0*nid(2) + 1.0d0
            if(ans /= mat%MAT%R%A(4*j-3)) stop "D1"
            ans = 5050.0d0*nid(1) + 5050.0d0*nid(2) + 3.0d0
            if(ans /= mat%MAT%R%A(4*j-2)) stop "D2"
            ans = 5050.0d0*nid(1) + 5050.0d0*nid(2) + 2.0d0
            if(ans /= mat%MAT%R%A(4*j-1)) stop "D3"
            ans = 5050.0d0*nid(1) + 5050.0d0*nid(2) + 4.0d0
            if(ans /= mat%MAT%R%A(4*j  )) stop "D4"
          endif
        endif
      enddo
    enddo

    call monolis_std_log_string("monolis_matrix_convert_to_symmetric_inner_R PASS")
    call monolis_std_log_string("monolis_matrix_convert_to_symmetric_outer_R PASS")

    call monolis_finalize(mat)
  end subroutine monolis_solver_parallel_R_test

end program main
