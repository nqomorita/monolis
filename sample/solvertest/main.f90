program main
  use mod_monolis
  implicit none

  call monolis_global_initialize()

  call monolis_solver_COCG_test()

  call monolis_global_finalize()

  contains

  subroutine monolis_solver_COCG_test()
    implicit none
    integer(kint) :: nnode, nelem, nbase_func, i, id, n_coef, eid(2), my_rank, ierr
    real(kdouble) :: tspmv, tcomm
    real(kdouble) :: val
    integer(kint), allocatable :: elem(:,:), global_eid(:)
    real(kdouble), allocatable :: coef(:), node(:,:)
    real(kdouble), allocatable :: a(:), b(:), b_th(:)
    type(monolis_structure) :: mat !> 疎行列変数
    character :: fname*128

    call monolis_initialize(mat, "./") !> 疎行列変数の初期化

    fname = monolis_get_input_filename("node.dat")
    call monolis_input_mesh_node(fname, nnode, node)

    fname = monolis_get_input_filename("elem.dat")
    call monolis_input_mesh_elem(fname, nelem, nbase_func, elem)

    if(monolis_global_commsize() > 1)then
      fname = monolis_get_input_filename("connectivity.id")
      call monolis_input_id(fname, global_eid)
    else
      allocate(global_eid(nelem), source = 0)
      do i = 1, nelem
        global_eid(i) = i
      enddo
    endif

    call monolis_get_nonzero_pattern(mat, nnode, 2, 1, nelem, elem)

    open(20, file = "coef.dat", status = "old")
      read(20,*) n_coef
      allocate(coef(n_coef), source = 0.0d0)
      do i = 1, n_coef
        read(20,*) coef(i)
      enddo
    close(20)

    do i = 1, nelem
      eid = elem(:,i)
      call get_coef_val(n_coef, coef, global_eid(i), val)
      if(eid(1) == eid(2))then
        call monolis_add_scalar_to_sparse_matrix(mat, eid(1), eid(2), 1, 1, val)
      else
        call monolis_add_scalar_to_sparse_matrix(mat, eid(1), eid(2), 1, 1, val)
        call monolis_add_scalar_to_sparse_matrix(mat, eid(2), eid(1), 1, 1, val)
      endif
    enddo

    allocate(a(nnode), source = 0.0d0)
    allocate(b(nnode), source = 0.0d0)
    allocate(b_th(nnode), source = 0.0d0)
    b_th = (1.0d0, 1.0d0)

    call monolis_matvec_product(mat, b_th, b)

    !> get comm free
!    write(100+monolis_global_myrank(),*)"mat%MAT%N", mat%MAT%N, mat%MAT%NP, mat%COM%internal_nnode
!    write(100+monolis_global_myrank(),*)"global_node_id"
!    write(100+monolis_global_myrank(),*)mat%COM%global_node_id

!    write(100+monolis_global_myrank(),*)"BEFORE"
!    write(100+monolis_global_myrank(),*)"recv_n_neib_pe", mat%com%recv_n_neib
!    write(100+monolis_global_myrank(),*)"recv_neib_pe", mat%com%recv_neib_pe
!    write(100+monolis_global_myrank(),*)"recv_index", mat%com%recv_index
!    write(100+monolis_global_myrank(),*)"recv_item", mat%com%recv_item
!    write(100+monolis_global_myrank(),*)"send_n_neib_pe", mat%com%send_n_neib
!    write(100+monolis_global_myrank(),*)"send_neib_pe", mat%com%send_neib_pe
!    write(100+monolis_global_myrank(),*)"send_index", mat%com%send_index
!    write(100+monolis_global_myrank(),*)"send_item", mat%com%send_item

    deallocate(mat%com%recv_neib_pe)
    deallocate(mat%com%recv_index)
    deallocate(mat%com%recv_item)
    deallocate(mat%com%send_neib_pe)
    deallocate(mat%com%send_index)
    deallocate(mat%com%send_item)

    !> get comm
    call monolis_com_get_comm_table(mat, mat%COM%internal_nnode, mat%MAT%NP, mat%COM%global_node_id)

!    write(100+monolis_global_myrank(),*)"AFTER"
!    write(100+monolis_global_myrank(),*)"recv_n_neib_pe", mat%com%recv_n_neib
!    write(100+monolis_global_myrank(),*)"recv_neib_pe", mat%com%recv_neib_pe
!    write(100+monolis_global_myrank(),*)"recv_index", mat%com%recv_index
!    write(100+monolis_global_myrank(),*)"recv_item", mat%com%recv_item
!    write(100+monolis_global_myrank(),*)"send_n_neib_pe", mat%com%send_n_neib
!    write(100+monolis_global_myrank(),*)"send_neib_pe", mat%com%send_neib_pe
!    write(100+monolis_global_myrank(),*)"send_index", mat%com%send_index
!    write(100+monolis_global_myrank(),*)"send_item", mat%com%send_item

    !> solver
    call monolis_param_set_method(mat, monolis_iter_CG)
    call monolis_param_set_precond(mat, monolis_prec_DIAG)
    call monolis_param_set_maxiter(mat, 200000)
    call monolis_param_set_tol(mat, 1.0d-8)
    call monolis_param_set_show_time(mat, .true.)
    call monolis_param_set_show_iterlog(mat, .false.)
    call monolis_param_set_show_summary(mat, .true.)

    call monolis_solve(mat, b, a)

    ! write(*,*)a
    call monolis_get_solution(mat%MAT, mat%MAT%X)
    call monolis_finalize(mat)
  end subroutine monolis_solver_COCG_test

  subroutine get_coef_val(n_coef, coef, eid, val)
    implicit none
    integer(kint) :: eid, i, n_coef
    real(kdouble) :: val
    real(kdouble) :: coef(:)

    !val = cmplx(coef(1,eid), coef(2,eid))
    val = coef(eid)
  end subroutine get_coef_val

end program main
