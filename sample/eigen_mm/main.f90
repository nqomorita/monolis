program main
  use mod_monolis
  implicit none
  integer(4) :: nnode, nelem, i
  integer(4) :: n_eigs, maxit, loglevel, n_bc
  real(8) :: tol, t1, t2, t3, t4

  integer(4), allocatable :: elem(:,:), i_bc(:)
  real(8), allocatable :: coef(:)
  real(8), allocatable :: val(:), vec(:,:)

  type(monolis_structure) :: mat !> 疎行列変数
  character :: finA*100, finBC*100

  call monolis_global_initialize()

  t1 = monolis_get_time()

  call get_input_arg(finA, finBC, n_eigs, maxit, loglevel, tol)

  write(*,*) "filename matrix A: ", trim(finA)
  write(*,*) "filename matrix B: ", trim(finBC)

  call monolis_input_mtx(finA, nnode, nelem, elem, coef)

  call input_BC(finBC, n_bc, i_bc)

  call convert_matrix_by_bc(nnode, nelem, elem, coef, n_bc, i_bc)

  if(n_eigs > nnode) n_eigs = nnode

  call monolis_initialize(mat, "./") !> 疎行列変数の初期化

  call monolis_get_nonzero_pattern(mat, nnode, 2, 1, nelem, elem)

  do i = 1, nelem
    if(elem(1,i) == elem(2,i))then
      call monolis_add_scalar_to_sparse_matrix(mat, elem(1,i), elem(2,i), 1, 1, coef(i))
    else
      call monolis_add_scalar_to_sparse_matrix(mat, elem(1,i), elem(2,i), 1, 1, coef(i))
      call monolis_add_scalar_to_sparse_matrix(mat, elem(2,i), elem(1,i), 1, 1, coef(i))
    endif
  enddo

  allocate(val(n_eigs), source = 0.0d0)
  allocate(vec(nnode,n_eigs), source = 0.0d0)

  call monolis_param_set_method(mat, monolis_iter_CG)
  call monolis_param_set_precond(mat, monolis_prec_MUMPS)
  call monolis_param_set_maxiter(mat, 100000)
  call monolis_param_set_tol(mat, 1.0d-8)
  call monolis_param_set_is_scaling(mat, .false.)
  call monolis_param_set_is_reordering(mat, .false.)
  call monolis_param_set_is_init_x(mat, .true.)
  call monolis_param_set_is_debug(mat, .true.)
  call monolis_param_set_show_iterlog(mat, .false.)
  call monolis_param_set_show_time(mat, .false.)
  call monolis_param_set_show_summary(mat, .false.)
  t2 = monolis_get_time()

  call monolis_eigen_inverted_standard_lanczos(mat, n_eigs, tol, maxit, val, vec)

  t3 = monolis_get_time()

  call output(nnode, n_eigs, val, vec)

  call monolis_finalize(mat) !> 疎行列変数の解放

  t4 = monolis_get_time()

  call monolis_global_finalize()

  write(*,"(a,1pe12.5)") "* input   time:", t2-t1
  write(*,"(a,1pe12.5)") "* Lanczos time:", t3-t2
  write(*,"(a,1pe12.5)") "* output  time:", t4-t3
  write(*,"(a,1pe12.5)") "* total   time:", t4-t1

contains

  subroutine output(N, n_eigs, eigen_val, eigen_vec)
    implicit none
    integer(4) :: N, n_eigs, i, j
    real(8) :: eigen_val(:), eigen_vec(:,:)
    character :: cnum*4

    call system('if [ ! -d output ]; then (echo "* create output"; mkdir -p output); fi')

    write(*,"(a)")"* eigen value"
    write(*,"(1pe12.5)")eigen_val

    open(20, file = "output/eigen_value.txt", status = "replace")
    write(20,"(i4)")n_eigs
    do i = 1, n_eigs
      write(20, "(i4,1pe12.4)")i, eigen_val(i)
    enddo
    close(20)

    do i = 1, n_eigs
      write(cnum,"(i0)")i
      open(20, file = "output/eigen_vector."//trim(cnum)//".txt", status = "replace")
      write(20,"(i8)")N
      do j = 1, N
        write(20, "(1pe12.4)")eigen_vec(j,i)
      enddo
      close(20)
    enddo
  end subroutine output

  subroutine convert_matrix_by_bc(nnode, nelem, elem, coef, n_bc, i_bc)
    integer(4) :: nnode, nelem, n_bc, i_bc(:)
    integer(4) :: nnode_o, nelem_o, in, i1, i2
    integer(4), allocatable :: elem(:,:), elem_o(:,:), perm(:)
    real(8), allocatable :: coef(:), coef_o(:)
    logical, allocatable :: is_use(:)

    nnode_o = nnode
    nelem_o = nelem
    allocate(elem_o(2,nelem), source = 0)
    allocate(coef_o(nelem), source = 0.0d0)
    elem_o = elem
    coef_o = coef
    deallocate(elem)
    deallocate(coef)

    allocate(is_use(nnode), source = .true.)
    do i = 1, n_bc
      in = i_bc(i)
      is_use(in) = .false.
    enddo

    allocate(perm(nnode), source = 0)
    in = 1
    do i = 1, nnode
      if(is_use(i))then
        perm(i) = in
        in = in + 1
      endif
    enddo

    nnode = nnode - n_bc

    nelem = 0
    do i = 1, nelem_o
      i1 = elem_o(1,i)
      i2 = elem_o(2,i)
      if(is_use(i1) .and. is_use(i2))then
        nelem = nelem + 1
      endif
    enddo

    allocate(elem(2,nelem), source = 0)
    allocate(coef(nelem), source = 0.0d0)

    in = 1
    do i = 1, nelem_o
      i1 = elem_o(1,i)
      i2 = elem_o(2,i)
      if(is_use(i1) .and. is_use(i2))then
        elem(1,in) = perm(elem_o(1,i))
        elem(2,in) = perm(elem_o(2,i))
        coef(in) = coef_o(i)
        in = in + 1
      endif
    enddo
  end subroutine convert_matrix_by_bc

  subroutine get_input_arg(finA, finBC, n_eigs, maxit, loglevel, tol)
    implicit none
    integer(4) :: count
    integer(4) :: n_eigs, maxit, loglevel
    real(8) :: tol
    character :: finA*100, finBC*100

    count = iargc()
    if(count == 2)then
      call getarg(1, finA)
      call getarg(2, finBC)
    else
      stop "please enter input file name of matrix A and matrix BC"
    endif

    open(10, file = "setting.dat", status = "old")
      read(10,*)n_eigs
      read(10,*)maxit
      read(10,*)tol
      read(10,*)loglevel
    close(10)

    write(*,"(a,i12)")     "* n_eigs  :", n_eigs
    write(*,"(a,i12)")     "* maxiter :", maxit
    write(*,"(a,1pe12.5)") "* tol     :", tol
    write(*,"(a,i12)")     "* loglevel:", loglevel
  end subroutine get_input_arg

  subroutine input_BC(fin, n_bc, i_bc)
    implicit none
    integer(4) :: n_bc, i, n, ndof, id, idof
    integer(4), allocatable :: i_bc(:)
    real(8) :: val
    character :: fin*100

    open(20, file = trim(fin), status = "old")
      read(20,*)n_bc, ndof
      allocate(i_bc(n_bc), source = 0)
      do i = 1, n_bc
        read(20,*) id, idof, val
        i_bc(i) = ndof*(id-1) + idof
      enddo
    close(20)
  end subroutine input_BC

end program main
