program main
  use mod_monolis
  implicit none
  type(monolis_structure) :: mat !> 疎行列変数
  integer(kint) :: nnode, nelem, nbase_func, i
  real(kdouble) :: resid
  integer(kint), allocatable :: elem(:,:)
  real(kdouble), allocatable :: node(:,:), coef(:), B(:), X(:)
  character :: fname*128

  call monolis_std_debug_log_flag(.true.)

  call monolis_global_initialize()

  fname = "node.dat"
  call monolis_input_mesh_node(fname, nnode, node)

  fname = "elem.dat"
  call monolis_input_mesh_elem(fname, nelem, nbase_func, elem)

  open(20, file = "coef.dat", status = "old")
    read(20,*) i
    if(i /= nelem) stop "ncoef /= nelem"
    allocate(coef(nelem), source = 0.0d0)
    do i = 1, nelem
      read(20,*) coef(i)
    enddo
  close(20)

  call monolis_initialize(mat, "./") !> 疎行列変数の初期化

  call monolis_get_nonzero_pattern(mat, nnode, nbase_func, 1, nelem, elem)

  do i = 1, nelem
    if(elem(1,i) == elem(2,i))then
      call monolis_add_scalar_to_sparse_matrix(mat, elem(1,i), elem(2,i), 1, 1, coef(i))
    else
      call monolis_add_scalar_to_sparse_matrix(mat, elem(1,i), elem(2,i), 1, 1, coef(i))
      call monolis_add_scalar_to_sparse_matrix(mat, elem(2,i), elem(1,i), 1, 1, coef(i))
    endif
  enddo

  allocate(B(nnode), source = 0.0d0)
  allocate(X(nnode), source = 0.0d0)

  X = 1.0d0
  call monolis_matvec_product(mat, X, B)
  X = 0.0d0

  call monolis_solve(mat, B, X)

  resid = 0.0d0
  do i = 1, nnode
    resid = resid + (X(i) - 1.0d0)**2
  enddo
  write(*,*)"resid: ", dsqrt(resid)

  call monolis_finalize(mat) !> 疎行列変数の解放

  call monolis_global_finalize()

end program main
