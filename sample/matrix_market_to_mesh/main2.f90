program main
  use mod_monolis
  implicit none
  type(monolis_structure) :: mat !> 疎行列変数
  integer(4) :: nnode, nelem, nbase_func, i
  integer(4), allocatable :: elem(:,:)
  real(kdouble), allocatable :: node(:,:), coef(:)
  character :: fname*128

  call monolis_set_debug(.true.)

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

  !call monolis_solve(mat)

  call monolis_finalize(mat) !> 疎行列変数の解放

  call monolis_global_finalize()

end program main
