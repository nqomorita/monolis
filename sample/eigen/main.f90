program main
  use mod_monolis
  implicit none
  integer(4) :: nnode, nelem, elem(2,4)

  type(monolis_structure) :: mat !> 疎行列変数

  call monolis_global_initialize()

  call monolis_initialize(mat) !> 疎行列変数の初期化

  nnode = 5
  nelem = 4
  elem(1,1) = 1; elem(2,1) = 2;
  elem(1,2) = 2; elem(2,2) = 3;
  elem(1,3) = 3; elem(2,3) = 4;
  elem(1,4) = 4; elem(2,4) = 5;

  call monolis_get_nonzero_pattern(mat, nnode, 2, 1, nelem, elem)

  !> 16  2  0  0  0
  !>  2 12  3  0  0
  !>  0  3  1  2  0
  !>  0  0  2 20  5
  !>  0  0  0  5 30
  call monolis_add_scalar_to_sparse_matrix(mat, 1, 1, 1, 1, 16.0d0)
  call monolis_add_scalar_to_sparse_matrix(mat, 1, 2, 1, 1,  2.0d0)
  call monolis_add_scalar_to_sparse_matrix(mat, 2, 1, 1, 1,  2.0d0)
  call monolis_add_scalar_to_sparse_matrix(mat, 2, 2, 1, 1, 12.0d0)
  call monolis_add_scalar_to_sparse_matrix(mat, 2, 3, 1, 1,  3.0d0)
  call monolis_add_scalar_to_sparse_matrix(mat, 3, 2, 1, 1,  3.0d0)
  call monolis_add_scalar_to_sparse_matrix(mat, 3, 3, 1, 1,  1.0d0)
  call monolis_add_scalar_to_sparse_matrix(mat, 3, 4, 1, 1,  2.0d0)
  call monolis_add_scalar_to_sparse_matrix(mat, 4, 3, 1, 1,  2.0d0)
  call monolis_add_scalar_to_sparse_matrix(mat, 4, 4, 1, 1, 20.0d0)
  call monolis_add_scalar_to_sparse_matrix(mat, 4, 5, 1, 1,  5.0d0)
  call monolis_add_scalar_to_sparse_matrix(mat, 5, 4, 1, 1,  5.0d0)
  call monolis_add_scalar_to_sparse_matrix(mat, 5, 5, 1, 1, 30.0d0)

  !mat%PRM%is_debug = .true.
  mat%PRM%is_reordering = .false.
  mat%PRM%show_iterlog = .false.
  mat%PRM%show_time = .false.
  call monolis_eigen_standard_lanczos(mat)

  call monolis_finalize(mat) !> 疎行列変数の解放

  call monolis_global_finalize()
end program main
