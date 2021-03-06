program main
  use mod_monolis
  implicit none
  integer(4) :: nnode, nelem, elem(2,5)

  type(monolis_structure) :: mat !> 疎行列変数

  call monolis_global_initialize()

  call monolis_initialize(mat) !> 疎行列変数の初期化

  nnode = 6
  nelem = 5
  elem(1,1) = 1; elem(2,1) = 2;
  elem(1,2) = 2; elem(2,2) = 3;
  elem(1,3) = 3; elem(2,3) = 4;
  elem(1,4) = 4; elem(2,4) = 5;
  elem(1,5) = 5; elem(2,5) = 6;

  call monolis_get_nonzero_pattern(mat, nnode, 2, 1, nelem, elem)

  !> 16  2  0  0  0
  !>  2 12  3  0  0
  !>  0  3  1  2  0
  !>  0  0  2 20  5
  !>  0  0  0  5 30
  !> call monolis_add_scalar_to_sparse_matrix(mat, 1, 1, 1, 1, 16.0d0)
  !> call monolis_add_scalar_to_sparse_matrix(mat, 1, 2, 1, 1,  2.0d0)
  !> call monolis_add_scalar_to_sparse_matrix(mat, 2, 1, 1, 1,  2.0d0)
  !> call monolis_add_scalar_to_sparse_matrix(mat, 2, 2, 1, 1, 12.0d0)
  !> call monolis_add_scalar_to_sparse_matrix(mat, 2, 3, 1, 1,  3.0d0)
  !> call monolis_add_scalar_to_sparse_matrix(mat, 3, 2, 1, 1,  3.0d0)
  !> call monolis_add_scalar_to_sparse_matrix(mat, 3, 3, 1, 1,  1.0d0)
  !> call monolis_add_scalar_to_sparse_matrix(mat, 3, 4, 1, 1,  2.0d0)
  !> call monolis_add_scalar_to_sparse_matrix(mat, 4, 3, 1, 1,  2.0d0)
  !> call monolis_add_scalar_to_sparse_matrix(mat, 4, 4, 1, 1, 20.0d0)
  !> call monolis_add_scalar_to_sparse_matrix(mat, 4, 5, 1, 1,  5.0d0)
  !> call monolis_add_scalar_to_sparse_matrix(mat, 5, 4, 1, 1,  5.0d0)
  !> call monolis_add_scalar_to_sparse_matrix(mat, 5, 5, 1, 1, 30.0d0)

  !>  2  1  0  0  0
  !>  1  2  1  0  0
  !>  0  1  2  2  0
  !>  0  0  1  2  1
  !>  0  0  0  1  2
  call monolis_add_scalar_to_sparse_matrix(mat, 1, 1, 1, 1,  2.0d0)
  call monolis_add_scalar_to_sparse_matrix(mat, 1, 2, 1, 1,  1.0d0)

  call monolis_add_scalar_to_sparse_matrix(mat, 2, 1, 1, 1,  1.0d0)
  call monolis_add_scalar_to_sparse_matrix(mat, 2, 2, 1, 1,  2.0d0)
  call monolis_add_scalar_to_sparse_matrix(mat, 2, 3, 1, 1,  1.0d0)

  call monolis_add_scalar_to_sparse_matrix(mat, 3, 2, 1, 1,  1.0d0)
  call monolis_add_scalar_to_sparse_matrix(mat, 3, 3, 1, 1,  2.0d0)
  call monolis_add_scalar_to_sparse_matrix(mat, 3, 4, 1, 1,  1.0d0)

  call monolis_add_scalar_to_sparse_matrix(mat, 4, 3, 1, 1,  1.0d0)
  call monolis_add_scalar_to_sparse_matrix(mat, 4, 4, 1, 1,  2.0d0)
  call monolis_add_scalar_to_sparse_matrix(mat, 4, 5, 1, 1,  1.0d0)

  call monolis_add_scalar_to_sparse_matrix(mat, 5, 4, 1, 1,  1.0d0)
  call monolis_add_scalar_to_sparse_matrix(mat, 5, 5, 1, 1,  2.0d0)
  call monolis_add_scalar_to_sparse_matrix(mat, 5, 6, 1, 1,  1.0d0)

  call monolis_add_scalar_to_sparse_matrix(mat, 6, 5, 1, 1,  1.0d0)
  call monolis_add_scalar_to_sparse_matrix(mat, 6, 6, 1, 1,  2.0d0)

  !> eigen_value
  !> 3.73205E+00
  !> 3.00000E+00
  !> 2.00000E+00
  !> 1.00000E+00
  !> 2.67949E-01

  !> e_mode
  !> 2.88675E-01 5.00000E-01 5.77350E-01 5.00000E-01 2.88675E-01
  !> 5.00000E-01 5.00000E-01-2.42335E-16-5.00000E-01-5.00000E-01
  !> 5.77350E-01-6.48664E-17-5.77350E-01 1.61183E-16 5.77350E-01
  !> 5.00000E-01-5.00000E-01-2.82388E-16 5.00000E-01-5.00000E-01
  !> 2.88675E-01-5.00000E-01 5.77350E-01-5.00000E-01 2.88675E-01

  !> 3.80194   3.24698   2.44504   1.55496   0.75302   0.198063
  !>
  !> 0.231921  0.417907  0.521121  0.521121  0.417907  0.231921
  !> 0.417907  0.521121  0.231921 -0.231921 -0.521121 -0.417906
  !> 0.521121  0.231921 -0.417907 -0.417907  0.231921  0.521121
  !> 0.521121 -0.231921 -0.417907  0.417907  0.231921 -0.521121
  !> 0.417906 -0.521121  0.231921  0.231921 -0.521121  0.417907
  !> 0.231921 -0.417907  0.521121 -0.521121  0.417907 -0.231921

  !mat%PRM%is_debug = .true.
  mat%PRM%show_iterlog = .true.
  mat%PRM%show_time = .false.
  mat%PRM%show_summary = .false.

  !mat%PRM%precond = monolis_prec_NONE
  !mat%PRM%precond = monolis_prec_JACOBI
  !mat%PRM%precond = monolis_prec_ILU
  mat%PRM%precond = monolis_prec_SOR
  !call monolis_eigen_inverted_lobpcg(mat, 1, 1.0d-6)
  call monolis_eigen_inverted_lobpcg(mat, 2, 1.0d-6)
  !call monolis_eigen_inverted_lobpcg(mat, 3, 1.0d-6)
  !call monolis_eigen_inverted_lobpcg(mat, 4, 1.0d-6)
  !call monolis_eigen_inverted_lobpcg(mat, 5, 1.0d-6)

  call monolis_finalize(mat) !> 疎行列変数の解放

  call monolis_global_finalize()
end program main
