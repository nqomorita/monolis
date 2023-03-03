!> 疎行列ベクトル積テストモジュール
module mod_monolis_matvec_test
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_matvec
  implicit none

contains

  subroutine monolis_matvec_test()
    implicit none

    call monolis_matvec_11_R_test()
    call monolis_matvec_11_C_test()
  end subroutine monolis_matvec_test

  subroutine monolis_matvec_11_R_test()
    implicit none
    type(monolis_structure) :: mat
    integer(kint) :: nnode, nelem, elem(2,4)
    real(kdouble) :: a(5), b(5), b_th(5), mat_dense(5,5)

    call monolis_std_log_string("monolis_set_scalar_to_sparse_matrix_R_test")

    call monolis_initialize(mat, "./")

    nnode = 5
    nelem = 4
    elem(1,1) = 1; elem(2,1) = 2;
    elem(1,2) = 2; elem(2,2) = 3;
    elem(1,3) = 3; elem(2,3) = 4;
    elem(1,4) = 4; elem(2,4) = 5;

    call monolis_get_nonzero_pattern(mat, nnode, 2, 1, nelem, elem)

    call monolis_add_scalar_to_sparse_matrix(mat, 1, 1, 1, 1, (2.0d0, 2.0d0))
    call monolis_add_scalar_to_sparse_matrix(mat, 1, 2, 1, 1, (1.0d0, 1.0d0))
    call monolis_add_scalar_to_sparse_matrix(mat, 2, 1, 1, 1, (1.0d0, 1.0d0))
    call monolis_add_scalar_to_sparse_matrix(mat, 2, 2, 1, 1, (2.0d0, 2.0d0))
    call monolis_add_scalar_to_sparse_matrix(mat, 2, 3, 1, 1, (1.0d0, 1.0d0))
    call monolis_add_scalar_to_sparse_matrix(mat, 3, 2, 1, 1, (1.0d0, 1.0d0))
    call monolis_add_scalar_to_sparse_matrix(mat, 3, 3, 1, 1, (2.0d0, 2.0d0))
    call monolis_add_scalar_to_sparse_matrix(mat, 3, 4, 1, 1, (1.0d0, 1.0d0))
    call monolis_add_scalar_to_sparse_matrix(mat, 4, 3, 1, 1, (1.0d0, 1.0d0))
    call monolis_add_scalar_to_sparse_matrix(mat, 4, 4, 1, 1, (2.0d0, 2.0d0))
    call monolis_add_scalar_to_sparse_matrix(mat, 4, 5, 1, 1, (1.0d0, 1.0d0))
    call monolis_add_scalar_to_sparse_matrix(mat, 5, 4, 1, 1, (1.0d0, 1.0d0))
    call monolis_add_scalar_to_sparse_matrix(mat, 5, 5, 1, 1, (2.0d0, 2.0d0))

    mat_dense = (0.0d0, 0.0d0)
    mat_dense(1, 1) = (2.0d0, 2.0d0)
    mat_dense(1, 2) = (1.0d0, 1.0d0)
    mat_dense(2, 1) = (1.0d0, 1.0d0)
    mat_dense(2, 2) = (2.0d0, 2.0d0)
    mat_dense(2, 3) = (1.0d0, 1.0d0)
    mat_dense(3, 2) = (1.0d0, 1.0d0)
    mat_dense(3, 3) = (2.0d0, 2.0d0)
    mat_dense(3, 4) = (1.0d0, 1.0d0)
    mat_dense(4, 3) = (1.0d0, 1.0d0)
    mat_dense(4, 4) = (2.0d0, 2.0d0)
    mat_dense(4, 5) = (1.0d0, 1.0d0)
    mat_dense(5, 4) = (1.0d0, 1.0d0)
    mat_dense(5, 5) = (2.0d0, 2.0d0)

    a(1) = (1.0d0, 1.0d0)
    a(2) = (1.0d0, 1.0d0)
    a(3) = (1.0d0, 1.0d0)
    a(4) = (1.0d0, 1.0d0)
    a(5) = (1.0d0, 1.0d0)

    call monolis_matvec_product(mat, a, b)

    b_th = matmul(mat_dense, a)

    call check_vector_complex_eq("monolis_matvec_product_test", b, b_th)

    call monolis_finalize(mat)

    call monolis_matvec_11_R(monoCOM, monoMAT, X, Y)
  end subroutine monolis_matvec_11_R_test

  subroutine monolis_matvec_11_C_test()
    implicit none
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    real(kdouble) :: X(4)
    real(kdouble) :: Y(4)
    integer(kint) :: NDOF

    integer(kint) :: nnode, nelem, elem(2,4)
    real(kdouble) :: tspmv, tcomm
    complex(kdouble) :: a(5), b(5), b_th(5), mat_dense(5,5)
    type(monolis_structure) :: mat

    call monolis_initialize(mat, "./")

    nnode = 5
    nelem = 4
    elem(1,1) = 1; elem(2,1) = 2;
    elem(1,2) = 2; elem(2,2) = 3;
    elem(1,3) = 3; elem(2,3) = 4;
    elem(1,4) = 4; elem(2,4) = 5;

    call monolis_get_nonzero_pattern(mat, nnode, 2, 1, nelem, elem)

    call monolis_add_scalar_to_sparse_matrix(mat, 1, 1, 1, 1, (2.0d0, 2.0d0))
    call monolis_add_scalar_to_sparse_matrix(mat, 1, 2, 1, 1, (1.0d0, 1.0d0))
    call monolis_add_scalar_to_sparse_matrix(mat, 2, 1, 1, 1, (1.0d0, 1.0d0))
    call monolis_add_scalar_to_sparse_matrix(mat, 2, 2, 1, 1, (2.0d0, 2.0d0))
    call monolis_add_scalar_to_sparse_matrix(mat, 2, 3, 1, 1, (1.0d0, 1.0d0))
    call monolis_add_scalar_to_sparse_matrix(mat, 3, 2, 1, 1, (1.0d0, 1.0d0))
    call monolis_add_scalar_to_sparse_matrix(mat, 3, 3, 1, 1, (2.0d0, 2.0d0))
    call monolis_add_scalar_to_sparse_matrix(mat, 3, 4, 1, 1, (1.0d0, 1.0d0))
    call monolis_add_scalar_to_sparse_matrix(mat, 4, 3, 1, 1, (1.0d0, 1.0d0))
    call monolis_add_scalar_to_sparse_matrix(mat, 4, 4, 1, 1, (2.0d0, 2.0d0))
    call monolis_add_scalar_to_sparse_matrix(mat, 4, 5, 1, 1, (1.0d0, 1.0d0))
    call monolis_add_scalar_to_sparse_matrix(mat, 5, 4, 1, 1, (1.0d0, 1.0d0))
    call monolis_add_scalar_to_sparse_matrix(mat, 5, 5, 1, 1, (2.0d0, 2.0d0))

    mat_dense = (0.0d0, 0.0d0)
    mat_dense(1, 1) = (2.0d0, 2.0d0)
    mat_dense(1, 2) = (1.0d0, 1.0d0)
    mat_dense(2, 1) = (1.0d0, 1.0d0)
    mat_dense(2, 2) = (2.0d0, 2.0d0)
    mat_dense(2, 3) = (1.0d0, 1.0d0)
    mat_dense(3, 2) = (1.0d0, 1.0d0)
    mat_dense(3, 3) = (2.0d0, 2.0d0)
    mat_dense(3, 4) = (1.0d0, 1.0d0)
    mat_dense(4, 3) = (1.0d0, 1.0d0)
    mat_dense(4, 4) = (2.0d0, 2.0d0)
    mat_dense(4, 5) = (1.0d0, 1.0d0)
    mat_dense(5, 4) = (1.0d0, 1.0d0)
    mat_dense(5, 5) = (2.0d0, 2.0d0)

    a(1) = (1.0d0, 1.0d0)
    a(2) = (1.0d0, 1.0d0)
    a(3) = (1.0d0, 1.0d0)
    a(4) = (1.0d0, 1.0d0)
    a(5) = (1.0d0, 1.0d0)

    call monolis_matvec_product(mat, a, b)

    b_th = matmul(mat_dense, a)

    call check_vector_complex_eq("monolis_matvec_product_test", b, b_th)

    call monolis_finalize(mat)

    call monolis_matvec_11_C(monoCOM, monoMAT, X, Y)
  end subroutine monolis_matvec_11_C_test

    !call monolis_matvec_11_R(monoCOM, monoMAT, X, Y)
    !call monolis_matvec_33_R(monoCOM, monoMAT, X, Y)
    !call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm)

end module mod_monolis_matvec_test
