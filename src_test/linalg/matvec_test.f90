!> 疎行列ベクトル積テストモジュール
module mod_monolis_matvec_test
  use mod_monolis
  use mod_monolis_utils
  use mod_monolis_matvec
  use mod_monolis_spmat_convert_dia
  use mod_monolis_spmat_convert_ell
  use mod_monolis_spmat_nonzero_pattern
  implicit none

contains

  subroutine monolis_matvec_test()
    implicit none

    if(monolis_mpi_get_global_comm_size() == 2) return

    call monolis_matvec_11_R_test()
    call monolis_matvec_11_C_test()

    call monolis_matvec_33_R_test()
    call monolis_matvec_33_C_test()

    call monolis_matvec_nn_R_test()
    call monolis_matvec_nn_C_test()

    call monolis_matvec_V_DIA_ELL_R_test()

    call monolis_std_global_log_string("monolis_matvec_product_R")
    call monolis_std_global_log_string("monolis_matvec_product_C")
    call monolis_std_global_log_string("monolis_matvec_product_main_R")
    call monolis_std_global_log_string("monolis_matvec_product_main_C")
    call monolis_std_global_log_string("monolis_matvec_V_R")
    call monolis_std_global_log_string("monolis_matvec_V_C")
    call monolis_std_global_log_string("monolis_matvec_DIA_nn_R")
    call monolis_std_global_log_string("monolis_matvec_ELL_nn_R")
    call monolis_std_global_log_string("monolis_residual_main_R")
    call monolis_std_global_log_string("monolis_residual_main_C")
  end subroutine monolis_matvec_test

  !> 可変ブロックの DIA / ELL 形式 SpMV を CSR 可変ブロック SpMV と突き合わせる
  subroutine monolis_matvec_V_DIA_ELL_R_test()
    implicit none
    type(monolis_structure) :: mat
    integer(kint) :: n_node, n_base, n_elem, elem(3,2)
    integer(kint) :: n_dof_list(4)
    integer(kint) :: N, NP, NNDOF, NPNDOF, k
    real(kdouble), allocatable :: X(:), Yref(:), Ydia(:), Yell(:)

    call monolis_std_global_log_string("monolis_matvec_DIA_V_R")
    call monolis_std_global_log_string("monolis_matvec_ELL_V_R")

    !# 4 節点・2 三角形要素、節点ごとに自由度 1,2,3,1 の可変ブロック
    n_node = 4
    n_base = 3
    n_elem = 2
    elem(1,1) = 1; elem(2,1) = 2; elem(3,1) = 3
    elem(1,2) = 2; elem(2,2) = 3; elem(3,2) = 4
    n_dof_list(1) = 1
    n_dof_list(2) = 2
    n_dof_list(3) = 3
    n_dof_list(4) = 1

    call monolis_initialize(mat)
    call monolis_get_nonzero_pattern_by_simple_mesh_V_R(mat, n_node, n_base, n_dof_list, n_elem, elem)

    N  = mat%MAT%N
    NP = mat%MAT%NP
    NNDOF  = mat%MAT%n_dof_index(N + 1)
    NPNDOF = mat%MAT%n_dof_index(NP + 1)

    !# 行列値・入力ベクトルを決定論的な値で埋める
    do k = 1, size(mat%MAT%R%A)
      mat%MAT%R%A(k) = dble(mod(k, 7) + 1)
    enddo

    call monolis_alloc_R_1d(X,    NPNDOF)
    call monolis_alloc_R_1d(Yref, NPNDOF)
    call monolis_alloc_R_1d(Ydia, NPNDOF)
    call monolis_alloc_R_1d(Yell, NPNDOF)
    do k = 1, NPNDOF
      X(k) = dble(mod(k, 5) + 1)
    enddo

    !# 参照解（CSR 可変ブロック SpMV）
    call monolis_matvec_V_R(N, mat%MAT%n_dof_list, mat%MAT%n_dof_index, mat%MAT%n_dof_index2, &
      mat%MAT%CSR%index, mat%MAT%CSR%item, mat%MAT%R%A, X, Yref)

    !# DIA 可変ブロック SpMV
    call monolis_convert_CSR_to_DIA_V_R(mat%MAT)
    call monolis_matvec_DIA_V_R(N, NP, mat%MAT%DIA%Ndiag, mat%MAT%DIA%offset, mat%MAT%DIA%Vptr, &
      mat%MAT%R%Adia, mat%MAT%n_dof_list, mat%MAT%n_dof_index, X, Ydia)
    call monolis_test_check_eq_R("monolis_matvec_DIA_V_R_test", Ydia(1:NNDOF), Yref(1:NNDOF))
    call monolis_dealloc_DIA_R(mat%MAT)

    !# ELL 可変ブロック SpMV
    call monolis_convert_CSR_to_ELL_V_R(mat%MAT)
    call monolis_matvec_ELL_V_R(N, NP, mat%MAT%ELL%Nmaxcol, mat%MAT%ELL%col, mat%MAT%ELL%Vptr, &
      mat%MAT%R%Aell, mat%MAT%n_dof_list, mat%MAT%n_dof_index, X, Yell)
    call monolis_test_check_eq_R("monolis_matvec_ELL_V_R_test", Yell(1:NNDOF), Yref(1:NNDOF))
    call monolis_dealloc_ELL_R(mat%MAT)

    call monolis_dealloc_R_1d(X)
    call monolis_dealloc_R_1d(Yref)
    call monolis_dealloc_R_1d(Ydia)
    call monolis_dealloc_R_1d(Yell)
    call monolis_finalize(mat)
  end subroutine monolis_matvec_V_DIA_ELL_R_test

  subroutine monolis_matvec_11_R_test()
    implicit none
    type(monolis_structure) :: mat
    type(monolis_com) :: com
    integer(kint) :: n_node, n_elem, elem(2,4)
    real(kdouble) :: a(5), b(5), b_th(5), mat_dense(5,5)

    call monolis_std_global_log_string("monolis_matvec_11_R")

    call monolis_initialize(mat)
    call monolis_com_initialize_by_self(com)

    n_node = 5

    n_elem = 4

    elem(1,1) = 1; elem(2,1) = 2;
    elem(1,2) = 2; elem(2,2) = 3;
    elem(1,3) = 3; elem(2,3) = 4;
    elem(1,4) = 4; elem(2,4) = 5;

    call monolis_get_nonzero_pattern_by_simple_mesh_R(mat, n_node, 2, 1, n_elem, elem)

    call monolis_add_scalar_to_sparse_matrix_R(mat, 1, 1, 1, 1, 2.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 1, 2, 1, 1, 1.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 2, 1, 1, 1, 1.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 2, 2, 1, 1, 2.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 2, 3, 1, 1, 3.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 3, 2, 1, 1, 1.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 3, 3, 1, 1, 2.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 3, 4, 1, 1, 4.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 4, 3, 1, 1, 1.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 4, 4, 1, 1, 2.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 4, 5, 1, 1, 5.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 5, 4, 1, 1, 1.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 5, 5, 1, 1, 2.0d0)

    mat_dense(1, 1) = 2.0d0
    mat_dense(1, 2) = 1.0d0
    mat_dense(2, 1) = 1.0d0
    mat_dense(2, 2) = 2.0d0
    mat_dense(2, 3) = 3.0d0
    mat_dense(3, 2) = 1.0d0
    mat_dense(3, 3) = 2.0d0
    mat_dense(3, 4) = 4.0d0
    mat_dense(4, 3) = 1.0d0
    mat_dense(4, 4) = 2.0d0
    mat_dense(4, 5) = 5.0d0
    mat_dense(5, 4) = 1.0d0
    mat_dense(5, 5) = 2.0d0

    a(1) = 1.0d0
    a(2) = 1.0d0
    a(3) = 1.0d0
    a(4) = 1.0d0
    a(5) = 1.0d0

    call monolis_matvec_product_R(mat, com, a, b)

    b_th = matmul(mat_dense, a)

    call monolis_test_check_eq_R("monolis_matvec_11_R_test", b, b_th)

    call monolis_finalize(mat)
  end subroutine monolis_matvec_11_R_test

  subroutine monolis_matvec_11_C_test()
    implicit none
    type(monolis_structure) :: mat
    type(monolis_com) :: com
    integer(kint) :: n_node, n_elem, elem(2,4)
    complex(kdouble) :: a(5), b(5), b_th(5), mat_dense(5,5)

    call monolis_std_global_log_string("monolis_matvec_11_C")

    call monolis_initialize(mat)
    call monolis_com_initialize_by_self(com)

    n_node = 5

    n_elem = 4

    elem(1,1) = 1; elem(2,1) = 2;
    elem(1,2) = 2; elem(2,2) = 3;
    elem(1,3) = 3; elem(2,3) = 4;
    elem(1,4) = 4; elem(2,4) = 5;

    call monolis_get_nonzero_pattern_by_simple_mesh_C(mat, n_node, 2, 1, n_elem, elem)

    call monolis_add_scalar_to_sparse_matrix_C(mat, 1, 1, 1, 1, (2.0d0, 2.0d0))
    call monolis_add_scalar_to_sparse_matrix_C(mat, 1, 2, 1, 1, (1.0d0, 1.0d0))
    call monolis_add_scalar_to_sparse_matrix_C(mat, 2, 1, 1, 1, (1.0d0, 1.0d0))
    call monolis_add_scalar_to_sparse_matrix_C(mat, 2, 2, 1, 1, (2.0d0, 2.0d0))
    call monolis_add_scalar_to_sparse_matrix_C(mat, 2, 3, 1, 1, (3.0d0, 3.0d0))
    call monolis_add_scalar_to_sparse_matrix_C(mat, 3, 2, 1, 1, (1.0d0, 1.0d0))
    call monolis_add_scalar_to_sparse_matrix_C(mat, 3, 3, 1, 1, (2.0d0, 2.0d0))
    call monolis_add_scalar_to_sparse_matrix_C(mat, 3, 4, 1, 1, (4.0d0, 4.0d0))
    call monolis_add_scalar_to_sparse_matrix_C(mat, 4, 3, 1, 1, (1.0d0, 1.0d0))
    call monolis_add_scalar_to_sparse_matrix_C(mat, 4, 4, 1, 1, (2.0d0, 2.0d0))
    call monolis_add_scalar_to_sparse_matrix_C(mat, 4, 5, 1, 1, (5.0d0, 5.0d0))
    call monolis_add_scalar_to_sparse_matrix_C(mat, 5, 4, 1, 1, (1.0d0, 1.0d0))
    call monolis_add_scalar_to_sparse_matrix_C(mat, 5, 5, 1, 1, (2.0d0, 2.0d0))

    mat_dense = (0.0d0, 0.0d0)
    mat_dense(1, 1) = (2.0d0, 2.0d0)
    mat_dense(1, 2) = (1.0d0, 1.0d0)
    mat_dense(2, 1) = (1.0d0, 1.0d0)
    mat_dense(2, 2) = (2.0d0, 2.0d0)
    mat_dense(2, 3) = (3.0d0, 3.0d0)
    mat_dense(3, 2) = (1.0d0, 1.0d0)
    mat_dense(3, 3) = (2.0d0, 2.0d0)
    mat_dense(3, 4) = (4.0d0, 4.0d0)
    mat_dense(4, 3) = (1.0d0, 1.0d0)
    mat_dense(4, 4) = (2.0d0, 2.0d0)
    mat_dense(4, 5) = (5.0d0, 5.0d0)
    mat_dense(5, 4) = (1.0d0, 1.0d0)
    mat_dense(5, 5) = (2.0d0, 2.0d0)

    a(1) = (1.0d0, 1.0d0)
    a(2) = (1.0d0, 1.0d0)
    a(3) = (1.0d0, 1.0d0)
    a(4) = (1.0d0, 1.0d0)
    a(5) = (1.0d0, 1.0d0)

    call monolis_matvec_product_C(mat, com, a, b)

    b_th = matmul(mat_dense, a)

    call monolis_test_check_eq_C("monolis_matvec_11_C_test", b, b_th)

    call monolis_finalize(mat)
  end subroutine monolis_matvec_11_C_test

  subroutine monolis_matvec_33_R_test()
    implicit none
    type(monolis_structure) :: mat
    type(monolis_com) :: com
    integer(kint) :: n_node, n_elem, elem(2,3)
    integer(kint) :: i1, i2, j1, j2
    real(kdouble) :: val
    real(kdouble) :: a(12), b(12), b_th(12), mat_dense(12,12)
    logical :: is_find

    call monolis_std_global_log_string("monolis_matvec_33_R")

    call monolis_initialize(mat)
    call monolis_com_initialize_by_self(com)

    n_node = 4

    n_elem = 3

    elem(1,1) = 1; elem(2,1) = 2;
    elem(1,2) = 2; elem(2,2) = 3;
    elem(1,3) = 3; elem(2,3) = 4;

    call monolis_get_nonzero_pattern_by_simple_mesh_R(mat, n_node, 2, 3, n_elem, elem)

    do i1 = 1, 3
      do i2 = 1, 3
      do j2 = 1, 3
        call random_number(val)
        val = val + 1.0d0
        call monolis_add_scalar_to_sparse_matrix_R(mat, elem(1,i1), elem(2,i1), i2, j2, val)
        call monolis_add_scalar_to_sparse_matrix_R(mat, elem(2,i1), elem(1,i1), j2, i2, val)
      enddo
      enddo
    enddo

    a = 1.0d0

    call monolis_matvec_product_R(mat, com, a, b)

    mat_dense = 0.0d0

    do i1 = 1, 4
    do j1 = 1, 4
      do i2 = 1, 3
      do j2 = 1, 3
        call monolis_get_scalar_from_sparse_matrix_R(mat, i1, j1, i2, j2, val, is_find)
        mat_dense(3*i1-3+i2, 3*j1-3+j2) = val
      enddo
      enddo
    enddo
    enddo

    b_th = matmul(mat_dense, a)

    call monolis_test_check_eq_R("monolis_matvec_33_R_test", b, b_th)

    call monolis_finalize(mat)
  end subroutine monolis_matvec_33_R_test

  subroutine monolis_matvec_33_C_test()
    implicit none
    type(monolis_structure) :: mat
    type(monolis_com) :: com
    integer(kint) :: n_node, n_elem, elem(2,3)
    integer(kint) :: i1, i2, j1, j2
    real(kdouble) :: v1, v2
    complex(kdouble) :: val
    complex(kdouble) :: a(12), b(12), b_th(12), mat_dense(12,12)
    logical :: is_find

    call monolis_std_global_log_string("monolis_matvec_33_C")

    call monolis_initialize(mat)
    call monolis_com_initialize_by_self(com)

    n_node = 4

    n_elem = 3

    elem(1,1) = 1; elem(2,1) = 2;
    elem(1,2) = 2; elem(2,2) = 3;
    elem(1,3) = 3; elem(2,3) = 4;

    call monolis_get_nonzero_pattern_by_simple_mesh_C(mat, n_node, 2, 3, n_elem, elem)

    do i1 = 1, 3
      do i2 = 1, 3
      do j2 = 1, 3
        call random_number(v1)
        call random_number(v2)
        val = cmplx(v1, v2) + (1.0d0, 1.0d0)
        call monolis_add_scalar_to_sparse_matrix_C(mat, elem(1,i1), elem(2,i1), i2, j2, val)
        call monolis_add_scalar_to_sparse_matrix_C(mat, elem(2,i1), elem(1,i1), j2, i2, val)
      enddo
      enddo
    enddo

    a = (1.0d0, 1.0d0)

    call monolis_matvec_product_C(mat, com, a, b)

    mat_dense = (0.0d0, 0.0d0)

    do i1 = 1, 4
    do j1 = 1, 4
      do i2 = 1, 3
      do j2 = 1, 3
        call monolis_get_scalar_from_sparse_matrix_C(mat, i1, j1, i2, j2, val, is_find)
        mat_dense(3*i1-3+i2, 3*j1-3+j2) = val
      enddo
      enddo
    enddo
    enddo

    b_th = matmul(mat_dense, a)

    call monolis_test_check_eq_C("monolis_matvec_33_C_test", b, b_th)

    call monolis_finalize(mat)
  end subroutine monolis_matvec_33_C_test

  subroutine monolis_matvec_nn_R_test()
    implicit none
    type(monolis_structure) :: mat
    type(monolis_com) :: com
    integer(kint) :: n_node, n_elem, elem(2,3)
    integer(kint) :: i1, i2, j1, j2
    real(kdouble) :: val
    real(kdouble) :: a(8), b(8), b_th(8), mat_dense(8,8)
    logical :: is_find

    call monolis_std_global_log_string("monolis_matvec_nn_R")

    call monolis_initialize(mat)
    call monolis_com_initialize_by_self(com)

    n_node = 4

    n_elem = 3

    elem(1,1) = 1; elem(2,1) = 2;
    elem(1,2) = 2; elem(2,2) = 3;
    elem(1,3) = 3; elem(2,3) = 4;

    call monolis_get_nonzero_pattern_by_simple_mesh_R(mat, n_node, 2, 2, n_elem, elem)

    do i1 = 1, 3
      do i2 = 1, 2
      do j2 = 1, 2
        call random_number(val)
        val = val + 1.0d0
        call monolis_add_scalar_to_sparse_matrix_R(mat, elem(1,i1), elem(2,i1), i2, j2, val)
        call monolis_add_scalar_to_sparse_matrix_R(mat, elem(2,i1), elem(1,i1), j2, i2, val)
      enddo
      enddo
    enddo

    a = 1.0d0

    call monolis_matvec_product_R(mat, com, a, b)

    mat_dense = 0.0d0

    do i1 = 1, 4
    do j1 = 1, 4
      do i2 = 1, 2
      do j2 = 1, 2
        call monolis_get_scalar_from_sparse_matrix_R(mat, i1, j1, i2, j2, val, is_find)
        mat_dense(2*i1-2+i2, 2*j1-2+j2) = val
      enddo
      enddo
    enddo
    enddo

    b_th = matmul(mat_dense, a)

    call monolis_test_check_eq_R("monolis_matvec_nn_R_test", b, b_th)

    call monolis_finalize(mat)
  end subroutine monolis_matvec_nn_R_test

  subroutine monolis_matvec_nn_C_test()
    implicit none
    type(monolis_structure) :: mat
    type(monolis_com) :: com
    integer(kint) :: n_node, n_elem, elem(2,3)
    integer(kint) :: i1, i2, j1, j2
    real(kdouble) :: v1, v2
    complex(kdouble) :: val
    complex(kdouble) :: a(8), b(8), b_th(8), mat_dense(8,8)
    logical :: is_find

    call monolis_std_global_log_string("monolis_matvec_nn_C")

    call monolis_initialize(mat)
    call monolis_com_initialize_by_self(com)

    n_node = 4

    n_elem = 3

    elem(1,1) = 1; elem(2,1) = 2;
    elem(1,2) = 2; elem(2,2) = 3;
    elem(1,3) = 3; elem(2,3) = 4;

    call monolis_get_nonzero_pattern_by_simple_mesh_C(mat, n_node, 2, 2, n_elem, elem)

    do i1 = 1, 3
      do i2 = 1, 2
      do j2 = 1, 2
        call random_number(v1)
        call random_number(v2)
        val = cmplx(v1, v2) + (1.0d0, 1.0d0)
        call monolis_add_scalar_to_sparse_matrix_C(mat, elem(1,i1), elem(2,i1), i2, j2, val)
        call monolis_add_scalar_to_sparse_matrix_C(mat, elem(2,i1), elem(1,i1), j2, i2, val)
      enddo
      enddo
    enddo

    a = (1.0d0, 1.0d0)

    call monolis_matvec_product_C(mat, com, a, b)

    mat_dense = (0.0d0, 0.0d0)

    do i1 = 1, 4
    do j1 = 1, 4
      do i2 = 1, 2
      do j2 = 1, 2
        call monolis_get_scalar_from_sparse_matrix_C(mat, i1, j1, i2, j2, val, is_find)
        mat_dense(2*i1-2+i2, 2*j1-2+j2) = val
      enddo
      enddo
    enddo
    enddo

    b_th = matmul(mat_dense, a)

    call monolis_test_check_eq_C("monolis_matvec_nn_C_test", b, b_th)

    call monolis_finalize(mat)
  end subroutine monolis_matvec_nn_C_test
end module mod_monolis_matvec_test
