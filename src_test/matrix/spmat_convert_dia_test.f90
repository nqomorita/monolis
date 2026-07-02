!> 疎行列 DIA 変換関数群テスト
module mod_monolis_spmat_convert_dia_test
  use mod_monolis
  use mod_monolis_spmat_convert_dia
  use mod_monolis_spmat_nonzero_pattern

  implicit none

contains

  subroutine monolis_spmat_convert_dia_test()
    implicit none
    type(monolis_structure) :: mat
    integer(kint) :: n_node, n_elem, elem(2,4)
    integer(kint) :: N, Ndiag
    real(kdouble) :: v

    call monolis_std_global_log_string("monolis_convert_CSR_to_DIA_R")
    call monolis_std_global_log_string("monolis_dealloc_DIA_R")

    call monolis_initialize(mat)

    n_node = 5
    n_elem = 4

    if(monolis_mpi_get_global_comm_size() == 2) return

    elem(1,1) = 1; elem(2,1) = 2;
    elem(1,2) = 2; elem(2,2) = 3;
    elem(1,3) = 3; elem(2,3) = 4;
    elem(1,4) = 4; elem(2,4) = 5;

    call monolis_get_nonzero_pattern_by_simple_mesh_R(mat, n_node, 2, 1, n_elem, elem)

    !# 三重対角行列を構築（オフセット -1, 0, +1）
    call monolis_add_scalar_to_sparse_matrix_R(mat, 1, 1, 1, 1, 2.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 1, 2, 1, 1, 1.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 2, 1, 1, 1, 1.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 2, 2, 1, 1, 2.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 2, 3, 1, 1, 1.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 3, 2, 1, 1, 1.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 3, 3, 1, 1, 2.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 3, 4, 1, 1, 1.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 4, 3, 1, 1, 1.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 4, 4, 1, 1, 2.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 4, 5, 1, 1, 1.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 5, 4, 1, 1, 1.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 5, 5, 1, 1, 2.0d0)

    mat%MAT%N = n_node

    call monolis_convert_CSR_to_DIA_R(mat%MAT)

    N = mat%MAT%N
    Ndiag = mat%MAT%DIA%Ndiag

    !# 対角線本数（-1, 0, +1 の 3 本）
    call monolis_test_check_eq_I1("monolis_spmat_convert_dia_test Ndiag", Ndiag, 3)

    !# オフセット
    call monolis_test_check_eq_I1("monolis_spmat_convert_dia_test offset 1", mat%MAT%DIA%offset(1), -1)
    call monolis_test_check_eq_I1("monolis_spmat_convert_dia_test offset 2", mat%MAT%DIA%offset(2),  0)
    call monolis_test_check_eq_I1("monolis_spmat_convert_dia_test offset 3", mat%MAT%DIA%offset(3),  1)

    !# 主対角（offset = 0, d = 2）の値はすべて 2.0
    v = mat%MAT%R%Adia((2-1)*N + 1)
    call monolis_test_check_eq_R1("monolis_spmat_convert_dia_test diag i1", v, 2.0d0)
    v = mat%MAT%R%Adia((2-1)*N + 3)
    call monolis_test_check_eq_R1("monolis_spmat_convert_dia_test diag i3", v, 2.0d0)

    !# 下三角（offset = -1, d = 1）：行 1 は範囲外パディングで 0
    v = mat%MAT%R%Adia((1-1)*N + 1)
    call monolis_test_check_eq_R1("monolis_spmat_convert_dia_test lower pad", v, 0.0d0)
    !# 下三角の行 2 は 1.0
    v = mat%MAT%R%Adia((1-1)*N + 2)
    call monolis_test_check_eq_R1("monolis_spmat_convert_dia_test lower i2", v, 1.0d0)

    !# 上三角（offset = +1, d = 3）：行 5 は範囲外パディングで 0
    v = mat%MAT%R%Adia((3-1)*N + 5)
    call monolis_test_check_eq_R1("monolis_spmat_convert_dia_test upper pad", v, 0.0d0)
    !# 上三角の行 1 は 1.0
    v = mat%MAT%R%Adia((3-1)*N + 1)
    call monolis_test_check_eq_R1("monolis_spmat_convert_dia_test upper i1", v, 1.0d0)

    call monolis_dealloc_DIA_R(mat%MAT)

    call monolis_test_check_eq_I1("monolis_spmat_convert_dia_test Ndiag after dealloc", mat%MAT%DIA%Ndiag, 0)

    call monolis_finalize(mat)

    call monolis_spmat_convert_dia_V_test()
  end subroutine monolis_spmat_convert_dia_test

  !> 可変ブロック DIA 変換テスト
  subroutine monolis_spmat_convert_dia_V_test()
    implicit none
    type(monolis_structure) :: mat
    integer(kint) :: n_node, n_base, n_elem, elem(3,2)
    integer(kint) :: n_dof_list(4)
    integer(kint) :: N, Ndiag, k
    real(kdouble) :: sumA, sumAdia

    call monolis_std_global_log_string("monolis_convert_CSR_to_DIA_V_R")

    if(monolis_mpi_get_global_comm_size() == 2) return

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

    !# 行列値を相異なる決定論的な値で埋める
    do k = 1, size(mat%MAT%R%A)
      mat%MAT%R%A(k) = dble(k)
    enddo

    call monolis_convert_CSR_to_DIA_V_R(mat%MAT)

    N = mat%MAT%N
    Ndiag = mat%MAT%DIA%Ndiag

    !# Vptr の整合性：先頭 0、末尾が値配列長
    call monolis_test_check_eq_I1("monolis_convert_CSR_to_DIA_V Vptr head", mat%MAT%DIA%Vptr(1), 0)
    call monolis_test_check_eq_I1("monolis_convert_CSR_to_DIA_V Vptr tail", &
      mat%MAT%DIA%Vptr(N*Ndiag + 1), size(mat%MAT%R%Adia))

    !# 値の総和は CSR と一致（DIA の余剰スロットは 0 のため）
    sumA = sum(mat%MAT%R%A)
    sumAdia = sum(mat%MAT%R%Adia)
    call monolis_test_check_eq_R1("monolis_convert_CSR_to_DIA_V sum", sumAdia, sumA)

    call monolis_dealloc_DIA_R(mat%MAT)
    call monolis_test_check_eq_I1("monolis_convert_CSR_to_DIA_V Ndiag after dealloc", mat%MAT%DIA%Ndiag, 0)

    call monolis_finalize(mat)
  end subroutine monolis_spmat_convert_dia_V_test

end module mod_monolis_spmat_convert_dia_test
