!> 疎行列積テストモジュール
module mod_monolis_matmat_test
  use mod_monolis
  use mod_monolis_utils
  implicit none

contains

  subroutine monolis_matmat_test()
    implicit none

    call monolis_matmat_symbolic_test()
    call monolis_matmat_value_nn_test()
    call monolis_matmat_nn_test()
  end subroutine monolis_matmat_test

  !> 疎行列積のシンボリック計算テスト
  subroutine monolis_matmat_symbolic_test()
    implicit none
    type(monolis_structure) :: A, B, C
    type(monolis_com) :: com
    integer(kint) :: i, n_node, n_dof, n_elem, n_base
    integer(kint), allocatable :: elem(:,:)

    call monolis_std_global_log_string("monolis_matmat_symbolic_test")

    !# 簡単な3節点の1次元要素メッシュを作成
    !# 1--2--3
    n_node = 3
    n_dof = 1
    n_elem = 2
    n_base = 2

    call monolis_initialize(A)
    call monolis_initialize(B)
    call monolis_initialize(C)
    !call monolis_com_initialize(com)

    call monolis_alloc_I_2d(elem, n_elem, n_base)
    elem(1,1) = 1; elem(1,2) = 2
    elem(2,1) = 2; elem(2,2) = 3

    !# 行列 A と B の非零パターンを作成
    call monolis_get_nonzero_pattern_by_simple_mesh_R(A, n_node, n_base, n_dof, n_elem, elem)
    call monolis_get_nonzero_pattern_by_simple_mesh_R(B, n_node, n_base, n_dof, n_elem, elem)

    !# シンボリック計算
    call monolis_matmat_symbolic(com, A%MAT, B%MAT, C%MAT, n_dof)

    !# 結果の検証
    call monolis_test_check_eq_I1("monolis_matmat_symbolic_test 1", C%MAT%N, n_node)
    call monolis_test_check_eq_I1("monolis_matmat_symbolic_test 2", C%MAT%NP, n_node)
    call monolis_test_check_eq_I1("monolis_matmat_symbolic_test 3", C%MAT%NDOF, n_dof)

    !# CSR%index の検証（各行の非零要素数）
    call monolis_test_check_eq_I1("monolis_matmat_symbolic_test 4", C%MAT%CSR%index(1), 0)
    call monolis_test_check_eq_I1("monolis_matmat_symbolic_test 5", C%MAT%CSR%index(2), 3)
    call monolis_test_check_eq_I1("monolis_matmat_symbolic_test 6", C%MAT%CSR%index(3), 6)
    call monolis_test_check_eq_I1("monolis_matmat_symbolic_test 7", C%MAT%CSR%index(4), 9)

    call monolis_test_check_eq_I1("monolis_matmat_symbolic_test 8a", C%MAT%CSR%item(1), 1)
    call monolis_test_check_eq_I1("monolis_matmat_symbolic_test 8b", C%MAT%CSR%item(2), 2)
    call monolis_test_check_eq_I1("monolis_matmat_symbolic_test 8c", C%MAT%CSR%item(3), 3)
    call monolis_test_check_eq_I1("monolis_matmat_symbolic_test 8d", C%MAT%CSR%item(4), 1)
    call monolis_test_check_eq_I1("monolis_matmat_symbolic_test 8e", C%MAT%CSR%item(5), 2)
    call monolis_test_check_eq_I1("monolis_matmat_symbolic_test 8f", C%MAT%CSR%item(6), 3)
    call monolis_test_check_eq_I1("monolis_matmat_symbolic_test 8g", C%MAT%CSR%item(7), 1)
    call monolis_test_check_eq_I1("monolis_matmat_symbolic_test 8h", C%MAT%CSR%item(8), 2)
    call monolis_test_check_eq_I1("monolis_matmat_symbolic_test 8i", C%MAT%CSR%item(9), 3)

    call monolis_finalize(A)
    call monolis_finalize(B)
    call monolis_finalize(C)
  end subroutine monolis_matmat_symbolic_test

  !> 疎行列積の数値計算テスト
  subroutine monolis_matmat_value_nn_test()
    implicit none
    type(monolis_structure) :: A, B, C
    type(monolis_com) :: com
    integer(kint) :: i, n_node, n_dof, n_elem, n_base
    integer(kint), allocatable :: elem(:,:)
    real(kdouble) :: expected

    call monolis_std_global_log_string("monolis_matmat_value_nn_test")

    !# 簡単な3節点の1次元要素メッシュを作成
    n_node = 3
    n_dof = 1
    n_elem = 2
    n_base = 2

    call monolis_initialize(A)
    call monolis_initialize(B)
    call monolis_initialize(C)
    !call monolis_com_initialize(com)

    call monolis_alloc_I_2d(elem, n_elem, n_base)
    elem(1,1) = 1; elem(1,2) = 2
    elem(2,1) = 2; elem(2,2) = 3

    !# 行列 A と B の非零パターンを作成
    call monolis_get_nonzero_pattern_by_simple_mesh_R(A, n_node, n_base, n_dof, n_elem, elem)
    call monolis_get_nonzero_pattern_by_simple_mesh_R(B, n_node, n_base, n_dof, n_elem, elem)

    !# A に値を設定（対角2、非対角1）
    A%MAT%R%A(1) = 2.0d0
    A%MAT%R%A(2) = 1.0d0
    A%MAT%R%A(3) = 1.0d0
    A%MAT%R%A(4) = 2.0d0
    A%MAT%R%A(5) = 1.0d0
    A%MAT%R%A(6) = 1.0d0
    A%MAT%R%A(7) = 2.0d0

    !# B に値を設定（単位行列）
    B%MAT%R%A(1) = 2.0d0
    B%MAT%R%A(2) = 1.0d0
    B%MAT%R%A(3) = 1.0d0
    B%MAT%R%A(4) = 2.0d0
    B%MAT%R%A(5) = 1.0d0
    B%MAT%R%A(6) = 1.0d0
    B%MAT%R%A(7) = 2.0d0

    !# シンボリック計算と数値計算
    call monolis_matmat_symbolic(com, A%MAT, B%MAT, C%MAT, n_dof)
    call monolis_matmat_value_nn(com, A%MAT, B%MAT, C%MAT, n_dof)

    call monolis_test_check_eq_R1("monolis_matmat_value_nn_test 1", C%MAT%R%A(1), 5.0d0)
    call monolis_test_check_eq_R1("monolis_matmat_value_nn_test 2", C%MAT%R%A(2), 4.0d0)
    call monolis_test_check_eq_R1("monolis_matmat_value_nn_test 3", C%MAT%R%A(3), 1.0d0)
    call monolis_test_check_eq_R1("monolis_matmat_value_nn_test 4", C%MAT%R%A(4), 4.0d0)
    call monolis_test_check_eq_R1("monolis_matmat_value_nn_test 5", C%MAT%R%A(5), 6.0d0)
    call monolis_test_check_eq_R1("monolis_matmat_value_nn_test 6", C%MAT%R%A(6), 4.0d0)
    call monolis_test_check_eq_R1("monolis_matmat_value_nn_test 7", C%MAT%R%A(7), 1.0d0)
    call monolis_test_check_eq_R1("monolis_matmat_value_nn_test 8", C%MAT%R%A(8), 4.0d0)
    call monolis_test_check_eq_R1("monolis_matmat_value_nn_test 9", C%MAT%R%A(9), 5.0d0)

    call monolis_finalize(A)
    call monolis_finalize(B)
    call monolis_finalize(C)
  end subroutine monolis_matmat_value_nn_test

  !> 疎行列積の統合テスト
  subroutine monolis_matmat_nn_test()
    implicit none
    type(monolis_structure) :: A, B, C
    type(monolis_com) :: com
    integer(kint) :: i, j, n_node, n_dof, n_elem, n_base
    integer(kint), allocatable :: elem(:,:)
    real(kdouble) :: tol

    call monolis_std_global_log_string("monolis_matmat_nn_test")

    !# 簡単な3節点の1次元要素メッシュを作成
    n_node = 3
    n_dof = 2
    n_elem = 2
    n_base = 2

    call monolis_initialize(A)
    call monolis_initialize(B)
    call monolis_initialize(C)
    !call monolis_com_initialize(com)

    call monolis_alloc_I_2d(elem, n_elem, n_base)
    elem(1,1) = 1; elem(1,2) = 2
    elem(2,1) = 1; elem(2,2) = 3

    call monolis_get_nonzero_pattern_by_simple_mesh_R(A, n_node, n_base, n_dof, n_elem, elem)
    call monolis_get_nonzero_pattern_by_simple_mesh_R(B, n_node, n_base, n_dof, n_elem, elem)

    A%MAT%R%A = 1.0d0
    A%MAT%R%A(1) = 1.0d0
    A%MAT%R%A(2) = 2.0d0
    A%MAT%R%A(3) = 3.0d0
    A%MAT%R%A(4) = 4.0d0
    B%MAT%R%A = 1.0d0

    call monolis_matmat_nn(com, A%MAT, B%MAT, C%MAT, n_dof)

    call monolis_test_check_eq_I1("monolis_matmat_nn_test 1", C%MAT%CSR%index(1), 0)
    call monolis_test_check_eq_I1("monolis_matmat_nn_test 2", C%MAT%CSR%index(2), 1)
    call monolis_test_check_eq_I1("monolis_matmat_nn_test 3", C%MAT%CSR%index(3), 3)
    call monolis_test_check_eq_I1("monolis_matmat_nn_test 4", C%MAT%CSR%index(4), 5)

    call monolis_test_check_eq_I1("monolis_matmat_nn_test 5", C%MAT%CSR%item(1), 1)
    call monolis_test_check_eq_I1("monolis_matmat_nn_test 6", C%MAT%CSR%item(2), 2)
    call monolis_test_check_eq_I1("monolis_matmat_nn_test 7", C%MAT%CSR%item(3), 3)
    call monolis_test_check_eq_I1("monolis_matmat_nn_test 8", C%MAT%CSR%item(4), 2)
    call monolis_test_check_eq_I1("monolis_matmat_nn_test 9", C%MAT%CSR%item(5), 3)

    call monolis_test_check_eq_R1("monolis_matmat_nn_test a1", C%MAT%R%A(1), 3.0d0)
    call monolis_test_check_eq_R1("monolis_matmat_nn_test a2", C%MAT%R%A(2), 3.0d0)
    call monolis_test_check_eq_R1("monolis_matmat_nn_test a3", C%MAT%R%A(3), 7.0d0)
    call monolis_test_check_eq_R1("monolis_matmat_nn_test a4", C%MAT%R%A(4), 7.0d0)
    call monolis_test_check_eq_R1("monolis_matmat_nn_test a5", C%MAT%R%A(5), 4.0d0)
    call monolis_test_check_eq_R1("monolis_matmat_nn_test a5", C%MAT%R%A(6), 4.0d0)
    call monolis_test_check_eq_R1("monolis_matmat_nn_test a5", C%MAT%R%A(7), 4.0d0)
    call monolis_test_check_eq_R1("monolis_matmat_nn_test a5", C%MAT%R%A(8), 4.0d0)

    call monolis_finalize(A)
    call monolis_finalize(B)
    call monolis_finalize(C)
  end subroutine monolis_matmat_nn_test

end module mod_monolis_matmat_test
