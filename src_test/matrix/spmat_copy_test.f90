!> 疎行列操作関数群
module mod_monolis_spmat_copy_test
  use mod_monolis
  use mod_monolis_spmat_copy
  use mod_monolis_spmat_nonzero_pattern

  implicit none

contains

  subroutine monolis_spmat_copy_test()
    implicit none
    call monolis_spmat_copy_test_R()
    call monolis_spmat_copy_test_C()
    call monolis_spmat_copy_test_V_R()
  end subroutine monolis_spmat_copy_test

  subroutine monolis_spmat_copy_test_R()
    implicit none
    type(monolis_structure) :: mat_in
    type(monolis_structure) :: mat_out
    integer(kint) :: n_node
    integer(kint) :: n_base
    integer(kint) :: ndof
    integer(kint) :: n_elem
    integer(kint) :: elem(2,3)

    call monolis_std_global_log_string("monolis_copy_mat_R")
    call monolis_std_global_log_string("monolis_copy_mat_nonzero_pattern_R")
    call monolis_std_global_log_string("monolis_copy_mat_nonzero_pattern_val_R")
    call monolis_std_global_log_string("monolis_copy_mat_value_R")
    call monolis_std_global_log_string("monolis_copy_mat_value_matrix_R")
    call monolis_std_global_log_string("monolis_copy_mat_value_rhs_R")
    call monolis_std_global_log_string("monolis_copy_mat_value_solution_R")

    call monolis_std_global_log_string("monolis_copy_mat_nonzero_pattern_CSR")
    call monolis_std_global_log_string("monolis_copy_mat_nonzero_pattern_CSC")
    call monolis_std_global_log_string("monolis_copy_mat_nonzero_pattern_SCSR")

    n_node = 4
    n_base = 2
    ndof = 2
    n_elem = 3
    elem(1,1) = 1; elem(2,1) = 2
    elem(1,2) = 2; elem(2,2) = 3
    elem(1,3) = 3; elem(2,3) = 4

    call monolis_initialize(mat_in)

    call monolis_get_nonzero_pattern_by_simple_mesh_R(mat_in, n_node, n_base, ndof, n_elem, elem)

    mat_in%MAT%R%A = 1.0d0
    mat_in%MAT%R%B = 2.0d0
    mat_in%MAT%R%X = 3.0d0

    call monolis_copy_mat_R(mat_in, mat_out)

    call monolis_test_check_eq_I1("monolis_copy_mat_R 1", mat_in%MAT%N, mat_out%MAT%N)
    call monolis_test_check_eq_I1("monolis_copy_mat_R 2", mat_in%MAT%NP, mat_out%MAT%NP)
    call monolis_test_check_eq_I1("monolis_copy_mat_R 3", mat_in%MAT%NDOF, mat_out%MAT%NDOF)

    call monolis_test_check_eq_I("monolis_copy_mat_R 4", mat_in%MAT%CSR%index, mat_out%MAT%CSR%index)
    call monolis_test_check_eq_I("monolis_copy_mat_R 5", mat_in%MAT%CSR%item, mat_out%MAT%CSR%item)

    call monolis_test_check_eq_I("monolis_copy_mat_R 6", mat_in%MAT%CSC%index, mat_out%MAT%CSC%index)
    call monolis_test_check_eq_I("monolis_copy_mat_R 7", mat_in%MAT%CSC%item, mat_out%MAT%CSC%item)
    call monolis_test_check_eq_I("monolis_copy_mat_R 8", mat_in%MAT%CSC%perm, mat_out%MAT%CSC%perm)

    call monolis_test_check_eq_R("monolis_copy_mat_R 9", mat_in%MAT%R%A, mat_out%MAT%R%A)
    call monolis_test_check_eq_R("monolis_copy_mat_R 10", mat_in%MAT%R%B, mat_out%MAT%R%B)
    call monolis_test_check_eq_R("monolis_copy_mat_R 11", mat_in%MAT%R%X, mat_out%MAT%R%X)

    call monolis_std_global_log_string("monolis_clear_mat_value_R")
    call monolis_std_global_log_string("monolis_clear_mat_value_matrix_R")
    call monolis_std_global_log_string("monolis_clear_mat_value_rhs_R")
    call monolis_std_global_log_string("monolis_clear_mat_value_solution_R")

    call monolis_clear_mat_value_R(mat_out)

    mat_in%MAT%R%A = 0.0d0
    mat_in%MAT%R%B = 0.0d0
    mat_in%MAT%R%X = 0.0d0

    call monolis_test_check_eq_R("monolis_copy_mat_R 12", mat_in%MAT%R%A, mat_out%MAT%R%A)
    call monolis_test_check_eq_R("monolis_copy_mat_R 13", mat_in%MAT%R%B, mat_out%MAT%R%B)
    call monolis_test_check_eq_R("monolis_copy_mat_R 14", mat_in%MAT%R%X, mat_out%MAT%R%X)
  end subroutine monolis_spmat_copy_test_R

  subroutine monolis_spmat_copy_test_C()
    implicit none
    type(monolis_structure) :: mat_in
    type(monolis_structure) :: mat_out
    integer(kint) :: n_node
    integer(kint) :: n_base
    integer(kint) :: ndof
    integer(kint) :: n_elem
    integer(kint) :: elem(2,3)

    call monolis_std_global_log_string("monolis_copy_mat_C")
    call monolis_std_global_log_string("monolis_copy_mat_nonzero_pattern_C")
    call monolis_std_global_log_string("monolis_copy_mat_nonzero_pattern_val_C")
    call monolis_std_global_log_string("monolis_copy_mat_value_C")
    call monolis_std_global_log_string("monolis_copy_mat_value_matrix_C")
    call monolis_std_global_log_string("monolis_copy_mat_value_rhs_C")
    call monolis_std_global_log_string("monolis_copy_mat_value_solution_C")

    n_node = 4
    n_base = 2
    ndof = 2
    n_elem = 3
    elem(1,1) = 1; elem(2,1) = 2
    elem(1,2) = 2; elem(2,2) = 3
    elem(1,3) = 3; elem(2,3) = 4

    call monolis_initialize(mat_in)

    call monolis_get_nonzero_pattern_by_simple_mesh_C(mat_in, n_node, n_base, ndof, n_elem, elem)

    mat_in%MAT%C%A = 1.0d0
    mat_in%MAT%C%B = 2.0d0
    mat_in%MAT%C%X = 3.0d0

    call monolis_copy_mat_C(mat_in, mat_out)

    call monolis_test_check_eq_I1("monolis_copy_mat_C 1", mat_in%MAT%N, mat_out%MAT%N)
    call monolis_test_check_eq_I1("monolis_copy_mat_C 2", mat_in%MAT%NP, mat_out%MAT%NP)
    call monolis_test_check_eq_I1("monolis_copy_mat_C 3", mat_in%MAT%NDOF, mat_out%MAT%NDOF)

    call monolis_test_check_eq_I("monolis_copy_mat_C 4", mat_in%MAT%CSR%index, mat_out%MAT%CSR%index)
    call monolis_test_check_eq_I("monolis_copy_mat_C 5", mat_in%MAT%CSR%item, mat_out%MAT%CSR%item)

    call monolis_test_check_eq_I("monolis_copy_mat_C 6", mat_in%MAT%CSC%index, mat_out%MAT%CSC%index)
    call monolis_test_check_eq_I("monolis_copy_mat_C 7", mat_in%MAT%CSC%item, mat_out%MAT%CSC%item)
    call monolis_test_check_eq_I("monolis_copy_mat_C 8", mat_in%MAT%CSC%perm, mat_out%MAT%CSC%perm)

    call monolis_test_check_eq_C("monolis_copy_mat_C 9", mat_in%MAT%C%A, mat_out%MAT%C%A)
    call monolis_test_check_eq_C("monolis_copy_mat_C 10", mat_in%MAT%C%B, mat_out%MAT%C%B)
    call monolis_test_check_eq_C("monolis_copy_mat_C 11", mat_in%MAT%C%X, mat_out%MAT%C%X)

    call monolis_std_global_log_string("monolis_clear_mat_value_C")
    call monolis_std_global_log_string("monolis_clear_mat_value_matrix_C")
    call monolis_std_global_log_string("monolis_clear_mat_value_rhs_C")
    call monolis_std_global_log_string("monolis_clear_mat_value_solution_C")

    call monolis_clear_mat_value_C(mat_out)

    mat_in%MAT%C%A = 0.0d0
    mat_in%MAT%C%B = 0.0d0
    mat_in%MAT%C%X = 0.0d0

    call monolis_test_check_eq_C("monolis_copy_mat_C 12", mat_in%MAT%C%A, mat_out%MAT%C%A)
    call monolis_test_check_eq_C("monolis_copy_mat_C 13", mat_in%MAT%C%B, mat_out%MAT%C%B)
    call monolis_test_check_eq_C("monolis_copy_mat_C 14", mat_in%MAT%C%X, mat_out%MAT%C%X)
  end subroutine monolis_spmat_copy_test_C

  subroutine monolis_spmat_copy_test_V_R()
    implicit none
    type(monolis_structure) :: mat_in
    type(monolis_structure) :: mat_out
    integer(kint) :: n_node
    integer(kint) :: n_base
    integer(kint) :: n_dof_list(4)
    integer(kint) :: n_elem
    integer(kint) :: elem(3,2)

    call monolis_std_global_log_string("monolis_copy_mat_R")

    ! 4節点、2要素のメッシュ設定
    n_node = 4
    n_base = 3
    n_elem = 2

    ! 三角形要素
    elem(1,1) = 1; elem(2,1) = 2; elem(3,1) = 3
    elem(1,2) = 2; elem(2,2) = 3; elem(3,2) = 4

    ! 各節点の自由度を異なる値に設定
    ! 節点1: 1自由度, 節点2: 2自由度, 節点3: 3自由度, 節点4: 1自由度
    n_dof_list(1) = 1
    n_dof_list(2) = 2
    n_dof_list(3) = 3
    n_dof_list(4) = 1

    call monolis_initialize(mat_in)

    ! monolis_get_nonzero_pattern_by_simple_mesh_V_Rで非ゼロ構造を決定
    call monolis_get_nonzero_pattern_by_simple_mesh_V_R(mat_in, n_node, n_base, n_dof_list, n_elem, elem)

    ! テスト用の値を設定
    mat_in%MAT%R%A = 1.0d0
    mat_in%MAT%R%B = 2.0d0
    mat_in%MAT%R%X = 3.0d0

    ! 行列をコピー
    call monolis_copy_mat_R(mat_in, mat_out)

    ! 基本パラメータのテスト
    call monolis_test_check_eq_I1("monolis_copy_mat_V_R 1", mat_in%MAT%N, mat_out%MAT%N)
    call monolis_test_check_eq_I1("monolis_copy_mat_V_R 2", mat_in%MAT%NP, mat_out%MAT%NP)
    call monolis_test_check_eq_I1("monolis_copy_mat_V_R 3", mat_in%MAT%NDOF, mat_out%MAT%NDOF)

    ! CSR構造のテスト
    call monolis_test_check_eq_I("monolis_copy_mat_V_R 4", mat_in%MAT%CSR%index, mat_out%MAT%CSR%index)
    call monolis_test_check_eq_I("monolis_copy_mat_V_R 5", mat_in%MAT%CSR%item, mat_out%MAT%CSR%item)

    ! CSC構造のテスト
    call monolis_test_check_eq_I("monolis_copy_mat_V_R 6", mat_in%MAT%CSC%index, mat_out%MAT%CSC%index)
    call monolis_test_check_eq_I("monolis_copy_mat_V_R 7", mat_in%MAT%CSC%item, mat_out%MAT%CSC%item)
    call monolis_test_check_eq_I("monolis_copy_mat_V_R 8", mat_in%MAT%CSC%perm, mat_out%MAT%CSC%perm)

    ! 自由度インデックスのテスト
    call monolis_test_check_eq_I("monolis_copy_mat_V_R 9", mat_in%MAT%n_dof_index, mat_out%MAT%n_dof_index)
    call monolis_test_check_eq_I("monolis_copy_mat_V_R 9", mat_in%MAT%n_dof_index2, mat_out%MAT%n_dof_index2)
    call monolis_test_check_eq_I("monolis_copy_mat_V_R 9", mat_in%MAT%n_dof_list, mat_out%MAT%n_dof_list)

    ! 行列値のテスト
    call monolis_test_check_eq_R("monolis_copy_mat_V_R 10", mat_in%MAT%R%A, mat_out%MAT%R%A)
    call monolis_test_check_eq_R("monolis_copy_mat_V_R 11", mat_in%MAT%R%B, mat_out%MAT%R%B)
    call monolis_test_check_eq_R("monolis_copy_mat_V_R 12", mat_in%MAT%R%X, mat_out%MAT%R%X)

    ! 値をクリアしてテスト
    call monolis_clear_mat_value_R(mat_out)

    mat_in%MAT%R%A = 0.0d0
    mat_in%MAT%R%B = 0.0d0
    mat_in%MAT%R%X = 0.0d0

    call monolis_test_check_eq_R("monolis_copy_mat_V_R 13", mat_in%MAT%R%A, mat_out%MAT%R%A)
    call monolis_test_check_eq_R("monolis_copy_mat_V_R 14", mat_in%MAT%R%B, mat_out%MAT%R%B)
    call monolis_test_check_eq_R("monolis_copy_mat_V_R 15", mat_in%MAT%R%X, mat_out%MAT%R%X)
  end subroutine monolis_spmat_copy_test_V_R
end module mod_monolis_spmat_copy_test
