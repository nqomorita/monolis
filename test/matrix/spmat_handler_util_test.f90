!> 疎行列操作関数群テストモジュール
module mod_monolis_spmat_handler_util_test
  use mod_monolis
  use mod_monolis_spmat_handler_util

  implicit none

contains

  subroutine monolis_spmat_handler_util_test()
    implicit none

    call monolis_set_scalar_to_sparse_matrix_main_R_test()
    call monolis_set_scalar_to_sparse_matrix_main_C_test()
    call monolis_set_block_to_sparse_matrix_main_R_test()
    call monolis_set_block_to_sparse_matrix_main_C_test()
  end subroutine monolis_spmat_handler_util_test

  subroutine monolis_set_scalar_to_sparse_matrix_main_R_test()
    implicit none
    !> index 配列
    integer(kint) :: index(5)
    !> item 配列
    integer(kint) :: item(10)
    !> 係数行列
    real(kdouble) :: A(40)
    !> ブロック自由度
    integer(kint) :: ndof
    !> 行番号
    integer(kint) :: ci
    !> 列番号
    integer(kint) :: cj
    !> ブロック中の行番号
    integer(kint) :: csub_i
    !> ブロック中の列番号
    integer(kint) :: csub_j
    !> 設定値
    real(kdouble) :: val
    logical :: is_find

    call monolis_std_log_string("monolis_set_scalar_to_sparse_matrix_main_R_test")

    ndof = 2

    index(1) = 0
    index(2) = 2
    index(3) = 5
    index(4) = 8
    index(5) = 10

    item(1) = 1
    item(2) = 2
    item(3) = 1
    item(4) = 2
    item(5) = 3
    item(6) = 2
    item(7) = 3
    item(8) = 4
    item(9) = 3
    item(10) = 4

    A = 0.0d0

    call monolis_set_scalar_to_sparse_matrix_main_R(index, item, A, ndof, 2, 2, 1, 1, 1.0d0)

    call monolis_test_check_eq_R1("monolis_set_scalar_to_sparse_matrix_main_R_test 1", A(13), 1.0d0)

    val = 0.0d0

    call monolis_get_scalar_from_sparse_matrix_main_R(index, item, A, ndof, 2, 2, 1, 1, val, is_find)

    call monolis_test_check_eq_R1("monolis_set_scalar_to_sparse_matrix_main_R_test 2", val, 1.0d0)

    call monolis_add_scalar_to_sparse_matrix_main_R(index, item, A, ndof, 2, 2, 1, 1, 1.0d0)

    val = 0.0d0

    call monolis_get_scalar_from_sparse_matrix_main_R(index, item, A, ndof, 2, 2, 1, 1, val, is_find)

    call monolis_test_check_eq_R1("monolis_set_scalar_to_sparse_matrix_main_R_test 3", val, 2.0d0)
  end subroutine monolis_set_scalar_to_sparse_matrix_main_R_test

  subroutine monolis_set_scalar_to_sparse_matrix_main_C_test()
    implicit none
    !> index 配列
    integer(kint) :: index(5)
    !> item 配列
    integer(kint) :: item(10)
    !> 係数行列
    complex(kdouble) :: A(40)
    !> ブロック自由度
    integer(kint) :: ndof
    !> 行番号
    integer(kint) :: ci
    !> 列番号
    integer(kint) :: cj
    !> ブロック中の行番号
    integer(kint) :: csub_i
    !> ブロック中の列番号
    integer(kint) :: csub_j
    !> 設定値
    complex(kdouble) :: val
    logical :: is_find

    call monolis_std_log_string("monolis_set_scalar_to_sparse_matrix_main_C_test")

    ndof = 2

    index(1) = 0
    index(2) = 2
    index(3) = 5
    index(4) = 8
    index(5) = 10

    item(1) = 1
    item(2) = 2
    item(3) = 1
    item(4) = 2
    item(5) = 3
    item(6) = 2
    item(7) = 3
    item(8) = 4
    item(9) = 3
    item(10) = 4

    A = (0.0d0, 0.0d0)

    call monolis_set_scalar_to_sparse_matrix_main_C(index, item, A, ndof, 2, 2, 1, 1, (1.0d0, 2.0d0))

    call monolis_test_check_eq_C1("monolis_set_scalar_to_sparse_matrix_main_C_test 1", A(13), (1.0d0, 2.0d0))

    val = (0.0d0, 0.0d0)

    call monolis_get_scalar_from_sparse_matrix_main_C(index, item, A, ndof, 2, 2, 1, 1, val, is_find)

    call monolis_test_check_eq_C1("monolis_set_scalar_to_sparse_matrix_main_C_test 2", val, (1.0d0, 2.0d0))

    call monolis_add_scalar_to_sparse_matrix_main_C(index, item, A, ndof, 2, 2, 1, 1, (1.0d0, 2.0d0))

    val = (0.0d0, 0.0d0)

    call monolis_get_scalar_from_sparse_matrix_main_C(index, item, A, ndof, 2, 2, 1, 1, val, is_find)

    call monolis_test_check_eq_C1("monolis_set_scalar_to_sparse_matrix_main_C_test 3", val, (2.0d0, 4.0d0))
  end subroutine monolis_set_scalar_to_sparse_matrix_main_C_test

  subroutine monolis_set_block_to_sparse_matrix_main_R_test()
    implicit none
    !> index 配列
    integer(kint) :: index(5)
    !> item 配列
    integer(kint) :: item(10)
    !> 係数行列
    real(kdouble) :: A(40)
    !> ブロック自由度
    integer(kint) :: ndof
    !> 行番号
    integer(kint) :: ci
    !> 列番号
    integer(kint) :: cj
    !> 設定値
    real(kdouble) :: val(2,2)
    integer(kint) :: e(1)

    call monolis_std_log_string("monolis_set_block_to_sparse_matrix_main_R_test")

    ndof = 2

    index(1) = 0
    index(2) = 2
    index(3) = 5
    index(4) = 8
    index(5) = 10

    item(1) = 1
    item(2) = 2
    item(3) = 1
    item(4) = 2
    item(5) = 3
    item(6) = 2
    item(7) = 3
    item(8) = 4
    item(9) = 3
    item(10) = 4

    A = 0.0d0

    val(1,1) = 1.0d0
    val(1,2) = 2.0d0
    val(2,1) = 3.0d0
    val(2,2) = 4.0d0

    call monolis_set_block_to_sparse_matrix_main_R(index, item, A, ndof, 2, 2, val)

    call monolis_test_check_eq_R1("monolis_set_block_to_sparse_matrix_main_R_test 1", A(13), 1.0d0)
    call monolis_test_check_eq_R1("monolis_set_block_to_sparse_matrix_main_R_test 2", A(14), 2.0d0)
    call monolis_test_check_eq_R1("monolis_set_block_to_sparse_matrix_main_R_test 3", A(15), 3.0d0)
    call monolis_test_check_eq_R1("monolis_set_block_to_sparse_matrix_main_R_test 4", A(16), 4.0d0)

    val(1,1) = 1.0d0
    val(1,2) = 2.0d0
    val(2,1) = 3.0d0
    val(2,2) = 4.0d0

    e(1) = 2

    call monolis_add_matrix_to_sparse_matrix_main_R(index, item, A, 1, 1, ndof, e, e, val)

    call monolis_test_check_eq_R1("monolis_set_block_to_sparse_matrix_main_R_test 1", A(13), 2.0d0)
    call monolis_test_check_eq_R1("monolis_set_block_to_sparse_matrix_main_R_test 2", A(14), 4.0d0)
    call monolis_test_check_eq_R1("monolis_set_block_to_sparse_matrix_main_R_test 3", A(15), 6.0d0)
    call monolis_test_check_eq_R1("monolis_set_block_to_sparse_matrix_main_R_test 4", A(16), 8.0d0)
  end subroutine monolis_set_block_to_sparse_matrix_main_R_test

    !call monolis_add_matrix_to_sparse_matrix_main_R(index, item, A, n1, n2, ndof, e1, e2, val)

  subroutine monolis_set_block_to_sparse_matrix_main_C_test()
    implicit none
    !> index 配列
    integer(kint) :: index(5)
    !> item 配列
    integer(kint) :: item(10)
    !> 係数行列
    complex(kdouble) :: A(40)
    !> ブロック自由度
    integer(kint) :: ndof
    !> 行番号
    integer(kint) :: ci
    !> 列番号
    integer(kint) :: cj
    !> 設定値
    complex(kdouble) :: val(2,2)
    integer(kint) :: e(1)

    call monolis_std_log_string("monolis_set_block_to_sparse_matrix_main_C_test")

    ndof = 2

    index(1) = 0
    index(2) = 2
    index(3) = 5
    index(4) = 8
    index(5) = 10

    item(1) = 1
    item(2) = 2
    item(3) = 1
    item(4) = 2
    item(5) = 3
    item(6) = 2
    item(7) = 3
    item(8) = 4
    item(9) = 3
    item(10) = 4

    A = 0.0d0

    val(1,1) = (1.0d0, 1.0d0)
    val(1,2) = (2.0d0, 2.0d0)
    val(2,1) = (3.0d0, 3.0d0)
    val(2,2) = (4.0d0, 4.0d0)

    call monolis_set_block_to_sparse_matrix_main_C(index, item, A, ndof, 2, 2, val)

    call monolis_test_check_eq_C1("monolis_set_block_to_sparse_matrix_main_C_test 1", A(13), (1.0d0, 1.0d0))
    call monolis_test_check_eq_C1("monolis_set_block_to_sparse_matrix_main_C_test 2", A(14), (2.0d0, 2.0d0))
    call monolis_test_check_eq_C1("monolis_set_block_to_sparse_matrix_main_C_test 3", A(15), (3.0d0, 3.0d0))
    call monolis_test_check_eq_C1("monolis_set_block_to_sparse_matrix_main_C_test 4", A(16), (4.0d0, 4.0d0))

    e(1) = 2

    call monolis_add_matrix_to_sparse_matrix_main_C(index, item, A, 1, 1, ndof, e, e, val)

    call monolis_test_check_eq_C1("monolis_set_block_to_sparse_matrix_main_C_test 1", A(13), (2.0d0, 2.0d0))
    call monolis_test_check_eq_C1("monolis_set_block_to_sparse_matrix_main_C_test 2", A(14), (4.0d0, 4.0d0))
    call monolis_test_check_eq_C1("monolis_set_block_to_sparse_matrix_main_C_test 3", A(15), (6.0d0, 6.0d0))
    call monolis_test_check_eq_C1("monolis_set_block_to_sparse_matrix_main_C_test 4", A(16), (8.0d0, 8.0d0))
  end subroutine monolis_set_block_to_sparse_matrix_main_C_test

    !call monolis_set_Dirichlet_bc_main_R(index, item, A, B, indexR, itemR, permA, &
    !call monolis_set_Dirichlet_bc_main_C(index, item, A, B, indexR, itemR, permA, &
    !call monolis_get_max_matrix_component_main_R(monoMAT, monoCOM, max_val)
    !call monolis_check_diagonal_zero_component_main_R(monoPRM, monoMAT)

end module mod_monolis_spmat_handler_util_test
