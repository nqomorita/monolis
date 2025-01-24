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
    call monolis_convert_sparse_matrix_to_dense_matrix_R_test()
  end subroutine monolis_spmat_handler_util_test

  subroutine monolis_set_scalar_to_sparse_matrix_main_R_test()
    implicit none
    integer(kint) :: index(5)
    integer(kint) :: item(10)
    real(kdouble) :: A(40)
    integer(kint) :: ndof
    integer(kint) :: ci
    integer(kint) :: cj
    integer(kint) :: csub_i
    integer(kint) :: csub_j
    real(kdouble) :: val
    logical :: is_find

    call monolis_std_global_log_string("monolis_set_scalar_to_sparse_matrix_main_R")
    call monolis_std_global_log_string("monolis_get_scalar_from_sparse_matrix_main_R")
    call monolis_std_global_log_string("monolis_add_scalar_to_sparse_matrix_main_R")

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
    integer(kint) :: index(5)
    integer(kint) :: item(10)
    complex(kdouble) :: A(40)
    integer(kint) :: ndof
    integer(kint) :: ci
    integer(kint) :: cj
    integer(kint) :: csub_i
    integer(kint) :: csub_j
    complex(kdouble) :: val
    logical :: is_find

    call monolis_std_global_log_string("monolis_set_scalar_to_sparse_matrix_main_C")
    call monolis_std_global_log_string("monolis_get_scalar_from_sparse_matrix_main_C")
    call monolis_std_global_log_string("monolis_add_scalar_to_sparse_matrix_main_C")

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
    integer(kint) :: index(5)
    integer(kint) :: item(10)
    real(kdouble) :: A(40)
    integer(kint) :: ndof
    integer(kint) :: ci
    integer(kint) :: cj
    real(kdouble) :: val(2,2)
    real(kdouble) :: B(5)
    integer(kint) :: e(1)

    call monolis_std_global_log_string("monolis_set_block_to_sparse_matrix_main_R")

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

  subroutine monolis_set_block_to_sparse_matrix_main_C_test()
    implicit none
    integer(kint) :: index(5)
    integer(kint) :: item(10)
    complex(kdouble) :: A(40)
    integer(kint) :: ndof
    integer(kint) :: ci
    integer(kint) :: cj
    complex(kdouble) :: val(2,2)
    integer(kint) :: e(1)

    call monolis_std_global_log_string("monolis_set_block_to_sparse_matrix_main_C")

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

  subroutine monolis_convert_sparse_matrix_to_dense_matrix_R_test
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: mat
    !> monolis 通信構造体
    type(monolis_com) :: com
    !> 疎行列を変換した密行列
    real(kdouble), allocatable :: dense(:,:)
    integer(kint) :: nnode, nelem, elem(2,3), n_dof
    integer(kint) :: i
    real(kdouble) :: ans(8)

    if(monolis_mpi_get_global_comm_size() /= 1) return

    call monolis_std_global_log_string("monolis_convert_sparse_matrix_to_dense_matrix_R_test")
    call monolis_std_global_log_string("monolis_convert_sparse_matrix_to_dense_matrix_R")

    call monolis_initialize(mat)
    call monolis_com_initialize_by_self(com)

    nnode = 4
    nelem = 3
    n_dof = 2

    elem(1,1) = 1; elem(2,1) = 2;
    elem(1,2) = 2; elem(2,2) = 3;
    elem(1,3) = 3; elem(2,3) = 4;

    call monolis_get_nonzero_pattern_by_simple_mesh_R(mat, nnode, 2, n_dof, nelem, elem)

    do i = 1, 4
      call monolis_add_scalar_to_sparse_matrix_R(mat, i, i, 1, 1, 10.0d0*dble(i) + 1.0d0)
      call monolis_add_scalar_to_sparse_matrix_R(mat, i, i, 2, 2, 10.0d0*dble(i) + 4.0d0)
      call monolis_add_scalar_to_sparse_matrix_R(mat, i, i, 1, 2, 10.0d0*dble(i) + 2.0d0)
      call monolis_add_scalar_to_sparse_matrix_R(mat, i, i, 2, 1, 10.0d0*dble(i) + 3.0d0)
    enddo

    do i = 1, 3
      call monolis_add_scalar_to_sparse_matrix_R(mat, i, i + 1, 1, 1, 100.0d0*dble(i) + 10.0d0*dble(i+1) + 1.0d0)
      call monolis_add_scalar_to_sparse_matrix_R(mat, i, i + 1, 2, 2, 100.0d0*dble(i) + 10.0d0*dble(i+1) + 4.0d0)
      call monolis_add_scalar_to_sparse_matrix_R(mat, i, i + 1, 1, 2, 100.0d0*dble(i) + 10.0d0*dble(i+1) + 2.0d0)
      call monolis_add_scalar_to_sparse_matrix_R(mat, i, i + 1, 2, 1, 100.0d0*dble(i) + 10.0d0*dble(i+1) + 3.0d0)
      call monolis_add_scalar_to_sparse_matrix_R(mat, i + 1, i, 1, 1, 100.0d0*dble(i+1) + 10.0d0*dble(i) + 1.0d0)
      call monolis_add_scalar_to_sparse_matrix_R(mat, i + 1, i, 2, 2, 100.0d0*dble(i+1) + 10.0d0*dble(i) + 4.0d0)
      call monolis_add_scalar_to_sparse_matrix_R(mat, i + 1, i, 1, 2, 100.0d0*dble(i+1) + 10.0d0*dble(i) + 2.0d0)
      call monolis_add_scalar_to_sparse_matrix_R(mat, i + 1, i, 2, 1, 100.0d0*dble(i+1) + 10.0d0*dble(i) + 3.0d0)
    enddo

    call monolis_convert_sparse_matrix_to_dense_matrix_R(mat%MAT, com, dense)

    ans = 0.0d0
    ans(1) = 11.0d0; ans(2) = 12.0d0; ans(3) = 121.0d0; ans(4) = 122.0d0;
    call monolis_test_check_eq_R("monolis_convert_sparse_matrix_to_dense_matrix_R_test 1", dense(1,:), ans)

    ans = 0.0d0
    ans(1) = 13.0d0; ans(2) = 14.0d0; ans(3) = 123.0d0; ans(4) = 124.0d0;
    call monolis_test_check_eq_R("monolis_convert_sparse_matrix_to_dense_matrix_R_test 2", dense(2,:), ans)

    ans = 0.0d0
    ans(1) = 211.0d0; ans(2) = 212.0d0; ans(3) = 21.0d0; ans(4) = 22.0d0; ans(5) = 231.0d0; ans(6) = 232.0d0;
    call monolis_test_check_eq_R("monolis_convert_sparse_matrix_to_dense_matrix_R_test 3", dense(3,:), ans)

    ans = 0.0d0
    ans(1) = 213.0d0; ans(2) = 214.0d0; ans(3) = 23.0d0; ans(4) = 24.0d0; ans(5) = 233.0d0; ans(6) = 234.0d0;
    call monolis_test_check_eq_R("monolis_convert_sparse_matrix_to_dense_matrix_R_test 4", dense(4,:), ans)

    ans = 0.0d0
    ans(3) = 321.0d0; ans(4) = 322.0d0; ans(5) = 31.0d0; ans(6) = 32.0d0; ans(7) = 341.0d0; ans(8) = 342.0d0;
    call monolis_test_check_eq_R("monolis_convert_sparse_matrix_to_dense_matrix_R_test 5", dense(5,:), ans)

    ans = 0.0d0
    ans(3) = 323.0d0; ans(4) = 324.0d0; ans(5) = 33.0d0; ans(6) = 34.0d0; ans(7) = 343.0d0; ans(8) = 344.0d0;
    call monolis_test_check_eq_R("monolis_convert_sparse_matrix_to_dense_matrix_R_test 6", dense(6,:), ans)

    ans = 0.0d0
    ans(5) = 431.0d0; ans(6) = 432.0d0; ans(7) = 41.0d0; ans(8) = 42.0d0;
    call monolis_test_check_eq_R("monolis_convert_sparse_matrix_to_dense_matrix_R_test 7", dense(7,:), ans)

    ans = 0.0d0
    ans(5) = 433.0d0; ans(6) = 434.0d0; ans(7) = 43.0d0; ans(8) = 44.0d0;
    call monolis_test_check_eq_R("monolis_convert_sparse_matrix_to_dense_matrix_R_test 8", dense(8,:), ans)

  end subroutine monolis_convert_sparse_matrix_to_dense_matrix_R_test
end module mod_monolis_spmat_handler_util_test
