!> 疎行列 ELL 変換関数群テスト
module mod_monolis_spmat_convert_ell_test
  use mod_monolis
  use mod_monolis_spmat_convert_ell
  use mod_monolis_spmat_nonzero_pattern

  implicit none

contains

  subroutine monolis_spmat_convert_ell_test()
    implicit none
    type(monolis_structure) :: mat
    integer(kint) :: n_node, n_elem, elem(2,4)
    integer(kint) :: N, Nmaxcol
    real(kdouble) :: v

    call monolis_std_global_log_string("monolis_convert_CSR_to_ELL_R")
    call monolis_std_global_log_string("monolis_dealloc_ELL_R")

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

    call monolis_convert_CSR_to_ELL_R(mat%MAT)

    N = mat%MAT%N
    Nmaxcol = mat%MAT%ELL%Nmaxcol

    !# 1 行あたりの最大非ゼロブロック数（行 2,3,4 で 3）
    call monolis_test_check_eq_I1("monolis_spmat_convert_ell_test Nmaxcol", Nmaxcol, 3)

    !# 列番号：行 1 はスロット 1,2 が列 1,2、スロット 3 はパディングで 0
    call monolis_test_check_eq_I1("monolis_spmat_convert_ell_test col i1 s1", mat%MAT%ELL%col((1-1)*N + 1), 1)
    call monolis_test_check_eq_I1("monolis_spmat_convert_ell_test col i1 s2", mat%MAT%ELL%col((2-1)*N + 1), 2)
    call monolis_test_check_eq_I1("monolis_spmat_convert_ell_test col i1 s3", mat%MAT%ELL%col((3-1)*N + 1), 0)

    !# 列番号：行 2 はスロット 1,2,3 が列 1,2,3
    call monolis_test_check_eq_I1("monolis_spmat_convert_ell_test col i2 s1", mat%MAT%ELL%col((1-1)*N + 2), 1)
    call monolis_test_check_eq_I1("monolis_spmat_convert_ell_test col i2 s2", mat%MAT%ELL%col((2-1)*N + 2), 2)
    call monolis_test_check_eq_I1("monolis_spmat_convert_ell_test col i2 s3", mat%MAT%ELL%col((3-1)*N + 2), 3)

    !# 値：行 1 スロット 1（主対角）は 2.0
    v = mat%MAT%R%Aell((1-1)*N + 1)
    call monolis_test_check_eq_R1("monolis_spmat_convert_ell_test val i1 s1", v, 2.0d0)
    !# 値：行 1 スロット 2（上三角）は 1.0
    v = mat%MAT%R%Aell((2-1)*N + 1)
    call monolis_test_check_eq_R1("monolis_spmat_convert_ell_test val i1 s2", v, 1.0d0)
    !# 値：行 1 スロット 3（パディング）は 0.0
    v = mat%MAT%R%Aell((3-1)*N + 1)
    call monolis_test_check_eq_R1("monolis_spmat_convert_ell_test val i1 s3 pad", v, 0.0d0)

    !# 値：行 2 スロット 1（下三角）は 1.0、スロット 2（主対角）は 2.0
    v = mat%MAT%R%Aell((1-1)*N + 2)
    call monolis_test_check_eq_R1("monolis_spmat_convert_ell_test val i2 s1", v, 1.0d0)
    v = mat%MAT%R%Aell((2-1)*N + 2)
    call monolis_test_check_eq_R1("monolis_spmat_convert_ell_test val i2 s2", v, 2.0d0)

    call monolis_dealloc_ELL_R(mat%MAT)

    call monolis_test_check_eq_I1("monolis_spmat_convert_ell_test Nmaxcol after dealloc", mat%MAT%ELL%Nmaxcol, 0)

    call monolis_finalize(mat)
  end subroutine monolis_spmat_convert_ell_test

end module mod_monolis_spmat_convert_ell_test
