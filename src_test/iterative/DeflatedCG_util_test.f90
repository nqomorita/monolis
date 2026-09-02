!> DeflatedCG 法テストモジュール
module mod_monolis_solver_DeflatedCG_util_test
  use mod_monolis
  use mod_monolis_solver_DeflatedCG_util

  implicit none

contains

  subroutine monolis_solver_DeflatedCG_util_test()
    implicit none
    write(*,*)"monolis_solver_DeflatedCG_util_test"

    call deflatedCG_E_initialize_tolerance_test()
    call deflatedCG_residual_replacement_test()
  end subroutine monolis_solver_DeflatedCG_util_test

  subroutine deflatedCG_E_initialize_tolerance_test()
    implicit none
    type(monolis_structure) :: mat
    type(monolis_prm) :: prm_coarse
    type(monolis_com) :: com, com_coarse
    type(monolis_mat) :: mat_coarse
    integer(kint) :: elem(2,1)
    real(kdouble) :: W(2,1)
    real(kdouble), allocatable :: AW(:,:), WtA(:,:)

    call monolis_std_global_log_string("deflatedCG_E_initialize")

    call monolis_initialize(mat)
    call monolis_com_initialize_by_self(com)

    elem(:,1) = (/1, 2/)
    call monolis_get_nonzero_pattern_by_simple_mesh_R(mat, 2, 2, 1, 1, elem)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 1, 1, 1, 1, 2.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 2, 2, 1, 1, 3.0d0)

    call monolis_set_tolerance(mat, 3.0d-7)
    mat%PRM%Iarray(monolis_prm_I_DCG_inner_method) = monolis_iter_CG
    mat%PRM%Iarray(monolis_prm_I_DCG_inner_prec) = monolis_prec_NONE
    mat%PRM%Iarray(monolis_prm_I_DCG_inner_max_iter) = 10

    W(:,1) = (/1.0d0, 0.0d0/)
    call deflatedCG_E_initialize(com, mat%PRM, mat%MAT, prm_coarse, com_coarse, mat_coarse, &
      & 1, 1, 2, W, AW, WtA)

    call monolis_test_check_eq_R1("deflatedCG_E_initialize relative tolerance", &
      & prm_coarse%Rarray(monolis_prm_R_tol), 3.0d-9)

    call monolis_dealloc_R_2d(AW)
    call monolis_dealloc_R_2d(WtA)
    call monolis_prm_finalize(prm_coarse)
    call monolis_mat_finalize(mat_coarse)
    call monolis_com_finalize(com_coarse)
    call monolis_com_finalize(com)
    call monolis_finalize(mat)
  end subroutine deflatedCG_E_initialize_tolerance_test

  subroutine deflatedCG_residual_replacement_test()
    implicit none
    type(monolis_structure) :: mat
    type(monolis_com) :: com
    integer(kint) :: elem(2,1)
    integer(kint), allocatable :: IPV_R(:)
    real(kdouble) :: X(2), B(2), R(2), R_ans(2), W(2,1)
    real(kdouble) :: tspmv, tcomm_spmv
    real(kdouble), allocatable :: WtW(:,:)

    call monolis_std_global_log_string("deflatedCG_residual_replacement")

    call monolis_initialize(mat)
    call monolis_com_initialize_by_self(com)

    elem(:,1) = (/1, 2/)
    call monolis_get_nonzero_pattern_by_simple_mesh_R(mat, 2, 2, 1, 1, elem)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 1, 1, 1, 1, 2.0d0)
    call monolis_add_scalar_to_sparse_matrix_R(mat, 2, 2, 1, 1, 3.0d0)

    X = (/1.0d0, 2.0d0/)
    B = (/5.0d0, 11.0d0/)
    R = 100.0d0
    W(:,1) = (/1.0d0, 0.0d0/)
    R_ans = (/0.0d0, 5.0d0/)
    tspmv = 0.0d0
    tcomm_spmv = 0.0d0

    call deflatedCG_residual_replacement_initialize(1, 2, W, WtW, IPV_R)
    call deflatedCG_residual_replacement( &
      & com, mat%MAT, X, B, 1, 2, 1, W, R, WtW, IPV_R, tspmv, tcomm_spmv)

    call monolis_test_check_eq_R("deflatedCG_residual_replacement true residual", R, R_ans)

    call monolis_dealloc_R_2d(WtW)
    call monolis_dealloc_I_1d(IPV_R)
    call monolis_com_finalize(com)
    call monolis_finalize(mat)
  end subroutine deflatedCG_residual_replacement_test

end module mod_monolis_solver_DeflatedCG_util_test
