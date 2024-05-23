!> LU 分解（nxn ブロック）
module mod_monolis_fact_LU_nn
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_fact_fillin
  use mod_monolis_fact_analysis

  implicit none

  integer(kint) :: n_fact_array
  integer(kint), allocatable :: fact_order(:)
  integer(kint), allocatable :: fact_array_index(:)
  integer(kint), allocatable :: add_location(:)
  real(kdouble), allocatable :: fact_array(:)

contains

  !> @ingroup prec
  !> 前処理生成：LU 前処理（nxn ブロック、実数型）
  subroutine monolis_fact_LU_nn_setup_R(monoMAT, monoPREC)
    implicit none
    !> [in] 行列構造体
    type(monolis_mat), target, intent(in) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), target, intent(inout) :: monoPREC
    logical :: is_asym = .false.
    logical :: is_fillin = .true.

    call monolis_matrix_get_fillin(monoMAT, monoPREC, is_asym, is_fillin)

    !write(*,*)"monoPREC%N", monoPREC%N
    !write(*,*)"monoPREC%NDOF", monoPREC%NDOF
    !write(*,*)"monoPREC%SCSR%indexU", monoPREC%SCSR%indexU
    !write(*,*)"monoPREC%SCSR%itemU", monoPREC%SCSR%itemU

    call monolis_matrix_alloc_with_fillin(monoPREC, is_asym)

    call monolis_matrix_get_factorize_order(monoPREC, fact_order)
    !write(*,*)"fact_order", fact_order

    call monolis_matrix_get_factorize_array(monoPREC, fact_order, n_fact_array, fact_array, fact_array_index)
    !write(*,*)"n_fact_array", n_fact_array
    !write(*,*)"fact_array", fact_array
    !write(*,*)"fact_array_index", fact_array_index

    call monolis_matrix_set_value_of_factorize_array(monoMAT, monoPREC, fact_order, n_fact_array, fact_array, fact_array_index)
    !write(*,*)"fact_array", fact_array

    call monolis_matrix_get_add_location(monoPREC, fact_order, fact_array_index, add_location)
    !write(*,*)"add_location", add_location
  end subroutine monolis_fact_LU_nn_setup_R

  !> @ingroup prec
  !> 前処理生成：LU 前処理（nxn ブロック、複素数型）
  subroutine monolis_fact_LU_nn_setup_C(monoMAT, monoPREC)
    implicit none
    !> [in] 行列構造体
    type(monolis_mat), target, intent(in) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), target, intent(inout) :: monoPREC

  end subroutine monolis_fact_LU_nn_setup_C

  !> @ingroup prec
  !> 前処理適用：LU 前処理（3x3 ブロック、実数型）
  subroutine monolis_fact_LU_nn_apply_R(monoMAT, monoPREC, X, Y)
    implicit none
    !> [in] 行列構造体
    type(monolis_mat), target, intent(in) :: monoMAT
    !> [in] 前処理構造体
    type(monolis_mat), target, intent(in) :: monoPREC
    real(kdouble) :: X(:), Y(:)

    Y = X
  end subroutine monolis_fact_LU_nn_apply_R

  !> 前処理適用：LU 前処理（3x3 ブロック、複素数型）
  subroutine monolis_fact_LU_nn_apply_C(monoMAT, monoPREC, X, Y)
    implicit none
    !> [in] 行列構造体
    type(monolis_mat), target, intent(in) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), target, intent(inout) :: monoPREC
    complex(kdouble) :: X(:), Y(:)

  end subroutine monolis_fact_LU_nn_apply_C

  !> @ingroup prec
  !> 前処理初期化：LU 前処理（nxn ブロック、実数型）
  subroutine monolis_fact_LU_nn_clear_R(monoPREC)
    implicit none
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC

    call monolis_pdealloc_R_1d(monoPREC%R%D)
  end subroutine monolis_fact_LU_nn_clear_R

  !> @ingroup prec
  !> 前処理初期化：LU 前処理（nxn ブロック、複素数型）
  subroutine monolis_fact_LU_nn_clear_C(monoPREC)
    implicit none
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC

    call monolis_pdealloc_C_1d(monoPREC%C%D)
  end subroutine monolis_fact_LU_nn_clear_C
end module mod_monolis_fact_LU_nn
