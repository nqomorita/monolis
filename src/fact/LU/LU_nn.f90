!> LU 分解（nxn ブロック）
module mod_monolis_fact_LU_nn
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_fact_fillin
  use mod_monolis_fact_analysis
  use mod_monolis_fact_factorize
  use mod_monolis_spmat_reorder

  implicit none

contains

  !> @ingroup prec
  !> 前処理生成：LU 前処理（nxn ブロック、実数型）
  subroutine monolis_fact_LU_nn_setup_R(monoMAT, monoLU)
    implicit none
    !> [in] 行列構造体
    type(monolis_mat), target, intent(in) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), target, intent(inout) :: monoLU
    type(monolis_mat) :: monoMAT_reorder
    real(kdouble) :: t(20)

    t(1) = monolis_get_time()

    !> analysis phase
    call monolis_matrix_reordering_fw_R(monoMAT, monoMAT_reorder)

    t(2) = monolis_get_time()
    write(*,"(a,1pe10.3)")"monolis_matrix_reordering_fw_R             ", t(2) - t(1)

    !monoMAT%REORDER%perm
    !monoMAT%REORDER%iperm

  end subroutine monolis_fact_LU_nn_setup_R

  !> @ingroup prec
  !> 前処理生成：LU 前処理（nxn ブロック、複素数型）
  subroutine monolis_fact_LU_nn_setup_C(monoMAT, monoLU)
    implicit none
    !> [in] 行列構造体
    type(monolis_mat), target, intent(in) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), target, intent(inout) :: monoLU

  end subroutine monolis_fact_LU_nn_setup_C

  !> @ingroup prec
  !> 前処理適用：LU 前処理（nxn ブロック、実数型）
  subroutine monolis_fact_LU_nn_apply_R(monoMAT, monoLU, Y, X)
    implicit none
    !> [in] 行列構造体
    type(monolis_mat), target, intent(in) :: monoMAT
    !> [in] 前処理構造体
    type(monolis_mat), target, intent(in) :: monoLU
    real(kdouble) :: X(:), Y(:)

  end subroutine monolis_fact_LU_nn_apply_R

  !> 前処理適用：LU 前処理（nxn ブロック、複素数型）
  subroutine monolis_fact_LU_nn_apply_C(monoMAT, monoLU, X, Y)
    implicit none
    !> [in] 行列構造体
    type(monolis_mat), target, intent(in) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), target, intent(inout) :: monoLU
    complex(kdouble) :: X(:), Y(:)

  end subroutine monolis_fact_LU_nn_apply_C

  !> @ingroup prec
  !> 前処理初期化：LU 前処理（nxn ブロック、実数型）
  subroutine monolis_fact_LU_nn_clear_R(monoLU)
    implicit none
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoLU

    call monolis_pdealloc_R_1d(monoLU%R%D)
  end subroutine monolis_fact_LU_nn_clear_R

  !> @ingroup prec
  !> 前処理初期化：LU 前処理（nxn ブロック、複素数型）
  subroutine monolis_fact_LU_nn_clear_C(monoLU)
    implicit none
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoLU

    call monolis_pdealloc_C_1d(monoLU%C%D)
  end subroutine monolis_fact_LU_nn_clear_C
end module mod_monolis_fact_LU_nn
