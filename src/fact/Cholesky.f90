!> Cholesky 分解（nxn ブロック）
module mod_monolis_fact_Cholesky_nn
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_fact_analysis
  use mod_monolis_fact_factorize_cholesky
  use mod_monolis_fact_solve_cholesky

  implicit none

contains

  !> @ingroup prec
  !> 前処理生成：Cholesky 前処理（nxn ブロック、実数型）
  subroutine monolis_fact_Cholesky_nn_setup_R(monoMAT, monoCHOL)
    implicit none
    !> [in] 行列構造体
    type(monolis_mat), target, intent(in) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), target, intent(inout) :: monoCHOL

    call monolis_fact_analysis(monoMAT, monoCHOL%LU)
    call monolis_fact_factorize_cholesky(monoMAT, monoCHOL%LU)
  end subroutine monolis_fact_Cholesky_nn_setup_R

  !> @ingroup prec
  !> 前処理生成：Cholesky 前処理（nxn ブロック、複素数型）
  subroutine monolis_fact_Cholesky_nn_setup_C(monoMAT, monoCHOL)
    implicit none
    !> [in] 行列構造体
    type(monolis_mat), target, intent(in) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), target, intent(inout) :: monoCHOL

  end subroutine monolis_fact_Cholesky_nn_setup_C

  !> @ingroup prec
  !> 前処理適用：Cholesky 前処理（nxn ブロック、実数型）
  subroutine monolis_fact_Cholesky_nn_apply_R(monoMAT, monoCHOL, Y, X)
    implicit none
    !> [in] 行列構造体
    type(monolis_mat), target, intent(in) :: monoMAT
    !> [in] 前処理構造体
    type(monolis_mat), target, intent(in) :: monoCHOL
    real(kdouble) :: X(:), Y(:)

    integer(kint) :: n

    n = monoCHOL%LU%N
    if (n <= 0) return
    X(1:n) = Y(1:n)
    call monolis_fact_solve_cholesky(monoCHOL%LU, X)
  end subroutine monolis_fact_Cholesky_nn_apply_R

  !> 前処理適用：Cholesky 前処理（nxn ブロック、複素数型）
  subroutine monolis_fact_Cholesky_nn_apply_C(monoMAT, monoCHOL, X, Y)
    implicit none
    !> [in] 行列構造体
    type(monolis_mat), target, intent(in) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), target, intent(inout) :: monoCHOL
    complex(kdouble) :: X(:), Y(:)

  end subroutine monolis_fact_Cholesky_nn_apply_C

  !> @ingroup prec
  !> 前処理初期化：Cholesky 前処理（nxn ブロック、実数型）
  subroutine monolis_fact_Cholesky_nn_clear_R(monoCHOL)
    implicit none
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoCHOL

    call monolis_mat_finalize_LU(monoCHOL%LU)
    call monolis_pdealloc_R_1d(monoCHOL%R%D)
  end subroutine monolis_fact_Cholesky_nn_clear_R

  !> @ingroup prec
  !> 前処理初期化：Cholesky 前処理（nxn ブロック、複素数型）
  subroutine monolis_fact_Cholesky_nn_clear_C(monoCHOL)
    implicit none
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoCHOL

    call monolis_pdealloc_C_1d(monoCHOL%C%D)
  end subroutine monolis_fact_Cholesky_nn_clear_C
end module mod_monolis_fact_Cholesky_nn
