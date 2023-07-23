!> ベクトル内積関数群
module mod_monolis_inner_product_wrap
  use mod_monolis_utils
  use mod_monolis_def_struc
  use mod_monolis_inner_product
  use iso_c_binding
  implicit none

contains

  !> @ingroup dev_linalg
  !> ベクトル内積（整数型）
  subroutine monolis_inner_product_I_c(n, n_dof, x, y, sum, comm)&
    & bind(c, name = "monolis_inner_product_I_c_main")
    implicit none
    !> [in] 計算点数
    integer(c_int), intent(in), value :: n
    !> [in] 計算点あたりの自由度
    integer(c_int), intent(in), value :: n_dof
    !> [in] ベクトル 1
    integer(c_int), intent(in) :: x(n*n_dof)
    !> [in] ベクトル 2
    integer(c_int), intent(in) :: y(n*n_dof)
    !> [out] 内積結果
    integer(c_int) :: sum
    !> [in] MPI コミュニケータ
    integer(c_int), intent(in), value :: comm
    type(monolis_com) :: monoCOM
    real(kdouble) :: tdotp, tcomm

    monoCOM%comm = comm
    call monolis_inner_product_main_I(monoCOM, n, n_dof, X, Y, sum, tdotp, tcomm)
  end subroutine monolis_inner_product_I_c

  !> @ingroup dev_linalg
  !> ベクトル内積（実数型）
  subroutine monolis_inner_product_R_c(n, n_dof, x, y, sum, comm)&
    & bind(c, name = "monolis_inner_product_R_c_main")
    implicit none
    !> [in] 計算点数
    integer(c_int), intent(in), value :: n
    !> [in] 計算点あたりの自由度
    integer(c_int), intent(in), value :: n_dof
    !> [in] ベクトル 1
    real(c_double), intent(in) :: x(n*n_dof)
    !> [in] ベクトル 2
    real(c_double), intent(in) :: y(n*n_dof)
    !> [out] 内積結果
    real(c_double) :: sum
    !> [in] MPI コミュニケータ
    integer(c_int), intent(in), value :: comm
    type(monolis_com) :: monoCOM
    real(kdouble) :: tdotp, tcomm

    monoCOM%comm = comm
    call monolis_inner_product_main_R(monoCOM, n, n_dof, X, Y, sum, tdotp, tcomm)
  end subroutine monolis_inner_product_R_c

  !> @ingroup dev_linalg
  !> ベクトル内積（複素数型）
  subroutine monolis_inner_product_C_c(n, n_dof, x, y, sum, comm)&
    & bind(c, name = "monolis_inner_product_C_c_main")
    implicit none
    !> [in] 計算点数
    integer(c_int), intent(in), value :: n
    !> [in] 計算点あたりの自由度
    integer(c_int), intent(in), value :: n_dof
    !> [in] ベクトル 1
    complex(c_double), intent(in) :: x(n*n_dof)
    !> [in] ベクトル 2
    complex(c_double), intent(in) :: y(n*n_dof)
    !> [out] 内積結果
    complex(c_double) :: sum
    !> [in] MPI コミュニケータ
    integer(c_int), intent(in), value :: comm
    type(monolis_com) :: monoCOM
    real(kdouble) :: tdotp, tcomm

    monoCOM%comm = comm
    call monolis_inner_product_main_C(monoCOM, n, n_dof, X, Y, sum, tdotp, tcomm)
  end subroutine monolis_inner_product_C_c
end module mod_monolis_inner_product_wrap