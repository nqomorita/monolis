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
  subroutine monolis_inner_product_I_c(m, x, y, sum, comm)&
    & bind(c, name = "monolis_inner_product_I_c_main")
    implicit none
    !> [in] 内部計算点数✕計算点が持つ自由度
    integer(kint_c), intent(in), value :: m
    !> [in] ベクトル 1
    integer(kint_c), intent(in) :: x(m)
    !> [in] ベクトル 2
    integer(kint_c), intent(in) :: y(m)
    !> [out] 内積結果
    integer(kint_c) :: sum
    !> [in] MPI コミュニケータ
    integer(kint_c), intent(in), value :: comm
    type(monolis_com) :: monoCOM
    real(kdouble) :: tdotp, tcomm

    monoCOM%comm = comm
    call monolis_inner_product_main_I(monoCOM, m, X, Y, sum, tdotp, tcomm)
  end subroutine monolis_inner_product_I_c

  !> @ingroup dev_linalg
  !> ベクトル内積（実数型）
  subroutine monolis_inner_product_R_c(m, x, y, sum, comm)&
    & bind(c, name = "monolis_inner_product_R_c_main")
    implicit none
    !> [in] 内部計算点数✕計算点が持つ自由度
    integer(kint_c), intent(in), value :: m
    !> [in] ベクトル 1
    real(c_double), intent(in) :: x(m)
    !> [in] ベクトル 2
    real(c_double), intent(in) :: y(m)
    !> [out] 内積結果
    real(c_double) :: sum
    !> [in] MPI コミュニケータ
    integer(kint_c), intent(in), value :: comm
    type(monolis_com) :: monoCOM
    real(kdouble) :: tdotp, tcomm

    monoCOM%comm = comm
    call monolis_inner_product_main_R(monoCOM, m, X, Y, sum, tdotp, tcomm)
  end subroutine monolis_inner_product_R_c

  !> @ingroup dev_linalg
  !> ベクトル内積（複素数型）
  subroutine monolis_inner_product_C_c(m, x, y, sum, comm)&
    & bind(c, name = "monolis_inner_product_C_c_main")
    implicit none
    !> [in] 内部計算点数✕計算点が持つ自由度
    integer(kint_c), intent(in), value :: m
    !> [in] ベクトル 1
    complex(c_double), intent(in) :: x(m)
    !> [in] ベクトル 2
    complex(c_double), intent(in) :: y(m)
    !> [out] 内積結果
    complex(c_double) :: sum
    !> [in] MPI コミュニケータ
    integer(kint_c), intent(in), value :: comm
    type(monolis_com) :: monoCOM
    real(kdouble) :: tdotp, tcomm

    monoCOM%comm = comm
    call monolis_inner_product_main_C(monoCOM, m, X, Y, sum, tdotp, tcomm)
  end subroutine monolis_inner_product_C_c
end module mod_monolis_inner_product_wrap