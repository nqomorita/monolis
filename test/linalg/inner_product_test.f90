!> ベクトル内積テストモジュール
module mod_monolis_linalg_test
  use mod_monolis
  implicit none

contains

  subroutine monolis_linalg_test
    implicit none

    !call monolis_inner_product_I(monolis, n, ndof, X, Y, sum)
    !call monolis_inner_product_main_I(monoCOM, n, ndof, X, Y, sum, tdotp, tcomm)
    !call monolis_inner_product_R(monolis, n, ndof, X, Y, sum)
    !call monolis_inner_product_main_R(monoCOM, n, ndof, X, Y, sum, tdotp, tcomm)
    !call monolis_inner_product_C(monolis, n, ndof, X, Y, sum)
    !call monolis_inner_product_main_C(monoCOM, n, ndof, X, Y, sum, tdotp, tcomm)
    !call monolis_inner_product_main_R_no_comm(monoCOM, n, ndof, X, Y, sum)
  end subroutine monolis_linalg_test
end module mod_monolis_linalg_test
