!> 疎行列ベクトル積テストモジュール
module mod_monolis_matvec_test
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_matvec
  implicit none

contains

  subroutine monolis_matvec_test()
    implicit none

    !call monolis_matvec_nn_R(monoCOM, monoMAT, X, Y, NDOF)
    !call monolis_matvec_11_R(monoCOM, monoMAT, X, Y)
    !call monolis_matvec_33_R(monoCOM, monoMAT, X, Y)
    !call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm)
  end subroutine monolis_matvec_test

end module mod_monolis_matvec_test
