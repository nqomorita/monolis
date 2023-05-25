!> ScaLAPACK ラッパーモジュール
module mod_monolis_scalapack_test
  use mod_monolis_utils
  use mod_monolis_lapack

  implicit none

contains

  subroutine monolis_scalapack_test()
    implicit none
    !> 行列の大きさ
    integer(kint) :: n
    !> 行列の対角成分
    real(kdouble) :: D(5)

    call monolis_std_global_log_string("monolis_scalapack_test")

  end subroutine monolis_scalapack_test

end module mod_monolis_scalapack_test
