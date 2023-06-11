!> ベクトル演算テストモジュール
module mod_monolis_vec_util_test
  use mod_monolis
  use mod_monolis_utils
  use mod_monolis_vec_util
  implicit none

contains

  subroutine monolis_vec_util_test()
    implicit none

    call monolis_vec_copy_test()
    call monolis_vec_AXPBY_test()
  end subroutine monolis_vec_util_test

  subroutine monolis_vec_copy_test()
    implicit none
    integer(kint) :: n
    integer(kint) :: ndof
    integer(kint) :: iX(4), iY(4)
    real(kdouble) :: rX(4), rY(4)
    complex(kdouble) :: cX(4), cY(4)

    call monolis_std_global_log_string("monolis_vec_copy_I")
    call monolis_std_global_log_string("monolis_vec_copy_R")
    call monolis_std_global_log_string("monolis_vec_copy_C")

    n = 2

    ndof = 2

    iX = 1

    iY = 0

    call monolis_vec_copy_I(n, ndof, iX, iY)

    call monolis_test_check_eq_I1("monolis_vec_copy_test 1a", iY(1), 1)
    call monolis_test_check_eq_I1("monolis_vec_copy_test 1b", iY(2), 1)
    call monolis_test_check_eq_I1("monolis_vec_copy_test 1c", iY(3), 1)
    call monolis_test_check_eq_I1("monolis_vec_copy_test 1d", iY(4), 1)

    rX = 1.0d0

    rY = 0.0d0

    call monolis_vec_copy_R(n, ndof, rX, rY)

    call monolis_test_check_eq_R1("monolis_vec_copy_test 2a", rY(1), 1.0d0)
    call monolis_test_check_eq_R1("monolis_vec_copy_test 2b", rY(2), 1.0d0)
    call monolis_test_check_eq_R1("monolis_vec_copy_test 2c", rY(3), 1.0d0)
    call monolis_test_check_eq_R1("monolis_vec_copy_test 2d", rY(4), 1.0d0)

    cX = (1.0d0, 1.0d0)

    cY = (0.0d0, 0.0d0)

    call monolis_vec_copy_C(n, ndof, cX, cY)

    call monolis_test_check_eq_C1("monolis_vec_copy_test 3a", cY(1), (1.0d0, 1.0d0))
    call monolis_test_check_eq_C1("monolis_vec_copy_test 3b", cY(2), (1.0d0, 1.0d0))
    call monolis_test_check_eq_C1("monolis_vec_copy_test 3c", cY(3), (1.0d0, 1.0d0))
    call monolis_test_check_eq_C1("monolis_vec_copy_test 3d", cY(4), (1.0d0, 1.0d0))
  end subroutine monolis_vec_copy_test

  subroutine monolis_vec_AXPBY_test()
    implicit none
    !> 自由度数
    integer(kint) :: n
    !> ブロックサイズ
    integer(kint) :: ndof
    !> ベクトル 1
    integer(kint) :: iX(4), iY(4), iZ(4), ia, ib
    real(kdouble) :: rX(4), rY(4), rZ(4), ra, rb
    complex(kdouble) :: cX(4), cY(4), cZ(4), ca, cb

    call monolis_std_global_log_string("monolis_vec_AXPBY_I")
    call monolis_std_global_log_string("monolis_vec_AXPBY_R")
    call monolis_std_global_log_string("monolis_vec_AXPBY_C")

    n = 2

    ndof = 2

    ia = 2

    ib = 2

    iX = 1

    iY = 2

    call monolis_vec_AXPBY_I(n, ndof, ia, iX, ib, iY, iZ)

    call monolis_test_check_eq_I1("monolis_vec_AXPBY_test 1a", iZ(1), 6)
    call monolis_test_check_eq_I1("monolis_vec_AXPBY_test 1b", iZ(2), 6)
    call monolis_test_check_eq_I1("monolis_vec_AXPBY_test 1c", iZ(3), 6)
    call monolis_test_check_eq_I1("monolis_vec_AXPBY_test 1d", iZ(4), 6)

    ra = 2.0d0

    rb = 2.0d0

    rX = 1.0d0

    rY = 2.0d0

    call monolis_vec_AXPBY_R(n, ndof, ra, rX, rb, rY, rZ)

    call monolis_test_check_eq_R1("monolis_vec_AXPBY_test 2a", rZ(1), 6.0d0)
    call monolis_test_check_eq_R1("monolis_vec_AXPBY_test 2b", rZ(2), 6.0d0)
    call monolis_test_check_eq_R1("monolis_vec_AXPBY_test 2c", rZ(3), 6.0d0)
    call monolis_test_check_eq_R1("monolis_vec_AXPBY_test 2d", rZ(4), 6.0d0)

    ca = (2.0d0, 0.0d0)

    cb = (2.0d0, 0.0d0)

    cX = (1.0d0, 1.0d0)

    cY = (2.0d0, 2.0d0)

    call monolis_vec_AXPBY_C(n, ndof, ca, cX, cb, cY, cZ)

    call monolis_test_check_eq_C1("monolis_vec_AXPBY_test 3a", cZ(1), (6.0d0, 6.0d0))
    call monolis_test_check_eq_C1("monolis_vec_AXPBY_test 3b", cZ(2), (6.0d0, 6.0d0))
    call monolis_test_check_eq_C1("monolis_vec_AXPBY_test 3c", cZ(3), (6.0d0, 6.0d0))
    call monolis_test_check_eq_C1("monolis_vec_AXPBY_test 3d", cZ(4), (6.0d0, 6.0d0))
  end subroutine monolis_vec_AXPBY_test

end module mod_monolis_vec_util_test
