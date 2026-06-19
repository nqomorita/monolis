!> 機械学習用パラメータ定義モジュールのテスト
module mod_monolis_def_opt_test
  use mod_monolis
  implicit none

contains

  !> @ingroup optimize
  !> mod_monolis_def_opt の定義値テスト (kdouble_ml が 32bit float であること)
  subroutine monolis_def_opt_test()
    implicit none
    real(kdouble_ml) :: x

    call monolis_test_check_eq_I1("def_opt kdouble_ml kind", kind(x), 4)
    call monolis_test_check_eq_I1("def_opt kdouble_ml bit size", storage_size(x), 32)
  end subroutine monolis_def_opt_test

end module mod_monolis_def_opt_test
