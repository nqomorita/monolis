!> Adam オプティマイザモジュールのテスト
module mod_monolis_opt_adam_test
  use mod_monolis
  implicit none

contains

  !> @ingroup optimize
  !> Adam 関連サブルーチンの一括テスト
  subroutine monolis_optimize_adam_test()
    implicit none
    type(monolis_opt_adam_state) :: st
    real(kdouble) :: W(2, 3), b(3), gW(2, 3), gb(3)
    integer(kint) :: t

    !> cover_check 用に対象サブルーチン名を全て登録
    call monolis_std_global_log_string("monolis_opt_adam_init")
    call monolis_std_global_log_string("monolis_opt_adam_free")
    call monolis_std_global_log_string("monolis_opt_adam_apply")

    W  = 1.0d0
    b  = 0.0d0
    gW = 0.5d0
    gb = -0.5d0

    call monolis_opt_adam_init(st, 2, 3)
    call monolis_test_check_eq_L1("adam_test m allocated", allocated(st%m), .true.)
    call monolis_test_check_eq_L1("adam_test vb allocated", allocated(st%vb), .true.)

    do t = 1, 5
      call monolis_opt_adam_apply(W, b, gW, gb, st, t, 1.0d-2)
    end do

    !> 勾配 gW > 0 のため W は減少、gb < 0 のため b は増加する
    call monolis_test_check_eq_L1("adam_test W finite", &
      all(W == W) .and. all(abs(W) < huge(0.0d0)), .true.)
    call monolis_test_check_eq_L1("adam_test W decreased", all(W < 1.0d0), .true.)
    call monolis_test_check_eq_L1("adam_test b increased", all(b > 0.0d0), .true.)

    call monolis_opt_adam_free(st)
    call monolis_test_check_eq_L1("adam_test m freed", allocated(st%m), .false.)
  end subroutine monolis_optimize_adam_test

end module mod_monolis_opt_adam_test
