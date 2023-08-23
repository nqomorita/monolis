!> 収束判定モジュール
module mod_monolis_converge
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_inner_product
  implicit none

contains

  !> @ingroup dev_linalg
  !> 収束判定閾値の設定
  subroutine monolis_set_converge_R(monoCOM, monoMAT, B, B2, is_converge, tdotp, tcomm)
    implicit none
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in,out] 行列構造体
    type(monolis_mat), intent(inout) :: monoMAT
    !> [in] 右辺ベクトル
    real(kdouble), intent(in) :: B(:)
    !> [out] L2 相対誤差の分母値
    real(kdouble), intent(out) :: B2
    !> [in,out] 収束判定フラグ
    logical, intent(inout) :: is_converge
    !> [in,out] 内積時間
    real(kdouble), intent(inout) :: tdotp
    !> [in,out] 通信時間
    real(kdouble), intent(inout) :: tcomm

    call monolis_std_debug_log_header("monolis_set_converge_R")

    is_converge = .false.
    call monolis_inner_product_main_R(monoCOM, monoMAT%N, monoMAT%NDOF, B, B, B2, tdotp, tcomm)

    if(B2 == 0.0d0)then
      if(monoCOM%my_rank == 0)then
        call monolis_std_error_string("norm of RHS vector B is 0.0")
      endif
      monoMAT%R%X = 0.0d0
      is_converge = .true.
    endif
  end subroutine monolis_set_converge_R

  !> @ingroup dev_linalg
  !> 収束判定閾値の設定
  subroutine monolis_set_converge_C(monoCOM, monoMAT, B, B2, is_converge, tdotp, tcomm)
    implicit none
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in,out] 行列構造体
    type(monolis_mat), intent(inout) :: monoMAT
    !> [in] 右辺ベクトル
    complex(kdouble), intent(in) :: B(:)
    !> [out] L2 相対誤差の分母値
    complex(kdouble), intent(out) :: B2
    !> [in,out] 収束判定フラグ
    logical, intent(inout) :: is_converge
    !> [in,out] 内積時間
    real(kdouble), intent(inout) :: tdotp
    !> [in,out] 通信時間
    real(kdouble), intent(inout) :: tcomm

    call monolis_std_debug_log_header("monolis_set_converge_C")

    is_converge = .false.
    call monolis_inner_product_main_C(monoCOM, monoMAT%N, monoMAT%NDOF, B, B, B2, tdotp, tcomm)

    if(B2 == 0.0d0)then
      if(monoCOM%my_rank == 0)then
        call monolis_std_error_string("norm of RHS vector B is 0.0")
      endif
      monoMAT%C%X = 0.0d0
      is_converge = .true.
    endif
  end subroutine monolis_set_converge_C

  !> @ingroup dev_linalg
  !> 収束の判定
  subroutine monolis_check_converge_R(monoPRM, monoCOM, monoMAT, R, B2, iter, is_converge, tdotp, tcomm)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: monoMAT
    !> [in] 残差ベクトル
    real(kdouble), intent(in) :: R(:)
    !> [in] L2 相対誤差の分母値
    real(kdouble), intent(in) :: B2
    !> [in] 反復回数
    integer(kint), intent(in) :: iter
    !> [in,out] 収束判定フラグ
    logical, intent(inout) :: is_converge
    !> [in,out] 内積時間
    real(kdouble), intent(inout) :: tdotp
    !> [in,out] 通信時間
    real(kdouble), intent(inout) :: tcomm
    real(kdouble) :: R2, resid

    call monolis_std_debug_log_header("monolis_check_converge_R")

    is_converge = .false.
    call monolis_inner_product_main_R(monoCOM, monoMAT%N, monoMAT%NDOF, R, R, R2, tdotp, tcomm)
    resid = dsqrt(R2/B2)

    monoPRM%Iarray(monolis_prm_I_cur_iter) = iter
    monoPRM%Rarray(monolis_prm_R_cur_resid) = resid

    if(monoCOM%my_rank == 0 .and. monoPRM%Iarray(monolis_prm_I_show_iterlog) == monolis_I_true)then
      write (*,"(i7, 1pe16.6)") iter, resid
    endif

    if(resid < monoPRM%Rarray(monolis_prm_R_tol))then
      is_converge = .true.
    endif

    if(iter == monoPRM%Iarray(monolis_prm_I_max_iter))then
      !is_converge = .true.
      if(monoPRM%Iarray(monolis_prm_I_is_measurement) == monolis_I_false)then
        call monolis_std_error_string("reached the maximum number of iterations")
        call monolis_std_error_stop()
      endif
    endif
  end subroutine monolis_check_converge_R

  !> @ingroup dev_linalg
  !> 収束の判定
  subroutine monolis_check_converge_C(monoPRM, monoCOM, monoMAT, R, B2, iter, is_converge, tdotp, tcomm)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: monoMAT
    !> [in] 残差ベクトル
    complex(kdouble), intent(in) :: R(:)
    !> [in] L2 相対誤差の分母値
    complex(kdouble), intent(in) :: B2
    !> [in] 反復回数
    integer(kint), intent(in) :: iter
    !> [in,out] 収束判定フラグ
    logical, intent(inout) :: is_converge
    !> [in,out] 内積時間
    real(kdouble), intent(inout) :: tdotp
    !> [in,out] 通信時間
    real(kdouble), intent(inout) :: tcomm
    complex(kdouble) :: R2, resid

    call monolis_std_debug_log_header("monolis_check_converge_C")

    is_converge = .false.
    call monolis_inner_product_main_C(monoCOM, monoMAT%N, monoMAT%NDOF, R, R, R2, tdotp, tcomm)
    resid = sqrt(R2/B2)

    monoPRM%Iarray(monolis_prm_I_cur_iter) = iter
    monoPRM%Rarray(monolis_prm_R_cur_resid) = resid

    if(monoCOM%my_rank == 0 .and. monoPRM%Iarray(monolis_prm_I_show_iterlog) == monolis_I_true)then
      write (*,"(i7, 1pe16.6)") iter, abs(resid)
    endif

    if(abs(resid) < monoPRM%Rarray(monolis_prm_R_tol))then
      is_converge = .true.
    endif

    if(iter == monoPRM%Iarray(monolis_prm_I_max_iter))then
      !is_converge = .true.
      if(monoPRM%Iarray(monolis_prm_I_is_measurement) == monolis_I_false)then
        call monolis_std_error_string("reached the maximum number of iterations")
        call monolis_std_error_stop()
      endif
    endif
  end subroutine monolis_check_converge_C
end module mod_monolis_converge