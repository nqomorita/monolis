!> 収束判定モジュール
module mod_monolis_converge
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  implicit none

contains

  !> @ingroup dev_linalg
  !> 収束判定閾値の設定
  subroutine monolis_set_converge_R(monoPRM, monoCOM, monoMAT, B, B2, is_converge, tdotp, tcomm)
    implicit none
    !> パラメータ構造体
    type(monolis_prm) :: monoPRM
    !> 通信テーブル構造体
    type(monolis_com) :: monoCOM
    !> 行列構造体
    type(monolis_mat) :: monoMAT
    !> 右辺ベクトル
    real(kdouble) :: B(:)
    !> L2 相対誤差の分母値
    real(kdouble) :: B2
    !> 収束判定フラグ
    logical :: is_converge
    !> 内積時間
    real(kdouble) :: tdotp
    !> 通信時間
    real(kdouble) :: tcomm

#ifdef DEBUG
    call monolis_std_debug_log_header("monolis_set_converge_R")
#endif

    is_converge = .false.
    call monolis_inner_product_R(monoCOM, monoMAT%N, monoMAT%NDOF, B, B, B2, tdotp, tcomm)

    if(B2 == 0.0d0)then
      if(monoCOM%my_rank == 0)then
        !write (*,"(a,1pe16.6)")" ** monolis warning: bnorm ", B2
      endif
      monoMAT%R%X = 0.0d0
      is_converge = .true.
    endif
  end subroutine monolis_set_converge_R

  !> @ingroup dev_linalg
  !> 収束の判定
  subroutine monolis_check_converge_R(monoPRM, monoCOM, monoMAT, R, B2, iter, is_converge, tdotp, tcomm)
    implicit none
    !> パラメータ構造体
    type(monolis_prm) :: monoPRM
    !> 通信テーブル構造体
    type(monolis_com) :: monoCOM
    !> 行列構造体
    type(monolis_mat) :: monoMAT
    !> 残差ベクトル
    real(kdouble) :: R(:)
    !> L2 相対誤差の分母値
    real(kdouble) :: B2
    !> 反復回数
    integer(kint) :: iter
    !> 収束判定フラグ
    logical :: is_converge
    !> 内積時間
    real(kdouble) :: tdotp
    !> 通信時間
    real(kdouble) :: tcomm
    real(kdouble) :: R2, resid

    is_converge = .false.
    call monolis_inner_product_R(monoCOM, monoMAT%N, monoMAT%NDOF, R, R, R2, tdotp, tcomm)
    resid = dsqrt(R2/B2)

    monoPRM%Iarray(monolis_prm_I_cur_iter) = iter
    monoPRM%Rarray(monolis_prm_R_cur_resid) = resid

    if(resid < monoPRM%Rarray(monolis_prm_R_tol))then
      is_converge = .true.
      if(monoCOM%my_rank == 0 .and. monoPRM%Iarray(monolis_prm_I_show_iterlog) == monolis_I_true)then
        !write (*,"(i7, 1pe16.6)") iter, resid
      endif
    endif

    if(iter == monoPRM%Iarray(monolis_prm_I_max_iter))then
      !is_converge = .true.
      if(monoPRM%Iarray(monolis_prm_I_is_measurement) == monolis_I_false)then
        !error stop "* monolis error: reached the maximum number of iterations"
      endif
    endif
  end subroutine monolis_check_converge_R
end module mod_monolis_converge