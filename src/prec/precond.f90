!> 前処理関連モジュール
module mod_monolis_precond
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_precond_diag
  use mod_monolis_precond_SOR

  implicit none

contains

  !> 前処理生成関数
  subroutine monolis_precond_setup_R(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> パラメータ構造体
    type(monolis_prm) :: monoPRM
    !> 通信テーブル構造体
    type(monolis_com) :: monoCOM
    !> 行列構造体
    type(monolis_mat) :: monoMAT
    !> 前処理構造体
    type(monolis_mat) :: monoPREC
    integer(kint) :: precond
    real(kdouble) :: t1, t2

#ifdef DEBUG
    call monolis_std_debug_log_header("monolis_precond_setup_R")
#endif

    t1 = monolis_get_time()

    precond = monoPRM%Iarray(monolis_prm_I_precond)

    if(precond == monolis_prec_DIAG)then
      call monolis_precond_diag_setup(monoPRM, monoCOM, monoMAT, monoPREC)
    !elseif(precond == monolis_prec_ILU)then
    !  call monolis_precond_ilu_setup(monoPRM, monoCOM, monoMAT)
    !elseif(precond == monolis_prec_JACOBI)then
    !  call monolis_precond_Jacobi_setup(monoPRM, monoCOM, monoMAT)
    elseif(precond == monolis_prec_SOR)then
      call monolis_precond_SOR_setup(monoPRM, monoCOM, monoMAT, monoPREC)
    !elseif(precond == monolis_prec_MUMPS)then
    !  call monolis_precond_MUMPS_setup(monoPRM, monoCOM, monoMAT)
    !elseif(precond == monolis_prec_MUMPS_LOCAL)then
    !  call monolis_precond_MUMPS_setup_local()
    !  call monolis_precond_MUMPS_setup(monoPRM, monoCOM, monoMAT)
    !elseif(precond == monolis_prec_MF)then
    !  call monolis_precond_MF_setup(monoPRM, monoCOM, monoMAT)
    endif

    t2 = monolis_get_time()
    !monoPRM%tprep = monoPRM%tprep + t2 - t1
  end subroutine monolis_precond_setup_R

  !> 前処理適用関数
  subroutine monolis_precond_apply_R(monoPRM, monoCOM, monoMAT, monoPREC, X, Y)
    implicit none
    !> パラメータ構造体
    type(monolis_prm) :: monoPRM
    !> 通信テーブル構造体
    type(monolis_com) :: monoCOM
    !> 行列構造体
    type(monolis_mat) :: monoMAT
    !> 前処理構造体
    type(monolis_mat) :: monoPREC
    integer(kint) :: i, precond
    real(kdouble) :: X(:), Y(:)
    real(kdouble) :: t1, t2

#ifdef DEBUG
    call monolis_std_debug_log_header("monolis_precond_apply_R")
#endif

    t1 = monolis_get_time()

    precond = monoPRM%Iarray(monolis_prm_I_precond)

    if(precond == monolis_prec_DIAG)then
      call monolis_precond_diag_apply(monoPRM, monoCOM, monoMAT, monoPREC, X, Y)
    !elseif(precond == monolis_prec_ILU)then
    !  call monolis_precond_ilu_apply(monoPRM, monoCOM, monoMAT, X, Y)
    !elseif(precond == monolis_prec_JACOBI)then
    !  call monolis_precond_Jacobi_apply(monoPRM, monoCOM, monoMAT, X, Y)
    elseif(precond == monolis_prec_SOR)then
      call monolis_precond_SOR_apply(monoPRM, monoCOM, monoMAT, monoPREC, X, Y)
    !elseif(precond == monolis_prec_MUMPS)then
    !  call monolis_precond_MUMPS_apply(monoPRM, monoCOM, monoMAT, X, Y)
    !elseif(precond == monolis_prec_MUMPS_LOCAL)then
    !  call monolis_precond_MUMPS_apply(monoPRM, monoCOM, monoMAT, X, Y)
    !elseif(precond == monolis_prec_MF)then
    !  call monolis_precond_MF_apply(monoPRM, monoCOM, monoMAT, X, Y)
    elseif(precond == monolis_prec_NONE)then
      do i = 1, monoMAT%CSR%N*monoMAT%CSR%NDOF
        Y(i) = X(i)
      enddo
    else
      stop "monolis_precond_apply_R"
    endif

    t2 = monolis_get_time()
    !monoPRM%tprec = monoPRM%tprec + t2 - t1
  end subroutine monolis_precond_apply_R

  !> 前処理初期化関数
  subroutine monolis_precond_clear_R(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> パラメータ構造体
    type(monolis_prm) :: monoPRM
    !> 通信テーブル構造体
    type(monolis_com) :: monoCOM
    !> 行列構造体
    type(monolis_mat) :: monoMAT
    !> 前処理構造体
    type(monolis_mat) :: monoPREC
    integer(kint) :: precond
    real(kdouble) :: t1, t2

#ifdef DEBUG
    call monolis_std_debug_log_header("monolis_precond_clear_R")
#endif

    t1 = monolis_get_time()

    precond = monoPRM%Iarray(monolis_prm_I_precond)

    if(precond == monolis_prec_DIAG)then
      call monolis_precond_diag_clear(monoPRM, monoCOM, monoMAT, monoPREC)
    !elseif(precond == monolis_prec_ILU)then
    !  call monolis_precond_ilu_clear(monoPRM, monoCOM, monoMAT)
    !elseif(precond == monolis_prec_JACOBI)then
    !  call monolis_precond_Jacobi_clear(monoPRM, monoCOM, monoMAT)
    !elseif(precond == monolis_prec_SOR)then
      call monolis_precond_SOR_clear(monoPRM, monoCOM, monoMAT, monoPREC)
    !elseif(precond == monolis_prec_MUMPS)then
    !  call monolis_precond_MUMPS_clear(monoPRM, monoCOM, monoMAT)
    !elseif(precond == monolis_prec_MUMPS_LOCAL)then
    !  call monolis_precond_MUMPS_clear(monoPRM, monoCOM, monoMAT)
    !elseif(precond == monolis_prec_MF)then
    !  call monolis_precond_MF_clear(monoPRM, monoCOM, monoMAT)
    endif

    t2 = monolis_get_time()
    !monoPRM%tprep = monoPRM%tprep + t2 - t1
  end subroutine monolis_precond_clear_R

end module mod_monolis_precond