!> 前処理関連モジュール
module mod_monolis_precond
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_precond_diag
  use mod_monolis_precond_SOR

  implicit none

contains

  !> @ingroup prec
  !> 前処理生成関数
  subroutine monolis_precond_setup(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> パラメータ構造体
    type(monolis_prm) :: monoPRM
    !> 通信テーブル構造体
    type(monolis_com) :: monoCOM
    !> 行列構造体
    type(monolis_mat) :: monoMAT
    !> 前処理構造体
    type(monolis_mat) :: monoPREC

    select case(monoPRM%Iarray(monolis_prm_I_method))
      case (monolis_iter_COCG)
        call monolis_precond_setup_C(monoPRM, monoCOM, monoMAT, monoPREC)
      case default
        call monolis_precond_setup_R(monoPRM, monoCOM, monoMAT, monoPREC)
    end select
  end subroutine monolis_precond_setup

  !> @ingroup prec
  !> 前処理初期化関数
  subroutine monolis_precond_clear(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> パラメータ構造体
    type(monolis_prm) :: monoPRM
    !> 通信テーブル構造体
    type(monolis_com) :: monoCOM
    !> 行列構造体
    type(monolis_mat) :: monoMAT
    !> 前処理構造体
    type(monolis_mat) :: monoPREC

    select case(monoPRM%Iarray(monolis_prm_I_method))
      case (monolis_iter_COCG)
        call monolis_precond_clear_C(monoPRM, monoCOM, monoMAT, monoPREC)
      case default
        call monolis_precond_clear_R(monoPRM, monoCOM, monoMAT, monoPREC)
    end select
  end subroutine monolis_precond_clear

  !> @ingroup prec
  !> 前処理生成関数（実数型）
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

    call monolis_std_debug_log_header("monolis_precond_setup_R")

    t1 = monolis_get_time()

    precond = monoPRM%Iarray(monolis_prm_I_precond)

    if(precond == monolis_prec_DIAG)then
      call monolis_precond_diag_setup_R(monoPRM, monoCOM, monoMAT, monoPREC)
    !elseif(precond == monolis_prec_ILU)then
    !  call monolis_precond_ilu_setup(monoPRM, monoCOM, monoMAT)
    !elseif(precond == monolis_prec_JACOBI)then
    !  call monolis_precond_Jacobi_setup(monoPRM, monoCOM, monoMAT)
    elseif(precond == monolis_prec_SOR)then
      call monolis_precond_SOR_setup_R(monoPRM, monoCOM, monoMAT, monoPREC)
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

  !> @ingroup prec
  !> 前処理生成関数（複素数型）
  subroutine monolis_precond_setup_C(monoPRM, monoCOM, monoMAT, monoPREC)
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

    call monolis_std_debug_log_header("monolis_precond_setup_C")

    t1 = monolis_get_time()

    precond = monoPRM%Iarray(monolis_prm_I_precond)

    if(precond == monolis_prec_DIAG)then
      call monolis_precond_diag_setup_C(monoPRM, monoCOM, monoMAT, monoPREC)
    !elseif(precond == monolis_prec_ILU)then
    !  call monolis_precond_ilu_setup(monoPRM, monoCOM, monoMAT)
    !elseif(precond == monolis_prec_JACOBI)then
    !  call monolis_precond_Jacobi_setup(monoPRM, monoCOM, monoMAT)
    elseif(precond == monolis_prec_SOR)then
      call monolis_precond_SOR_setup_C(monoPRM, monoCOM, monoMAT, monoPREC)
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
  end subroutine monolis_precond_setup_C

  !> @ingroup prec
  !> 前処理適用関数（実数型）
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

    call monolis_std_debug_log_header("monolis_precond_apply_R")

    t1 = monolis_get_time()

    precond = monoPRM%Iarray(monolis_prm_I_precond)

    if(precond == monolis_prec_DIAG)then
      call monolis_precond_diag_apply_R(monoPRM, monoCOM, monoMAT, monoPREC, X, Y)
    !elseif(precond == monolis_prec_ILU)then
    !  call monolis_precond_ilu_apply(monoPRM, monoCOM, monoMAT, X, Y)
    !elseif(precond == monolis_prec_JACOBI)then
    !  call monolis_precond_Jacobi_apply(monoPRM, monoCOM, monoMAT, X, Y)
    elseif(precond == monolis_prec_SOR)then
      call monolis_precond_SOR_apply_R(monoPRM, monoCOM, monoMAT, monoPREC, X, Y)
    !elseif(precond == monolis_prec_MUMPS)then
    !  call monolis_precond_MUMPS_apply(monoPRM, monoCOM, monoMAT, X, Y)
    !elseif(precond == monolis_prec_MUMPS_LOCAL)then
    !  call monolis_precond_MUMPS_apply(monoPRM, monoCOM, monoMAT, X, Y)
    !elseif(precond == monolis_prec_MF)then
    !  call monolis_precond_MF_apply(monoPRM, monoCOM, monoMAT, X, Y)
    elseif(precond == monolis_prec_NONE)then
      do i = 1, monoMAT%N*monoMAT%NDOF
        Y(i) = X(i)
      enddo
    else
      write(*,*) "precond", precond
      stop "monolis_precond_apply_R"
    endif

    t2 = monolis_get_time()
    !monoPRM%tprec = monoPRM%tprec + t2 - t1
  end subroutine monolis_precond_apply_R

  !> @ingroup prec
  !> 前処理適用関数（複素数型）
  subroutine monolis_precond_apply_C(monoPRM, monoCOM, monoMAT, monoPREC, X, Y)
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
    complex(kdouble) :: X(:), Y(:)
    real(kdouble) :: t1, t2

    call monolis_std_debug_log_header("monolis_precond_apply_C")

    t1 = monolis_get_time()

    precond = monoPRM%Iarray(monolis_prm_I_precond)

    if(precond == monolis_prec_DIAG)then
      call monolis_precond_diag_apply_C(monoPRM, monoCOM, monoMAT, monoPREC, X, Y)
    !elseif(precond == monolis_prec_ILU)then
    !  call monolis_precond_ilu_apply(monoPRM, monoCOM, monoMAT, X, Y)
    !elseif(precond == monolis_prec_JACOBI)then
    !  call monolis_precond_Jacobi_apply(monoPRM, monoCOM, monoMAT, X, Y)
    elseif(precond == monolis_prec_SOR)then
      call monolis_precond_SOR_apply_C(monoPRM, monoCOM, monoMAT, monoPREC, X, Y)
    !elseif(precond == monolis_prec_MUMPS)then
    !  call monolis_precond_MUMPS_apply(monoPRM, monoCOM, monoMAT, X, Y)
    !elseif(precond == monolis_prec_MUMPS_LOCAL)then
    !  call monolis_precond_MUMPS_apply(monoPRM, monoCOM, monoMAT, X, Y)
    !elseif(precond == monolis_prec_MF)then
    !  call monolis_precond_MF_apply(monoPRM, monoCOM, monoMAT, X, Y)
    elseif(precond == monolis_prec_NONE)then
      do i = 1, monoMAT%N*monoMAT%NDOF
        Y(i) = X(i)
      enddo
    else
      write(*,*) "precond", precond
      stop "monolis_precond_apply_C"
    endif

    t2 = monolis_get_time()
    !monoPRM%tprec = monoPRM%tprec + t2 - t1
  end subroutine monolis_precond_apply_C

  !> @ingroup prec
  !> 前処理初期化関数（実数型）
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

    call monolis_std_debug_log_header("monolis_precond_clear_R")

    t1 = monolis_get_time()

    precond = monoPRM%Iarray(monolis_prm_I_precond)

    if(precond == monolis_prec_DIAG)then
      call monolis_precond_diag_clear_R(monoPRM, monoCOM, monoMAT, monoPREC)
    !elseif(precond == monolis_prec_ILU)then
    !  call monolis_precond_ilu_clear(monoPRM, monoCOM, monoMAT)
    !elseif(precond == monolis_prec_JACOBI)then
    !  call monolis_precond_Jacobi_clear(monoPRM, monoCOM, monoMAT)
    elseif(precond == monolis_prec_SOR)then
      call monolis_precond_SOR_clear_R(monoPRM, monoCOM, monoMAT, monoPREC)
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

  !> @ingroup prec
  !> 前処理初期化関数（複素数型）
  subroutine monolis_precond_clear_C(monoPRM, monoCOM, monoMAT, monoPREC)
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

    call monolis_std_debug_log_header("monolis_precond_clear_C")

    t1 = monolis_get_time()

    precond = monoPRM%Iarray(monolis_prm_I_precond)

    if(precond == monolis_prec_DIAG)then
      call monolis_precond_diag_clear_C(monoPRM, monoCOM, monoMAT, monoPREC)
    !elseif(precond == monolis_prec_ILU)then
    !  call monolis_precond_ilu_clear(monoPRM, monoCOM, monoMAT)
    !elseif(precond == monolis_prec_JACOBI)then
    !  call monolis_precond_Jacobi_clear(monoPRM, monoCOM, monoMAT)
    elseif(precond == monolis_prec_SOR)then
      call monolis_precond_SOR_clear_C(monoPRM, monoCOM, monoMAT, monoPREC)
    !elseif(precond == monolis_prec_MUMPS)then
    !  call monolis_precond_MUMPS_clear(monoPRM, monoCOM, monoMAT)
    !elseif(precond == monolis_prec_MUMPS_LOCAL)then
    !  call monolis_precond_MUMPS_clear(monoPRM, monoCOM, monoMAT)
    !elseif(precond == monolis_prec_MF)then
    !  call monolis_precond_MF_clear(monoPRM, monoCOM, monoMAT)
    endif

    t2 = monolis_get_time()
    !monoPRM%tprep = monoPRM%tprep + t2 - t1
  end subroutine monolis_precond_clear_C
end module mod_monolis_precond