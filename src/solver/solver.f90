!> 線形ソルバモジュール
module mod_monolis_solve
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_solver_CG
  use mod_monolis_solver_BiCGSTAB
  use mod_monolis_solver_BiCGSTAB_noprec
  use mod_monolis_solver_GropCG
  use mod_monolis_solver_PipeCG
  use mod_monolis_solver_PipeCR
  use mod_monolis_solver_PipeBiCGSTAB
  use mod_monolis_solver_PipeBiCGSTAB_noprec
  use mod_monolis_solver_DeflatedCG
  use mod_monolis_solver_COCG
  use mod_monolis_precond

  implicit none

contains

  !> @ingroup solver
  !> 線形ソルバ関数（実数型）
  subroutine monolis_solve_R(monolis, monoCOM, B, X)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 右辺ベクトル
    real(kdouble), intent(in) :: B(:)
    !> [in,out] 解ベクトル
    real(kdouble), intent(inout) :: X(:)

    call monolis_set_RHS_R(monolis%MAT, B)

    call monolis_set_initial_solution_R(monolis%MAT, X)

    if(monoCOM%comm_size > 1) monolis%MAT%N = monoCOM%n_internal_vertex

    call monolis_solve_main_R(monolis%PRM, monoCOM, monolis%MAT, monolis%PREC)

    call monolis_get_solution_R(monolis%MAT, X)
  end subroutine monolis_solve_R

  !> @ingroup solver
  !> 線形ソルバ関数（複素数型）
  subroutine monolis_solve_C(monolis, monoCOM, B, X)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 右辺ベクトル
    complex(kdouble), intent(in) :: B(:)
    !> [in,out] 解ベクトル
    complex(kdouble), intent(inout) :: X(:)

    call monolis_set_RHS_C(monolis%MAT, B)

    call monolis_set_initial_solution_C(monolis%MAT, X)

    if(monoCOM%comm_size > 1) monolis%MAT%N = monoCOM%n_internal_vertex

    call monolis_solve_main_C(monolis%PRM, monoCOM, monolis%MAT, monolis%PREC)

    call monolis_get_solution_C(monolis%MAT, X)
  end subroutine monolis_solve_C

  !> @ingroup dev_solver
  !> 線形ソルバ関数（メイン関数）
  subroutine monolis_solve_main_R(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in,out] 行列構造体
    type(monolis_mat), intent(inout) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC

    call monolis_std_debug_log_header("monolis_solve_main_R")

    call monolis_timer_initialize(monoPRM)

    !call monolis_check_input_param(monoCOM, monoMAT)

    !call monolis_matrix_reordering_fw_R()

    call monolis_precond_setup(monoPRM, monoCOM, monoMAT, monoPREC)

    call monolis_solver_select_R(monoPRM, monoCOM, monoMAT, monoPREC)

    call monolis_precond_clear(monoPRM, monoCOM, monoMAT, monoPREC)

    !call monolis_matrix_reordering_bk_R()

    call monolis_timer_finalize(monoPRM, monoCOM)
  end subroutine monolis_solve_main_R

  !> @ingroup dev_solver
  !> 線形ソルバ関数（メイン関数）
  subroutine monolis_solve_main_C(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in,out] 行列構造体
    type(monolis_mat), intent(inout) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC

    call monolis_std_debug_log_header("monolis_solve_main_C")

    call monolis_timer_initialize(monoPRM)

    !call monolis_check_input_param(monoCOM, monoMAT)

    call monolis_precond_setup(monoPRM, monoCOM, monoMAT, monoPREC)

    call monolis_solver_select_C(monoPRM, monoCOM, monoMAT, monoPREC)

    call monolis_precond_clear(monoPRM, monoCOM, monoMAT, monoPREC)

    call monolis_timer_finalize(monoPRM, monoCOM)
  end subroutine monolis_solve_main_C

  !> @ingroup dev_solver
  !> 線形ソルバの選択関数
  subroutine monolis_solver_select_R(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in,out] 行列構造体
    type(monolis_mat), intent(inout) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC

    call monolis_std_debug_log_header("monolis_solver_select_R")

    if(monoPRM%Iarray(monolis_prm_I_show_summary) == monolis_I_true .and. monoCOM%my_rank == 0)then
      write(*,"(a)") &
      & "** MONOLIS solver: "//trim(monolis_str_iter(monoPRM%Iarray(monolis_prm_I_method)))//&
      & ", prec: "//trim(monolis_str_prec(monoPRM%Iarray(monolis_prm_I_precond)))
    endif

    select case(monoPRM%Iarray(monolis_prm_I_method))
      case (monolis_iter_CG)
        call monolis_solver_CG(monoPRM, monoCOM, monoMAT, monoPREC)

      case (monolis_iter_BiCGSTAB)
        call monolis_solver_BiCGSTAB(monoPRM, monoCOM, monoMAT, monoPREC)

      case (monolis_iter_BiCGSTAB_noprec)
        call monolis_solver_BiCGSTAB_noprec(monoPRM, monoCOM, monoMAT, monoPREC)

      case (monolis_iter_GropCG)
        call monolis_solver_GropCG(monoPRM, monoCOM, monoMAT, monoPREC)

      case (monolis_iter_PipeCG)
        call monolis_solver_PipeCG(monoPRM, monoCOM, monoMAT, monoPREC)

      case (monolis_iter_PipeCR)
        call monolis_solver_PipeCR(monoPRM, monoCOM, monoMAT, monoPREC)

      !case (monolis_iter_CABiCGSTAB_noprec)
      !  call monolis_solver_CABiCGSTAB_noprec(monoPRM, monoCOM, monoMAT)

      case (monolis_iter_PipeBiCGSTAB)
        call monolis_solver_PipeBiCGSTAB(monoPRM, monoCOM, monoMAT, monoPREC)

      case (monolis_iter_PipeBiCGSTAB_noprec)
        call monolis_solver_PipeBiCGSTAB_noprec(monoPRM, monoCOM, monoMAT, monoPREC)

      case (monolis_iter_DeflatedCG1)
        call monolis_solver_DeflatedCG1(monoPRM, monoCOM, monoMAT, monoPREC)

      case (monolis_iter_DeflatedCG2)
        call monolis_solver_DeflatedCG2(monoPRM, monoCOM, monoMAT, monoPREC)

      case (monolis_iter_ADeflatedCG2)
        call monolis_solver_ADeflatedCG2(monoPRM, monoCOM, monoMAT, monoPREC)

      !case (monolis_iter_GMRES)
      !  call monolis_solver_GMRES(monoPRM, monoCOM, monoMAT)

      case default
        call monolis_std_error_string("monolis_solver_select_R")
        call monolis_std_error_string("please select a linear solver for real system")
        call monolis_std_error_stop()
    end select
  end subroutine monolis_solver_select_R

  !> @ingroup dev_solver
  !> 線形ソルバの選択関数
  subroutine monolis_solver_select_C(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in,out] 行列構造体
    type(monolis_mat), intent(inout) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC

    call monolis_std_debug_log_header("monolis_solver_select_R")

    if(monoPRM%Iarray(monolis_prm_I_show_summary) == monolis_I_true .and. monoCOM%my_rank == 0)then
      write(*,"(a)") &
      & "** MONOLIS solver: "//trim(monolis_str_iter(monoPRM%Iarray(monolis_prm_I_method)))//&
      & ", prec: "//trim(monolis_str_prec(monoPRM%Iarray(monolis_prm_I_precond)))
    endif

    select case(monoPRM%Iarray(monolis_prm_I_method))
      case (monolis_iter_COCG)
        call monolis_solver_COCG(monoPRM, monoCOM, monoMAT, monoPREC)

      case default
        call monolis_std_error_string("monolis_solver_select_C")
        call monolis_std_error_string("please select a linear solver for complex system")
        call monolis_std_error_stop()
    end select
  end subroutine monolis_solver_select_C
end module mod_monolis_solve
