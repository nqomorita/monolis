!> 線形ソルバモジュール
module mod_monolis_solve
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_solver_CG
  use mod_monolis_solver_BiCGSTAB
  use mod_monolis_solver_BiCGSTAB_noprec
  use mod_monolis_solver_GropCG
  use mod_monolis_solver_COCG
  use mod_monolis_precond

  implicit none

contains

  !> 線形ソルバ関数（実数型）
  subroutine monolis_solve_R(monolis, B, X)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> 右辺ベクトル
    real(kdouble) :: B(:)
    !> 解ベクトル
    real(kdouble) :: X(:)

    call monolis_set_RHS_R(monolis%MAT, B)

    call monolis_set_initial_solution_R(monolis%MAT, X)

    call monolis_solve_main(monolis%PRM, monolis%COM, monolis%MAT, monolis%PREC)

    call monolis_get_solution_R(monolis%MAT, X)
  end subroutine monolis_solve_R

  !> 線形ソルバ関数（複素数型）
  subroutine monolis_solve_C(monolis, B, X)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> 右辺ベクトル
    complex(kdouble) :: B(:)
    !> 解ベクトル
    complex(kdouble) :: X(:)

    call monolis_set_RHS_C(monolis%MAT, B)

    call monolis_set_initial_solution_C(monolis%MAT, X)

    call monolis_solve_main(monolis%PRM, monolis%COM, monolis%MAT, monolis%PREC)

    call monolis_get_solution_C(monolis%MAT, X)
  end subroutine monolis_solve_C

  !> 線形ソルバ関数（メイン関数）
  subroutine monolis_solve_main(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> パラメータ構造体
    type(monolis_prm) :: monoPRM
    !> 通信テーブル構造体
    type(monolis_com) :: monoCOM
    !> 行列構造体
    type(monolis_mat) :: monoMAT
    !> 前処理構造体
    type(monolis_mat) :: monoPREC

    call monolis_std_debug_log_header("monolis_solve_main")

    call monolis_timer_initialize(monoPRM, monoCOM)

    call monolis_check_input_param(monoPRM, monoCOM, monoMAT)

    call monolis_precond_setup(monoPRM, monoCOM, monoMAT, monoPREC)

    call monolis_solver(monoPRM, monoCOM, monoMAT, monoPREC)

    call monolis_precond_clear(monoPRM, monoCOM, monoMAT, monoPREC)

    call monolis_timer_finalize(monoPRM, monoCOM)
  end subroutine monolis_solve_main

  subroutine monolis_solver(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> パラメータ構造体
    type(monolis_prm) :: monoPRM
    !> 通信テーブル構造体
    type(monolis_com) :: monoCOM
    !> 行列構造体
    type(monolis_mat) :: monoMAT
    !> 前処理構造体
    type(monolis_mat) :: monoPREC

    call monolis_std_debug_log_header("monolis_solver")

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

      !case (monolis_iter_PipeCG)
      !  call monolis_solver_PipeCG(monoPRM, monoCOM, monoMAT)

      !case (monolis_iter_PipeCR)
      !  call monolis_solver_PipeCR(monoPRM, monoCOM, monoMAT)

      !case (monolis_iter_CABiCGSTAB_noprec)
      !  call monolis_solver_CABiCGSTAB_noprec(monoPRM, monoCOM, monoMAT)

      !case (monolis_iter_PipeBiCGSTAB)
      !  call monolis_solver_PipeBiCGSTAB(monoPRM, monoCOM, monoMAT)

      !case (monolis_iter_PipeBiCGSTAB_noprec)
      !  call monolis_solver_PipeBiCGSTAB_noprec(monoPRM, monoCOM, monoMAT)

      !case (monolis_iter_SOR)
      !  call monolis_solver_SOR(monoPRM, monoCOM, monoMAT)

      !case (monolis_iter_IR)
      !  call monolis_solver_IR(monoPRM, monoCOM, monoMAT)

      !case (monolis_iter_GMRES)
      !  call monolis_solver_GMRES(monoPRM, monoCOM, monoMAT)

      case (monolis_iter_COCG)
        call monolis_solver_COCG(monoPRM, monoCOM, monoMAT, monoPREC)
    end select
  end subroutine monolis_solver

end module mod_monolis_solve
