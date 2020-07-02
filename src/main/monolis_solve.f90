module mod_monolis_solve
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_scaling
  use mod_monolis_precond
  use mod_monolis_reorder
  use mod_monolis_util
  use mod_monolis_solver_CG
  use mod_monolis_solver_BiCGSTAB
  use mod_monolis_solver_BiCGSTAB_noprec
  use mod_monolis_solver_GropCG
  use mod_monolis_solver_PipeCG
  use mod_monolis_solver_PipeCR
  use mod_monolis_solver_CABiCGSTAB_noprec
  use mod_monolis_solver_PipeBiCGSTAB
  use mod_monolis_solver_PipeBiCGSTAB_noprec
  use mod_monolis_solver_SOR
  use mod_monolis_solver_IR
  implicit none

contains

  subroutine monolis_solve(monolis, B, X)
    implicit none
    type(monolis_structure) :: monolis
    real(kdouble) :: B(:), X(:)
    integer(kint) :: i

    do i = 1, monolis%MAT%NP*monolis%MAT%NDOF
      monolis%MAT%B(i) = B(i)
      monolis%MAT%X(i) = X(i)
    enddo

    call monolis_solve_(monolis%PRM, monolis%COM, monolis%MAT)

    do i = 1, monolis%MAT%NP*monolis%MAT%NDOF
      X(i) = monolis%MAT%X(i)
    enddo
  end subroutine monolis_solve

  subroutine monolis_solve_(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_com) :: monoCOM_reorder
    type(monolis_mat) :: monoMAT
    type(monolis_mat) :: monoMAT_reorder

#ifdef DEBUG
    monoPRM%is_debug = .true.
#endif

    call monolis_timer_initialize(monoPRM)
    call monolis_check_diagonal(monoPRM, monoMAT)
    call monolis_reorder_matrix_fw(monoPRM, monoCOM, monoCOM_reorder, monoMAT, monoMAT_reorder)
    call monolis_scaling_fw(monoPRM, monoCOM_reorder, monoMAT_reorder)
    call monolis_precond_setup(monoPRM, monoCOM_reorder, monoMAT_reorder)
    call monolis_solver(monoPRM, monoCOM_reorder, monoMAT_reorder)
    call monolis_precond_clear(monoPRM, monoCOM_reorder, monoMAT_reorder)
    call monolis_scaling_bk(monoPRM, monoCOM_reorder, monoMAT_reorder)
    call monolis_reorder_matrix_bk(monoPRM, monoCOM_reorder, monoMAT_reorder, monoMAT)
    call monolis_timer_finalize(monoPRM, monoCOM)
  end subroutine monolis_solve_

  subroutine monolis_solve_test(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_com) :: monoCOM_reorder
    type(monolis_mat) :: monoMAT
    type(monolis_mat) :: monoMAT_reorder
    integer(kint) :: i, j

    call monolis_check_diagonal(monoPRM, monoMAT)
    do i = 1, 9
      monoPRM%method = i
      do j = 1, 4
        monoPRM%precond = j
        call monolis_timer_initialize(monoPRM)
        call monolis_reorder_matrix_fw(monoPRM, monoCOM, monoCOM_reorder, monoMAT, monoMAT_reorder)
        call monolis_scaling_fw(monoPRM, monoCOM_reorder, monoMAT_reorder)
        call monolis_precond_setup(monoPRM, monoCOM_reorder, monoMAT_reorder)
        call monolis_solver(monoPRM, monoCOM_reorder, monoMAT_reorder)
        call monolis_precond_clear(monoPRM, monoCOM_reorder, monoMAT_reorder)
        call monolis_scaling_bk(monoPRM, monoCOM_reorder, monoMAT_reorder)
        call monolis_reorder_matrix_bk(monoPRM, monoCOM_reorder, monoMAT_reorder, monoMAT)
        call monolis_timer_finalize(monoPRM, monoCOM)
      enddo
    enddo
  end subroutine monolis_solve_test

  subroutine monolis_solver(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    if(monoPRM%is_debug) call monolis_debug_header("monolis_solver")

    if(monoPRM%show_summary .and. monoCOM%myrank == 0) write(*,"(a)")" ** monolis solver: "// &
    & trim(monolis_str_iter(monoPRM%method))//", prec: "//trim(monolis_str_prec(monoPRM%precond))

    select case(monoPRM%method)
      case (monolis_iter_CG)
        call monolis_solver_CG(monoPRM, monoCOM, monoMAT)

      case (monolis_iter_BiCGSTAB)
        call monolis_solver_BiCGSTAB(monoPRM, monoCOM, monoMAT)

      case (monolis_iter_BiCGSTAB_noprec)
        call monolis_solver_BiCGSTAB_noprec(monoPRM, monoCOM, monoMAT)

      case (monolis_iter_GropCG)
        call monolis_solver_GropCG(monoPRM, monoCOM, monoMAT)

      case (monolis_iter_PipeCG)
        call monolis_solver_PipeCG(monoPRM, monoCOM, monoMAT)

      case (monolis_iter_PipeCR)
        call monolis_solver_PipeCR(monoPRM, monoCOM, monoMAT)

      case (monolis_iter_CABiCGSTAB_noprec)
        call monolis_solver_CABiCGSTAB_noprec(monoPRM, monoCOM, monoMAT)

      case (monolis_iter_PipeBiCGSTAB)
        call monolis_solver_PipeBiCGSTAB(monoPRM, monoCOM, monoMAT)

      case (monolis_iter_PipeBiCGSTAB_noprec)
        call monolis_solver_PipeBiCGSTAB_noprec(monoPRM, monoCOM, monoMAT)

      case (monolis_iter_SOR)
        call monolis_solver_SOR(monoPRM, monoCOM, monoMAT)

      case (monolis_iter_IR)
        call monolis_solver_IR(monoPRM, monoCOM, monoMAT)

    end select
  end subroutine monolis_solver

end module mod_monolis_solve
