module mod_monolis_iterative
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_solver_CG
  use mod_monolis_solver_BiCGSTAB
  use mod_monolis_solver_BiCGSTAB_noprec
  use mod_monolis_solver_GropCG
  use mod_monolis_solver_PipeCG
  use mod_monolis_solver_PipeCR
  use mod_monolis_solver_CABiCGSTAB_noprec
  use mod_monolis_solver_SOR
  implicit none

contains

  subroutine monolis_iterative(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    select case(monoPRM%method)
      case (monolis_iter_CG)
        if(monoCOM%myrank == 0) write(*,"(a)")" ** monolis_solver_CG"
        call monolis_solver_CG(monoPRM, monoCOM, monoMAT)

      case (monolis_iter_BiCGSTAB)
        if(monoCOM%myrank == 0) write(*,"(a)")" ** monolis_solver_BiCGSTAB"
        call monolis_solver_BiCGSTAB(monoPRM, monoCOM, monoMAT)

      case (monolis_iter_BiCGSTAB_noprec)
        if(monoCOM%myrank == 0) write(*,"(a)")" ** monolis_solver_BiCGSTAB_noprec"
        call monolis_solver_BiCGSTAB_noprec(monoPRM, monoCOM, monoMAT)

      case (monolis_iter_GropCG)
        if(monoCOM%myrank == 0) write(*,"(a)")" ** monolis_solver_GropCG"
        call monolis_solver_GropCG(monoPRM, monoCOM, monoMAT)

      case (monolis_iter_PipeCG)
        if(monoCOM%myrank == 0) write(*,"(a)")" ** monolis_solver_PipeCG"
        call monolis_solver_PipeCG(monoPRM, monoCOM, monoMAT)

      case (monolis_iter_PipeCR)
        if(monoCOM%myrank == 0) write(*,"(a)")" ** monolis_solver_PipeCR"
        call monolis_solver_PipeCR(monoPRM, monoCOM, monoMAT)

      case (monolis_iter_CABiCGSTAB_noprec)
        if(monoCOM%myrank == 0) write(*,"(a)")" ** monolis_solver_CABiCGSTAB_noprec"
        call monolis_solver_CABiCGSTAB_noprec(monoPRM, monoCOM, monoMAT)

      case (monolis_iter_SOR)
        if(monoCOM%myrank == 0) write(*,"(a)")" ** monolis_solver_SOR"
        call monolis_solver_SOR(monoPRM, monoCOM, monoMAT)

      !case (monolis_iter_IR)
      !  if(monoCOM%myrank == 0) write(*,"(a)")" ** monolis_solver_IR"
      !  call monolis_solver_IR(monoPRM, monoCOM, monoMAT)
    end select

  end subroutine monolis_iterative

end module mod_monolis_iterative