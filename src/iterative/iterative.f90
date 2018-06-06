module mod_monolis_iterative
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_solver_CG
  use mod_monolis_solver_BiCGSTAB
  use mod_monolis_solver_BiCGSTAB_noprec
  use mod_monolis_solver_GropCG
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
    end select

  end subroutine monolis_iterative

end module mod_monolis_iterative