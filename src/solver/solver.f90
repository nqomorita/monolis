!> 線形ソルバモジュール
module mod_monolis_solve
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_solver_CG
  use mod_monolis_solver_BiCGSTAB
  use mod_monolis_solver_COCG

  implicit none

contains

  subroutine monolis_solve(monolis, B, X)
    implicit none
    type(monolis_structure) :: monolis
    real(kdouble) :: B(:), X(:)
    integer(kint) :: i

    !call monolis_set_RHS(monolis%MAT, B)

    !call monolis_set_initial_solution(monolis%MAT, X)

    !call monolis_set_initial_comm(monolis%COM, monolis%MAT)

    call monolis_solve_main(monolis%PRM, monolis%COM, monolis%MAT, monolis%PREC)

    !call monolis_get_solution(monolis%MAT, X)
  end subroutine monolis_solve

  subroutine monolis_set_RHS(monoMAT, B)
    implicit none
    type(monolis_mat) :: monoMAT
    real(kdouble) :: B(:)
    integer(kint) :: i

    do i = 1, monoMAT%NP*monoMAT%NDOF
      monoMAT%R%B(i) = B(i)
    enddo
  end subroutine monolis_set_RHS

  subroutine monolis_set_initial_solution(monoMAT, X)
    implicit none
    type(monolis_mat) :: monoMAT
    real(kdouble) :: X(:)
    integer(kint) :: i

    do i = 1, monoMAT%NP*monoMAT%NDOF
      monoMAT%R%X(i) = X(i)
    enddo
  end subroutine monolis_set_initial_solution

  subroutine monolis_get_solution(monoMAT, X)
    implicit none
    type(monolis_mat) :: monoMAT
    real(kdouble) :: X(:)
    integer(kint) :: i

    do i = 1, monoMAT%NP*monoMAT%NDOF
      X(i) = monoMAT%R%X(i)
    enddo
  end subroutine monolis_get_solution

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

    !call monolis_timer_initialize(monoPRM, monoCOM)
    !call monolis_check_diagonal(monoPRM, monoMAT)
    !call monolis_reorder_matrix_fw(monoPRM, monoCOM, monoCOM_reorder, monoMAT, monoMAT_reorder)
    !call monolis_scaling_fw(monoPRM, monoCOM_reorder, monoMAT_reorder)
    !call monolis_precond_setup(monoPRM, monoCOM_reorder, monoMAT_reorder)
    call monolis_solver(monoPRM, monoCOM, monoMAT, monoPREC)
    !call monolis_precond_clear(monoPRM, monoCOM_reorder, monoMAT_reorder)
    !call monolis_scaling_bk(monoPRM, monoCOM_reorder, monoMAT_reorder)
    !call monolis_reorder_matrix_bk(monoPRM, monoCOM_reorder, monoMAT_reorder, monoMAT)
    !call monolis_timer_finalize(monoPRM, monoCOM)
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

!    if(monoPRM%is_debug) call monolis_std_debug_log_header("monolis_solver v0.0.0")

!    if(monoPRM%show_summary .and. monoCOM%my_rank == 0) write(*,"(a)")" ** monolis solver: "// &
!    & trim(monolis_str_iter(monoPRM%method))//", prec: "//trim(monolis_str_prec(monoPRM%precond))

    select case(monoPRM%Iarray(monolis_prm_I_method))
      case (monolis_iter_CG)
        call monolis_solver_CG(monoPRM, monoCOM, monoMAT, monoPREC)

      case (monolis_iter_BiCGSTAB)
        call monolis_solver_BiCGSTAB(monoPRM, monoCOM, monoMAT, monoPREC)

      !case (monolis_iter_BiCGSTAB_noprec)
      !  call monolis_solver_BiCGSTAB_noprec(monoPRM, monoCOM, monoMAT)

      !case (monolis_iter_GropCG)
      !  call monolis_solver_GropCG(monoPRM, monoCOM, monoMAT)

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

  subroutine monolis_timer_finalize(monoPRM, monoCOM)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    real(kdouble) :: t1, time(6), t_max, t_min, t_avg, t_sd
    logical :: is_output

!    call monolis_std_debug_log_header("monolis_timer_finalize")
!
!    call monolis_barrier_(monoCOM%comm)
!    t1 = monolis_get_time()
!    monoPRM%tsol = t1 - monoPRM%tsol
!
!    if(monoPRM%show_summary .and. monoCOM%my_rank == 0)then
!      write(*,"(a,i10)")" ** monolis converge iter:", monoPRM%curiter
!      write(*,"(a,1p4e10.3)")" ** monolis rel. residual:", monoPRM%curresid
!    endif
!
!    is_output = monoPRM%show_summary .or. monoPRM%show_time
!    if(is_output .and. monoCOM%my_rank == 0)then
!      write(*,"(a,1p4e10.3)")" ** monolis solution time:", monoPRM%tsol
!    endif
!
!    if(monoPRM%show_time)then
!      time(1) = monoPRM%tprep
!      time(2) = monoPRM%tspmv
!      time(3) = monoPRM%tdotp
!      time(4) = monoPRM%tprec
!      time(5) = monoPRM%tcomm_dotp
!      time(6) = monoPRM%tcomm_spmv
!      call monolis_allreduce_R(6, time, monolis_mpi_sum, monoCOM%comm)
!      time = time/dble(monolis_mpi_global_comm_size())
!
!      if(monoCOM%my_rank == 0)then
!        write(*,"(a,1p4e10.3)")"  - solution/prepost time:", time(1)
!        write(*,"(a,1p4e10.3)")"  - solution/SpMV    time:", time(2)
!        write(*,"(a,1p4e10.3)")"  - solution/inner p time:", time(3)
!        write(*,"(a,1p4e10.3)")"  - solution/precond time:", time(4)
!        write(*,"(a,1p4e10.3)")"  - (comm time/inner p)  :", time(5)
!        write(*,"(a,1p4e10.3)")"  - (comm time/spmv)     :", time(6)
!      endif
!    endif
!
!    if(monoPRM%show_time_statistics)then
!      if(monoCOM%my_rank == 0) write(*,"(a,1p4e10.3)")" ** monolis solution time statistics"
!      if(monoCOM%my_rank == 0) write(*,"(a,1p4e10.3)")"                           max       min       average   std"
!
!      time(1) = monoPRM%tprep
!      call monolis_time_statistics (monoCOM, time(1), t_max, t_min, t_avg, t_sd)
!      if(monoCOM%my_rank == 0) write(*,"(a,1p4e10.3)")"  - solution/prepost time:", t_max, t_min, t_avg, t_sd
!
!      time(2) = monoPRM%tspmv
!      call monolis_time_statistics (monoCOM, time(2), t_max, t_min, t_avg, t_sd)
!      if(monoCOM%my_rank == 0) write(*,"(a,1p4e10.3)")"  - solution/SpMV    time:", t_max, t_min, t_avg, t_sd
!
!      time(3) = monoPRM%tdotp
!      call monolis_time_statistics (monoCOM, time(3), t_max, t_min, t_avg, t_sd)
!      if(monoCOM%my_rank == 0) write(*,"(a,1p4e10.3)")"  - solution/inner p time:", t_max, t_min, t_avg, t_sd
!
!      time(4) = monoPRM%tprec
!      call monolis_time_statistics (monoCOM, time(4), t_max, t_min, t_avg, t_sd)
!      if(monoCOM%my_rank == 0) write(*,"(a,1p4e10.3)")"  - solution/precond time:", t_max, t_min, t_avg, t_sd
!
!      time(5) = monoPRM%tcomm_dotp
!      call monolis_time_statistics (monoCOM, time(5), t_max, t_min, t_avg, t_sd)
!      if(monoCOM%my_rank == 0) write(*,"(a,1p4e10.3)")"  - (comm time/inner p)  :", t_max, t_min, t_avg, t_sd
!
!      time(6) = monoPRM%tcomm_spmv
!      call monolis_time_statistics (monoCOM, time(6), t_max, t_min, t_avg, t_sd)
!      if(monoCOM%my_rank == 0) write(*,"(a,1p4e10.3)")"  - (comm time/spmv)     :", t_max, t_min, t_avg, t_sd
!    endif
!
!    !> get average time
!    time(1) = monoPRM%tprep
!    time(2) = monoPRM%tspmv
!    time(3) = monoPRM%tdotp
!    time(4) = monoPRM%tprec
!    time(5) = monoPRM%tcomm_dotp
!    time(6) = monoPRM%tcomm_spmv
!
!    call monolis_allreduce_R(6, time, monolis_mpi_sum, monoCOM%comm)
!    time = time/dble(monolis_mpi_global_comm_size())
!
!    monoPRM%tprep = time(1)
!    monoPRM%tspmv = time(2)
!    monoPRM%tdotp = time(3)
!    monoPRM%tprec = time(4)
!    monoPRM%tcomm_dotp = time(5)
!    monoPRM%tcomm_spmv = time(6)
  end subroutine monolis_timer_finalize

end module mod_monolis_solve
