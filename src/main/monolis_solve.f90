module mod_monolis_solve
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_scaling
  use mod_monolis_precond
  use mod_monolis_reorder
  use mod_monolis_util
  use mod_monolis_util_debug
  use mod_monolis_util_time
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
  use mod_monolis_solver_GMRES
  implicit none

contains

  subroutine monolis_solve(monolis, B, X)
    implicit none
    type(monolis_structure) :: monolis
    real(kdouble) :: B(:), X(:)
    integer(kint) :: i

    call monolis_set_RHS(monolis%MAT, B)
    call monolis_set_initial_solution(monolis%MAT, X)
    call monolis_set_initial_comm(monolis%COM, monolis%MAT)

    call monolis_solve_(monolis%PRM, monolis%COM, monolis%MAT)

    call monolis_get_solution(monolis%MAT, X)
  end subroutine monolis_solve

  subroutine monolis_set_RHS(monoMAT, B)
    implicit none
    type(monolis_mat) :: monoMAT
    real(kdouble) :: B(:)
    integer(kint) :: i

    do i = 1, monoMAT%NP*monoMAT%NDOF
      monoMAT%B(i) = B(i)
    enddo
  end subroutine monolis_set_RHS

  subroutine monolis_set_initial_solution(monoMAT, X)
    implicit none
    type(monolis_mat) :: monoMAT
    real(kdouble) :: X(:)
    integer(kint) :: i

    do i = 1, monoMAT%NP*monoMAT%NDOF
      monoMAT%X(i) = X(i)
    enddo
  end subroutine monolis_set_initial_solution

  subroutine monolis_set_initial_comm(monoCOM, monoMAT)
    implicit none
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    if(monolis_global_commsize() > 1)then
      monoMAT%N = monoCOM%internal_nnode
    endif
  end subroutine monolis_set_initial_comm

  subroutine monolis_get_solution(monoMAT, X)
    implicit none
    type(monolis_mat) :: monoMAT
    real(kdouble) :: X(:)
    integer(kint) :: i

    do i = 1, monoMAT%NP*monoMAT%NDOF
      X(i) = monoMAT%X(i)
    enddo
  end subroutine monolis_get_solution

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

    call monolis_timer_initialize(monoPRM, monoCOM)
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
        call monolis_timer_initialize(monoPRM, monoCOM)
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

    if(monoPRM%is_debug) call monolis_debug_header("monolis_solver v0.0.0")

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

      case (monolis_iter_GMRES)
        call monolis_solver_GMRES(monoPRM, monoCOM, monoMAT)

    end select
  end subroutine monolis_solver

  !> c interface
  subroutine monolis_solve_c(N, NP, NZ, NDOF, A, X, B, index, item, &
    myrank, comm, commsize, &
    recv_n_neib, recv_nitem, recv_neib_pe, recv_index, recv_item, &
    send_n_neib, send_nitem, send_neib_pe, send_index, send_item, &
    method, precond, maxiter, tol, &
    iterlog, timelog_statistics, timelog, summary, &
    is_check_diag, is_measurement, is_init_x, curiter, curresid, time) &
    & bind(c, name = "monolis_solve_c_main")
    implicit none
    type(monolis_structure) :: monolis
    integer(c_int), intent(in), value :: N, NP, NZ, NDOF
    integer(c_int), intent(in), value :: myrank, comm, commsize
    integer(c_int), intent(in), value :: recv_n_neib, send_n_neib, recv_nitem, send_nitem
    integer(c_int), intent(in), value :: method, precond, maxiter
    integer(c_int), intent(in), value :: iterlog, timelog, timelog_statistics, summary
    integer(c_int), intent(in), value :: is_check_diag, is_measurement, is_init_x
    integer(c_int), intent(in), target :: index(0:NP)
    integer(c_int), intent(in), target :: item(NZ)
    integer(c_int), intent(in), target :: recv_neib_pe(recv_n_neib)
    integer(c_int), intent(in), target :: recv_index(0:recv_n_neib), recv_item(recv_nitem)
    integer(c_int), intent(in), target :: send_neib_pe(send_n_neib)
    integer(c_int), intent(in), target :: send_index(0:send_n_neib), send_item(send_nitem)
    integer(c_int), intent(out), target :: curiter
    real(c_double), intent(in), value :: tol
    real(c_double), intent(in), target :: A(NDOF*NDOF*NZ)
    real(c_double), intent(in), target :: X(NDOF*NP)
    real(c_double), intent(in), target :: B(NDOF*NP)
    real(c_double), intent(out), target :: curresid
    real(c_double), intent(out), target :: time(7)

    !> for monoMAT
    monolis%MAT%N = N
    monolis%MAT%NP = NP
    monolis%MAT%NZ = NZ
    monolis%MAT%NDOF = NDOF
    monolis%MAT%A => A
    monolis%MAT%X => X
    monolis%MAT%B => B
    monolis%MAT%index => index
    monolis%MAT%item => item

    !> for monoCOM
    monoliS%COM%myrank = myrank
    monoliS%COM%comm = comm
    monoliS%COM%commsize = commsize
    monoliS%COM%recv_n_neib = recv_n_neib
    monoliS%COM%recv_neib_pe => recv_neib_pe
    monoliS%COM%recv_index => recv_index
    monoliS%COM%recv_item => recv_item
    monoliS%COM%send_n_neib = send_n_neib
    monoliS%COM%send_neib_pe => send_neib_pe
    monoliS%COM%send_index => send_index
    monoliS%COM%send_item => send_item

    !> for monoPRM
    monolis%PRM%method = method
    monolis%PRM%precond = precond
    monolis%PRM%maxiter = maxiter
    monolis%PRM%tol = tol

    monolis%PRM%is_scaling    = .false.
    monolis%PRM%is_reordering = .false.
    monolis%PRM%is_init_x     = .true.
    monolis%PRM%is_debug      = .false.
    monolis%PRM%is_check_diag = .true.
    monolis%PRM%is_measurement= .false.
    monolis%PRM%show_iterlog  = .false.
    monolis%PRM%show_time     = .false.
    monolis%PRM%show_summary  = .false.

    if(iterlog == 1) monolis%PRM%show_iterlog  = .true.
    if(timelog == 1) monolis%PRM%show_time     = .true.
    if(timelog_statistics == 1) monolis%PRM%show_time_statistics = .true.
    if(summary == 1) monolis%PRM%show_summary  = .true.
    if(is_check_diag == 0) monolis%PRM%is_check_diag= .false.
    if(is_init_x == 0) monolis%PRM%is_init_x= .false.
    if(is_measurement == 1) monolis%PRM%is_measurement= .true.

    call monolis_solve_(monolis%PRM, monolis%COM, monolis%MAT)

    curiter = monolis%PRM%curiter
    curresid = monolis%PRM%curresid
    time(1) = monolis%PRM%tsol
    time(2) = monolis%PRM%tprep
    time(3) = monolis%PRM%tspmv
    time(4) = monolis%PRM%tdotp
    time(5) = monolis%PRM%tprec
    time(6) = monolis%PRM%tcomm_dotp
    time(7) = monolis%PRM%tcomm_spmv
  end subroutine monolis_solve_c

  subroutine monolis_timer_finalize(monoPRM, monoCOM)
    use mod_monolis_linalg_com
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    real(kdouble) :: t1, time(6), t_max, t_min, t_avg, t_sd
    logical :: is_output

    call monolis_debug_header("monolis_timer_finalize")

    call monolis_barrier_(monoCOM%comm)
    t1 = monolis_get_time()
    monoPRM%tsol = t1 - monoPRM%tsol

    if(monoPRM%show_summary .and. monoCOM%myrank == 0)then
      write(*,"(a,i10)")" ** monolis converge iter:", monoPRM%curiter
      write(*,"(a,1p4e10.3)")" ** monolis rel. residual:", monoPRM%curresid
    endif

    is_output = monoPRM%show_summary .or. monoPRM%show_time
    if(is_output .and. monoCOM%myrank == 0)then
      write(*,"(a,1p4e10.3)")" ** monolis solution time:", monoPRM%tsol
    endif

    if(monoPRM%show_time)then
      time(1) = monoPRM%tprep
      time(2) = monoPRM%tspmv
      time(3) = monoPRM%tdotp
      time(4) = monoPRM%tprec
      time(5) = monoPRM%tcomm_dotp
      time(6) = monoPRM%tcomm_spmv
      call monolis_allreduce_R(6, time, monolis_sum, monoCOM%comm)
      time = time/dble(monolis_global_commsize())

      if(monoCOM%myrank == 0)then
        write(*,"(a,1p4e10.3)")"  - solution/prepost time:", time(1)
        write(*,"(a,1p4e10.3)")"  - solution/SpMV    time:", time(2)
        write(*,"(a,1p4e10.3)")"  - solution/inner p time:", time(3)
        write(*,"(a,1p4e10.3)")"  - solution/precond time:", time(4)
        write(*,"(a,1p4e10.3)")"  - (comm time/inner p)  :", time(5)
        write(*,"(a,1p4e10.3)")"  - (comm time/spmv)     :", time(6)
      endif
    endif

    if(monoPRM%show_time_statistics)then
      if(monoCOM%myrank == 0) write(*,"(a,1p4e10.3)")" ** monolis solution time statistics"
      if(monoCOM%myrank == 0) write(*,"(a,1p4e10.3)")"                           max       min       average   std"

      time(1) = monoPRM%tprep
      call monolis_time_statistics (monoCOM, time(1), t_max, t_min, t_avg, t_sd)
      if(monoCOM%myrank == 0) write(*,"(a,1p4e10.3)")"  - solution/prepost time:", t_max, t_min, t_avg, t_sd

      time(2) = monoPRM%tspmv
      call monolis_time_statistics (monoCOM, time(2), t_max, t_min, t_avg, t_sd)
      if(monoCOM%myrank == 0) write(*,"(a,1p4e10.3)")"  - solution/SpMV    time:", t_max, t_min, t_avg, t_sd

      time(3) = monoPRM%tdotp
      call monolis_time_statistics (monoCOM, time(3), t_max, t_min, t_avg, t_sd)
      if(monoCOM%myrank == 0) write(*,"(a,1p4e10.3)")"  - solution/inner p time:", t_max, t_min, t_avg, t_sd

      time(4) = monoPRM%tprec
      call monolis_time_statistics (monoCOM, time(4), t_max, t_min, t_avg, t_sd)
      if(monoCOM%myrank == 0) write(*,"(a,1p4e10.3)")"  - solution/precond time:", t_max, t_min, t_avg, t_sd

      time(5) = monoPRM%tcomm_dotp
      call monolis_time_statistics (monoCOM, time(5), t_max, t_min, t_avg, t_sd)
      if(monoCOM%myrank == 0) write(*,"(a,1p4e10.3)")"  - (comm time/inner p)  :", t_max, t_min, t_avg, t_sd

      time(6) = monoPRM%tcomm_spmv
      call monolis_time_statistics (monoCOM, time(6), t_max, t_min, t_avg, t_sd)
      if(monoCOM%myrank == 0) write(*,"(a,1p4e10.3)")"  - (comm time/spmv)     :", t_max, t_min, t_avg, t_sd
    endif

    !> get average time
    time(1) = monoPRM%tprep
    time(2) = monoPRM%tspmv
    time(3) = monoPRM%tdotp
    time(4) = monoPRM%tprec
    time(5) = monoPRM%tcomm_dotp
    time(6) = monoPRM%tcomm_spmv

    call monolis_allreduce_R(6, time, monolis_sum, monoCOM%comm)
    time = time/dble(monolis_global_commsize())

    monoPRM%tprep = time(1)
    monoPRM%tspmv = time(2)
    monoPRM%tdotp = time(3)
    monoPRM%tprec = time(4)
    monoPRM%tcomm_dotp = time(5)
    monoPRM%tcomm_spmv = time(6)
  end subroutine monolis_timer_finalize

  subroutine monolis_time_statistics (monoCOM, time, t_max, t_min, t_avg, t_sd)
    implicit none
    type(monolis_com) :: monoCOM
    real(kdouble), intent(in) :: time
    real(kdouble), intent(out) :: t_max, t_min, t_avg, t_sd
    real(kdouble) :: tmp
    integer(kint) :: np

    np = monolis_global_commsize()

    t_max = time
    call monolis_allreduce_R1(t_max, monolis_max, monoCOM%comm)

    t_min = time
    call monolis_allreduce_R1(t_min, monolis_min, monoCOM%comm)

    t_avg = time
    call monolis_allreduce_R1(t_avg, monolis_sum, monoCOM%comm)
    t_avg = t_avg / np

    tmp = time*time
    call monolis_allreduce_R1(tmp, monolis_sum, monoCOM%comm)
    tmp = tmp / np

    if(np == 1)then
      t_sd = 0.0d0
    else
      t_sd = dsqrt(tmp - t_avg*t_avg)
    endif
  end subroutine monolis_time_statistics
end module mod_monolis_solve
