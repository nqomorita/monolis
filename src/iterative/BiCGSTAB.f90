!> BiCGStab 法モジュール
module mod_monolis_solver_BiCGSTAB
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_precond
  use mod_monolis_matvec
  use mod_monolis_inner_product
  use mod_monolis_vec_util
  use mod_monolis_converge

  implicit none

contains

  !> @ingroup solver
  !> BiCGStab 法
  subroutine monolis_solver_BiCGSTAB(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in,out] 行列構造体
    type(monolis_mat), target, intent(inout) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC
    integer(kint) :: N, NP, NDOF, NNDOF
    integer(kint) :: i, iter, iter_RR
    real(kdouble) :: alpha, beta, rho, rho1, c2, omega
    real(kdouble) :: B2, CG(2)
    real(kdouble) :: tspmv, tdotp, tcomm_spmv, tcomm_dotp
    real(kdouble), allocatable :: R(:), RT(:), P(:), PT(:), S(:), ST(:), T(:), V(:)
    real(kdouble), pointer :: B(:), X(:)
    logical :: is_converge

    call monolis_std_debug_log_header("monolis_solver_BiCGSTAB")

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NNDOF = N*NDOF
    X => monoMAT%R%X
    B => monoMAT%R%B
    iter_RR = 200
    omega = 0.0d0

    tspmv = monoPRM%Rarray(monolis_R_time_spmv)
    tcomm_spmv = monoPRM%Rarray(monolis_R_time_comm_spmv)
    tdotp = monoPRM%Rarray(monolis_R_time_dotp)
    tcomm_dotp = monoPRM%Rarray(monolis_R_time_comm_dotp)

    if(monoPRM%Iarray(monolis_prm_I_is_init_x) == monolis_I_true)then
      X = 0.0d0
    endif

    call monolis_alloc_R_1d(R , NDOF*NP)
    call monolis_alloc_R_1d(RT, NDOF*NP)
    call monolis_alloc_R_1d(P , NDOF*NP)
    call monolis_alloc_R_1d(PT, NDOF*NP)
    call monolis_alloc_R_1d(S , NDOF*NP)
    call monolis_alloc_R_1d(ST, NDOF*NP)
    call monolis_alloc_R_1d(T , NDOF*NP)
    call monolis_alloc_R_1d(V , NDOF*NP)

    call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
    call monolis_set_converge_R(monoCOM, monoMAT, R, B2, is_converge, tdotp, tcomm_dotp)
    if(is_converge) return

    call monolis_vec_copy_R(N, NDOF, R, RT)

    do iter = 1, monoPRM%Iarray(monolis_prm_I_max_iter)
      call monolis_inner_product_main_R(monoCOM, N, NDOF, R, RT, rho, tdotp, tcomm_dotp)

      if(1 < iter)then
        beta = (rho/rho1) * (alpha/omega)
        do i = 1, NNDOF
          P(i) = R(i) + beta * (P(i) - omega * V(i))
        enddo
      else
        call monolis_vec_copy_R(N, NDOF, R, P)
      endif

      call monolis_precond_apply_R(monoPRM, monoCOM, monoMAT, monoPREC, P, PT)
      call monolis_matvec_product_main_R(monoCOM, monoMAT, PT, V, tspmv, tcomm_spmv)
      call monolis_inner_product_main_R(monoCOM, N, NDOF, RT, V, c2, tdotp, tcomm_dotp)

      alpha = rho / c2
      call monolis_vec_AXPBY_R(N, NDOF, -alpha, V, 1.0d0, R, S)

      call monolis_precond_apply_R(monoPRM, monoCOM, monoMAT, monoPREC, S, ST)
      call monolis_matvec_product_main_R(monoCOM, monoMAT, ST, T, tspmv, tcomm_spmv)

      call monolis_inner_product_main_R_no_comm(N, NDOF, T, S, CG(1))
      call monolis_inner_product_main_R_no_comm(N, NDOF, T, T, CG(2))
      call monolis_allreduce_R(2, CG, monolis_mpi_sum, monoCOM%comm)

      if(CG(2) /= 0.0d0)then
        omega = CG(1) / CG(2)
      else
        omega = 0.0d0
      endif

      do i = 1, NNDOF
        X(i) = X(i) + alpha*PT(i) + omega*ST(i)
      enddo

      if(mod(iter, iter_RR) == 0)then
        call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
      else
        call monolis_vec_AXPBY_R(N, NDOF, -omega, T, 1.0d0, S, R)
      endif

      call monolis_check_converge_R(monoPRM, monoCOM, monoMAT, R, B2, iter, is_converge, tdotp, tcomm_dotp)
      if(is_converge) exit

      rho1 = rho
    enddo

    call monolis_mpi_update_R(monoCOM, NDOF, X, tcomm_spmv)

    monoPRM%Rarray(monolis_R_time_spmv) = tspmv
    monoPRM%Rarray(monolis_R_time_comm_spmv) = tcomm_spmv
    monoPRM%Rarray(monolis_R_time_dotp) = tdotp
    monoPRM%Rarray(monolis_R_time_comm_dotp) = tcomm_dotp

    call monolis_dealloc_R_1d(R )
    call monolis_dealloc_R_1d(RT)
    call monolis_dealloc_R_1d(P )
    call monolis_dealloc_R_1d(PT)
    call monolis_dealloc_R_1d(S )
    call monolis_dealloc_R_1d(ST)
    call monolis_dealloc_R_1d(T )
    call monolis_dealloc_R_1d(V )
  end subroutine monolis_solver_BiCGSTAB

end module mod_monolis_solver_BiCGSTAB
