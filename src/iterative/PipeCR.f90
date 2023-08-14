!> PipeCR 法モジュール
module mod_monolis_solver_PipeCR
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
  !> PipeCR 法
  subroutine monolis_solver_PipeCR(monoPRM, monoCOM, monoMAT, monoPREC)
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
    real(kdouble) :: tol, resid, R2, B2, U2
    real(kdouble) :: alpha, alpha1, beta, gamma, gamma1, delta, phi, utol
    real(kdouble) :: buf(3), CG(3)
    real(kdouble) :: tspmv, tdotp, tcomm_spmv, tcomm_dotp
    real(kdouble), allocatable :: R(:), U(:), V(:), Q(:), P(:), Z(:), L(:), M(:), S(:)
    real(kdouble), pointer :: B(:), X(:)
    logical :: is_converge

    call monolis_std_debug_log_header("monolis_solver_PipeCR")

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NNDOF = N*NDOF
    X => monoMAT%R%X
    B => monoMAT%R%B
    iter_RR = 50

    tspmv = monoPRM%Rarray(monolis_R_time_spmv)
    tcomm_spmv = monoPRM%Rarray(monolis_R_time_comm_spmv)
    tdotp = monoPRM%Rarray(monolis_R_time_dotp)
    tcomm_dotp = monoPRM%Rarray(monolis_R_time_comm_dotp)

    if(monoPRM%Iarray(monolis_prm_I_is_init_x) == monolis_I_true)then
      X = 0.0d0
    endif

    tol = monoPRM%Rarray(monolis_prm_R_tol)

    call monolis_alloc_R_1d(R, NDOF*NP)
    call monolis_alloc_R_1d(U, NDOF*NP)
    call monolis_alloc_R_1d(V, NDOF*NP)
    call monolis_alloc_R_1d(Q, NDOF*NP)
    call monolis_alloc_R_1d(P, NDOF*NP)
    call monolis_alloc_R_1d(Z, NDOF*NP)
    call monolis_alloc_R_1d(L, NDOF*NP)
    call monolis_alloc_R_1d(M, NDOF*NP)
    call monolis_alloc_R_1d(S, NDOF*NP)

    call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
    call monolis_inner_product_main_R(monoCOM, N, NDOF, R, R, B2, tdotp, tcomm_dotp)
    call monolis_inner_product_main_R(monoCOM, N, NDOF, R, R, R2, tdotp, tcomm_dotp)
    call monolis_precond_apply_R(monoPRM, monoCOM, monoMAT, monoPREC, R, U)
    call monolis_inner_product_main_R(monoCOM, N, NDOF, U, U, U2, tdotp, tcomm_dotp)
    call monolis_matvec_product_main_R(monoCOM, monoMAT, U, V, tspmv, tcomm_spmv)

    phi  = dsqrt(R2/U2)
    utol = tol/phi

    do iter = 1, monoPRM%Iarray(monolis_prm_I_max_iter)
      call monolis_precond_apply_R(monoPRM, monoCOM, monoMAT, monoPREC, V, M)

      call monolis_inner_product_main_R_no_comm(N, NDOF, V, U, CG(1))
      call monolis_inner_product_main_R_no_comm(N, NDOF, V, M, CG(2))
      call monolis_inner_product_main_R_no_comm(N, NDOF, U, U, CG(3))
      call monolis_allreduce_R(3, CG, monolis_mpi_sum, monoCOM%comm)

      call monolis_matvec_product_main_R(monoCOM, monoMAT, M, L, tspmv, tcomm_spmv)

      gamma = CG(1)
      delta = CG(2)
      U2    = CG(3)

      resid = dsqrt(U2/B2)
      !if(monoCOM%my_rank == 0 .and. monoPRM%show_iterlog) write (*,"(i7, 1pe16.6)") iter, resid*phi
      if(resid <= utol) exit

      if(1 < iter)then
        beta  = gamma*gamma1
        alpha = gamma/(delta-beta*gamma*alpha1)
      else
        beta  = 0.0d0
        alpha = gamma/delta
      endif

      call monolis_vec_AXPBY_R(N, NDOF, beta, Z, 1.0d0, L, Z)
      call monolis_vec_AXPBY_R(N, NDOF, beta, Q, 1.0d0, M, Q)
      call monolis_vec_AXPBY_R(N, NDOF, beta, P, 1.0d0, U, P)
      call monolis_vec_AXPBY_R(N, NDOF,alpha, P, 1.0d0, X, X)

      if(mod(iter, iter_RR) == 0)then
        call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
        call monolis_precond_apply_R(monoPRM, monoCOM, monoMAT, monoPREC, R, U)
        call monolis_matvec_product_main_R(monoCOM, monoMAT, U, V, tspmv, tcomm_spmv)
      else
        call monolis_vec_AXPBY_R(N, NDOF,-alpha, Q, 1.0d0, U, U)
        call monolis_vec_AXPBY_R(N, NDOF,-alpha, Z, 1.0d0, V, V)
      endif

      gamma1 = 1.0d0/gamma
      alpha1 = 1.0d0/alpha
    enddo

    call monolis_mpi_update_R(monoCOM, NDOF, X, tcomm_spmv)

    monoPRM%Rarray(monolis_R_time_spmv) = tspmv
    monoPRM%Rarray(monolis_R_time_comm_spmv) = tcomm_spmv
    monoPRM%Rarray(monolis_R_time_dotp) = tdotp
    monoPRM%Rarray(monolis_R_time_comm_dotp) = tcomm_dotp

    call monolis_dealloc_R_1d(R)
    call monolis_dealloc_R_1d(U)
    call monolis_dealloc_R_1d(V)
    call monolis_dealloc_R_1d(Q)
    call monolis_dealloc_R_1d(P)
    call monolis_dealloc_R_1d(Z)
    call monolis_dealloc_R_1d(L)
    call monolis_dealloc_R_1d(M)
    call monolis_dealloc_R_1d(S)
  end subroutine monolis_solver_PipeCR

end module mod_monolis_solver_PipeCR
