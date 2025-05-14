!> COCG 法モジュール
module mod_monolis_solver_COCG
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
  !> COCG 法
  subroutine monolis_solver_COCG(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in,out] 行列構造体
    type(monolis_mat), target, intent(inout) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC
    integer(kint) :: N, NP, NDOF
    integer(kint) :: iter, iter_RR
    real(kdouble) :: tspmv, tdotp, tcomm_spmv, tcomm_dotp
    complex(kdouble) :: alpha, beta, rho, rho1, omega, B2
    complex(kdouble), allocatable :: R(:), Z(:), Q(:), P(:)
    complex(kdouble), pointer :: B(:), X(:)
    logical :: is_converge

#ifdef DEBUG
    call monolis_std_debug_log_header("monolis_solver_COCG")
#endif

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    X => monoMAT%C%X
    B => monoMAT%C%B
    iter_RR = 200

    tspmv = monoPRM%Rarray(monolis_R_time_spmv)
    tcomm_spmv = monoPRM%Rarray(monolis_R_time_comm_spmv)
    tdotp = monoPRM%Rarray(monolis_R_time_dotp)
    tcomm_dotp = monoPRM%Rarray(monolis_R_time_comm_dotp)

    if(monoPRM%Iarray(monolis_prm_I_is_init_x) == monolis_I_true)then
      X = 0.0d0
    endif

    call monolis_alloc_C_1d(R, NDOF*NP)
    call monolis_alloc_C_1d(Z, NDOF*NP)
    call monolis_alloc_C_1d(Q, NDOF*NP)
    call monolis_alloc_C_1d(P, NDOF*NP)

    call monolis_residual_main_C(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
    call monolis_set_converge_C(monoCOM, monoMAT, R, B2, is_converge, tdotp, tcomm_dotp)
    if(is_converge) return

    do iter = 1, monoPRM%Iarray(monolis_prm_I_max_iter)
      call monolis_precond_apply_C(monoPRM, monoCOM, monoMAT, monoPREC, R, Z)
      call monolis_inner_product_main_C(monoCOM, N*NDOF, R, Z, rho, tdotp, tcomm_dotp)

      if(1 < iter)then
        beta = rho/rho1
        call monolis_vec_AXPBY_C(N*NDOF, beta, P, (1.0d0, 0.0d0), Z, P)
      else
        call monolis_vec_copy_C(N*NDOF, Z, P)
      endif

      call monolis_matvec_product_main_C(monoCOM, monoMAT, P, Q, tspmv, tcomm_spmv)
      call monolis_inner_product_main_C(monoCOM, N*NDOF, P, Q, omega, tdotp, tcomm_dotp)
      alpha = rho/omega

      call monolis_vec_AXPBY_C(N*NDOF, alpha, P, (1.0d0, 0.0d0), X, X)

      if(mod(iter, iter_RR) == 0)then
        call monolis_residual_main_C(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
      else
        call monolis_vec_AXPBY_C(N*NDOF, -alpha, Q, (1.0d0, 0.0d0), R, R)
      endif

      call monolis_check_converge_C(monoPRM, monoCOM, monoMAT, R, B2, iter, is_converge, tdotp, tcomm_dotp)
      if(is_converge) exit

      rho1 = rho
    enddo

    call monolis_mpi_update_C_wrapper(monoCOM, NDOF, monoMAT%n_dof_index, X, tcomm_spmv)

    monoPRM%Rarray(monolis_R_time_spmv) = tspmv
    monoPRM%Rarray(monolis_R_time_comm_spmv) = tcomm_spmv
    monoPRM%Rarray(monolis_R_time_dotp) = tdotp
    monoPRM%Rarray(monolis_R_time_comm_dotp) = tcomm_dotp

    call monolis_dealloc_C_1d(R)
    call monolis_dealloc_C_1d(Z)
    call monolis_dealloc_C_1d(Q)
    call monolis_dealloc_C_1d(P)
  end subroutine monolis_solver_COCG

end module mod_monolis_solver_COCG
