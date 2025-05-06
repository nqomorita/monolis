!> BiCGSAFE 法モジュール
module mod_monolis_solver_BiCGSAFE
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
  !> BiCGSAFE 法
  subroutine monolis_solver_BiCGSAFE(monoPRM, monoCOM, monoMAT, monoPREC)
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
    real(kdouble) :: alpha, beta, zeta, eta, temp, B2, r1, r2
    real(kdouble) :: v_r, y_r, v_v, y_y, v_y, denom
    real(kdouble) :: tspmv, tdotp, tcomm_spmv, tcomm_dotp
    logical :: is_converge
    real(kdouble), allocatable :: R(:), Z(:), P(:), AP(:), U(:), Y(:), R0(:), V(:), T(:)
    real(kdouble), pointer :: B(:), X(:)

    call monolis_std_debug_log_header("monolis_solver_BiCGSAFE")

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NNDOF = N*NDOF
    X => monoMAT%R%X
    B => monoMAT%R%B
    iter_RR = 200

    tspmv = monoPRM%Rarray(monolis_R_time_spmv)
    tcomm_spmv = monoPRM%Rarray(monolis_R_time_comm_spmv)
    tdotp = monoPRM%Rarray(monolis_R_time_dotp)
    tcomm_dotp = monoPRM%Rarray(monolis_R_time_comm_dotp)

    if(monoPRM%Iarray(monolis_prm_I_is_init_x) == monolis_I_true)then
      X = 0.0d0
    endif

    call monolis_alloc_R_1d(R, NDOF*NP)
    call monolis_alloc_R_1d(Z, NDOF*NP)
    call monolis_alloc_R_1d(P, NDOF*NP)
    call monolis_alloc_R_1d(AP, NDOF*NP)
    call monolis_alloc_R_1d(U, NDOF*NP)
    call monolis_alloc_R_1d(Y, NDOF*NP)
    call monolis_alloc_R_1d(R0, NDOF*NP)
    call monolis_alloc_R_1d(V, NDOF*NP)
    call monolis_alloc_R_1d(T, NDOF*NP)

    call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
    call monolis_set_converge_R(monoCOM, monoMAT, R, B2, is_converge, tdotp, tcomm_dotp)
    if(is_converge) return

    call monolis_vec_copy_R(N, NDOF, R, R0)

    beta = 0.0d0

    do iter = 1, monoPRM%Iarray(monolis_prm_I_max_iter)
      call monolis_precond_apply_R(monoPRM, monoCOM, monoMAT, monoPREC, R, Z)

      call monolis_vec_AXPBY_R(N, NDOF, -1.0d0, U, 1.0d0, P, T)
      call monolis_vec_AXPBY_R(N, NDOF, beta, T, 1.0d0, Z, P)

!      call monolis_matvec_product_main_R(monoCOM, monoMAT, P, AP, tspmv, tcomm_spmv)

      call monolis_inner_product_main_R(monoCOM, N, NDOF, R0, R, r1, tdotp, tcomm_dotp)
      call monolis_inner_product_main_R(monoCOM, N, NDOF, R0, AP, r2, tdotp, tcomm_dotp)
      alpha = r1 / r2

      call monolis_matvec_product_main_R(monoCOM, monoMAT, Z, V, tspmv, tcomm_spmv)

      if (iter == 1) then
        call monolis_inner_product_main_R(monoCOM, N, NDOF, V, R, v_r, tdotp, tcomm_dotp)
        call monolis_inner_product_main_R(monoCOM, N, NDOF, V, V, v_v, tdotp, tcomm_dotp)
        
        if (abs(v_v) < 1.0d-16) then
          zeta = 1.0d0
        else
          zeta = v_r / v_v
        endif
        eta = 0.0d0
      else
        call monolis_inner_product_main_R(monoCOM, N, NDOF, Y, Y, y_y, tdotp, tcomm_dotp)
        call monolis_inner_product_main_R(monoCOM, N, NDOF, V, R, v_r, tdotp, tcomm_dotp)
        call monolis_inner_product_main_R(monoCOM, N, NDOF, Y, R, y_r, tdotp, tcomm_dotp)
        call monolis_inner_product_main_R(monoCOM, N, NDOF, V, V, v_v, tdotp, tcomm_dotp)
        call monolis_inner_product_main_R(monoCOM, N, NDOF, Y, V, y_v, tdotp, tcomm_dotp)
        call monolis_inner_product_main_R(monoCOM, N, NDOF, V, Y, v_y, tdotp, tcomm_dotp)
        
        denom = v_v * y_y - v_y * y_v
        
        if (abs(denom) < 1.0d-16) then
          zeta = 1.0d0
          eta = 0.0d0
        else
          zeta = (y_y * v_r - y_r * v_y) / denom
          eta = (v_v * y_r - y_v * v_r) / denom
        endif
      endif

      call monolis_vec_AXPBY_R(N, NDOF, zeta, AP, eta, Y, Y)
      call monolis_precond_apply_R(monoPRM, monoCOM, monoMAT, monoPREC, Y, U)
      if (iter > 1) then
        call monolis_vec_AXPBY_R(N, NDOF, eta*beta, U, 1.0d0, U, U)
      endif

      call monolis_vec_AXPBY_R(N, NDOF, zeta, Z, eta, Z, Z)
      call monolis_vec_AXPBY_R(N, NDOF, -alpha, U, 1.0d0, Z, Z)

      call monolis_matvec_product_main_R(monoCOM, monoMAT, U, Y, tspmv, tcomm_spmv)
      call monolis_vec_AXPBY_R(N, NDOF, zeta, V, eta, Y, Y)
      call monolis_vec_AXPBY_R(N, NDOF, -alpha, Y, 1.0d0, Y, Y)

      call monolis_vec_AXPBY_R(N, NDOF, alpha, P, 1.0d0, X, X)
      call monolis_vec_AXPBY_R(N, NDOF, 1.0d0, Z, 1.0d0, X, X)

      if(mod(iter, iter_RR) == 0)then
        call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
      else
        call monolis_vec_AXPBY_R(N, NDOF, -alpha, AP, 1.0d0, R, R)
        call monolis_vec_AXPBY_R(N, NDOF, -1.0d0, Y, 1.0d0, R, R)
      endif

      temp = r0_r
      call monolis_inner_product_main_R(monoCOM, N, NDOF, R0, R, r0_r, tdotp, tcomm_dotp)
      
      if (abs(temp) < 1.0d-16 .or. abs(zeta) < 1.0d-16) then
        beta = 0.0d0
      else
        beta = (alpha/zeta) * (r0_r/temp)
      endif

      call monolis_check_converge_R(monoPRM, monoCOM, monoMAT, R, B2, iter, is_converge, tdotp, tcomm_dotp)
      if(is_converge) exit
    enddo

    call monolis_mpi_update_R(monoCOM, NDOF, X, tcomm_spmv)

    monoPRM%Rarray(monolis_R_time_spmv) = tspmv
    monoPRM%Rarray(monolis_R_time_comm_spmv) = tcomm_spmv
    monoPRM%Rarray(monolis_R_time_dotp) = tdotp
    monoPRM%Rarray(monolis_R_time_comm_dotp) = tcomm_dotp

    call monolis_dealloc_R_1d(R)
    call monolis_dealloc_R_1d(Z)
    call monolis_dealloc_R_1d(P)
    call monolis_dealloc_R_1d(AP)
    call monolis_dealloc_R_1d(U)
    call monolis_dealloc_R_1d(Y)
    call monolis_dealloc_R_1d(R0)
    call monolis_dealloc_R_1d(V)
    call monolis_dealloc_R_1d(T)
  end subroutine monolis_solver_BiCGSAFE

end module mod_monolis_solver_BiCGSAFE
