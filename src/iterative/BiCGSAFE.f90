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
    integer(kint) :: iter, iter_RR
    real(kdouble) :: alpha, beta, zeta, eta, B2, r1, r2, rho
    real(kdouble) :: y_y, y_r, y_v, v_v, v_r, denom, C(5)
    real(kdouble) :: tspmv, tdotp, tcomm_spmv, tcomm_dotp, ths
    logical :: is_converge
    real(kdouble), allocatable :: R(:), P(:), Z(:), U(:), Y(:), R0(:), V(:), T1(:), T2(:)
    real(kdouble), allocatable :: AP(:), MR(:), AU(:)
    real(kdouble), pointer :: B(:), X(:)

    call monolis_std_debug_log_header("monolis_solver_BiCGSAFE")

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NNDOF = N*NDOF
    X => monoMAT%R%X
    B => monoMAT%R%B
    iter_RR = 200
    ths   = 1.0d-15

    tspmv = monoPRM%Rarray(monolis_R_time_spmv)
    tcomm_spmv = monoPRM%Rarray(monolis_R_time_comm_spmv)
    tdotp = monoPRM%Rarray(monolis_R_time_dotp)
    tcomm_dotp = monoPRM%Rarray(monolis_R_time_comm_dotp)

    if(monoPRM%Iarray(monolis_prm_I_is_init_x) == monolis_I_true)then
      X = 0.0d0
    endif

    call monolis_alloc_R_1d(P, NDOF*NP)
    call monolis_alloc_R_1d(U, NDOF*NP)
    call monolis_alloc_R_1d(Z, NDOF*NP)
    call monolis_alloc_R_1d(Y, NDOF*NP)
    call monolis_alloc_R_1d(R, NDOF*NP)
    call monolis_alloc_R_1d(V, NDOF*NP)
    call monolis_alloc_R_1d(T1, NDOF*NP)
    call monolis_alloc_R_1d(T2, NDOF*NP)
    call monolis_alloc_R_1d(AP, NDOF*NP)
    call monolis_alloc_R_1d(MR, NDOF*NP)
    call monolis_alloc_R_1d(AU, NDOF*NP)
    call monolis_alloc_R_1d(R0, NDOF*NP)

    call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
    call monolis_set_converge_R(monoCOM, monoMAT, R, B2, is_converge, tdotp, tcomm_dotp)
    if(is_converge) return

    call monolis_vec_copy_R(N*NDOF, R, R0)
    call monolis_inner_product_main_R(monoCOM, N*NDOF, R0, R, r1, tdotp, tcomm_dotp)

    beta = 0.0d0

    do iter = 1, monoPRM%Iarray(monolis_prm_I_max_iter)
      call monolis_precond_apply_R(monoPRM, monoCOM, monoMAT, monoPREC, R, MR)

      call monolis_vec_AXPBY_R(N*NDOF, 1.0d0, P, -1.0d0, U, T1)
      call monolis_vec_AXPBY_R(N*NDOF, beta, T1, 1.0d0, MR, P)

      call monolis_matvec_product_main_R(monoCOM, monoMAT, MR, V, tspmv, tcomm_spmv)
      call monolis_vec_AXPBY_R(N*NDOF, 1.0d0, AP, -1.0d0, AU, T1)
      call monolis_vec_AXPBY_R(N*NDOF, beta, T1, 1.0d0, V, AP)

      call monolis_inner_product_main_R(monoCOM, N*NDOF, R0, AP, r2, tdotp, tcomm_dotp)
      alpha = r1 / r2

      if (iter == 1) then
        call monolis_inner_product_main_R(monoCOM, N*NDOF, V, R, v_r, tdotp, tcomm_dotp)
        call monolis_inner_product_main_R(monoCOM, N*NDOF, V, V, v_v, tdotp, tcomm_dotp)
        
        if (abs(v_v) == 0.0d0) then
          zeta = 1.0d0
        else
          zeta = v_r / v_v
        endif
        eta = 0.0d0
      else
        call monolis_inner_product_main_R_no_comm(N*NDOF, Y, Y, C(1))
        call monolis_inner_product_main_R_no_comm(N*NDOF, Y, R, C(2))
        call monolis_inner_product_main_R_no_comm(N*NDOF, Y, V, C(3))
        call monolis_inner_product_main_R_no_comm(N*NDOF, V, V, C(4))
        call monolis_inner_product_main_R_no_comm(N*NDOF, V, R, C(5))
        call monolis_allreduce_R(5, C, monolis_mpi_sum, monoCOM%comm)
        y_y = C(1)
        y_r = C(2)
        y_v = C(3)
        v_v = C(4)
        v_r = C(5)

        denom = v_v * y_y - y_v * y_v

        if (abs(denom) == 0.0d0) then
          zeta = 1.0d0
          eta = 0.0d0
        else
          zeta = (y_y * v_r - y_r * y_v) / denom
          eta  = (v_v * y_r - y_v * v_r) / denom
        endif
      endif

      call monolis_vec_AXPBY_R(N*NDOF, zeta, AP, eta, Y, T1)
      call monolis_precond_apply_R(monoPRM, monoCOM, monoMAT, monoPREC, T1, T2)
      call monolis_vec_AXPBY_R(N*NDOF, eta*beta, U, 1.0d0, T2, U)

      call monolis_vec_AXPBY_R(N*NDOF, eta, Z, -alpha, U, T1)
      call monolis_vec_AXPBY_R(N*NDOF, zeta, MR, 1.0d0, T1, Z)

      call monolis_matvec_product_main_R(monoCOM, monoMAT, U, AU, tspmv, tcomm_spmv)
      call monolis_vec_AXPBY_R(N*NDOF, eta, Y, -alpha, AU, T1)
      call monolis_vec_AXPBY_R(N*NDOF, zeta, V, 1.0d0, T1, Y)

      call monolis_vec_AXPBY_R(N*NDOF, alpha, P, 1.0d0, Z, T1)
      call monolis_vec_AXPBY_R(N*NDOF, 1.0d0, T1, 1.0d0, X, X)

      if(mod(iter, iter_RR) == 0)then
        call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
      else
        call monolis_vec_AXPBY_R(N*NDOF, -alpha, AP, -1.0d0, Y, T1)
        call monolis_vec_AXPBY_R(N*NDOF, 1.0d0, T1, 1.0d0, R, R)
      endif

      call monolis_check_converge_R(monoPRM, monoCOM, monoMAT, R, B2, iter, is_converge, tdotp, tcomm_dotp)
      if(is_converge) exit

      rho = r1
      call monolis_inner_product_main_R(monoCOM, N*NDOF, R0, R, r1, tdotp, tcomm_dotp)
      
      if (abs(zeta) == 0.0d0 .or. abs(rho) == 0.0d0) then
        beta = 0.0d0
      else
        beta = (alpha/zeta) * (r1/rho)
      endif
    enddo

    call monolis_mpi_update_R(monoCOM, NDOF, X, tcomm_spmv)

    monoPRM%Rarray(monolis_R_time_spmv) = tspmv
    monoPRM%Rarray(monolis_R_time_comm_spmv) = tcomm_spmv
    monoPRM%Rarray(monolis_R_time_dotp) = tdotp
    monoPRM%Rarray(monolis_R_time_comm_dotp) = tcomm_dotp

    call monolis_dealloc_R_1d(R)
    call monolis_dealloc_R_1d(Z)
    call monolis_dealloc_R_1d(P)
    call monolis_dealloc_R_1d(U)
    call monolis_dealloc_R_1d(Y)
    call monolis_dealloc_R_1d(V)
    call monolis_dealloc_R_1d(T1)
    call monolis_dealloc_R_1d(T2)
    call monolis_dealloc_R_1d(AP)
    call monolis_dealloc_R_1d(MR)
    call monolis_dealloc_R_1d(AU)
    call monolis_dealloc_R_1d(R0)
  end subroutine monolis_solver_BiCGSAFE

end module mod_monolis_solver_BiCGSAFE
