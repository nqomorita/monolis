!> COCG 法モジュール
module mod_monolis_solver_COCG
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_precond
  use mod_monolis_matvec
  use mod_monolis_linalg
  use mod_monolis_converge

  implicit none

contains

  !> CG 法
  subroutine monolis_solver_COCG(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> パラメータ構造体
    type(monolis_prm) :: monoPRM
    !> 通信テーブル構造体
    type(monolis_com) :: monoCOM
    !> 行列構造体
    type(monolis_mat), target :: monoMAT
    !> 前処理構造体
    type(monolis_mat) :: monoPREC
    integer(kint) :: N, NP, NDOF, NNDOF
    integer(kint) :: i, iter, iter_RR
    real(kdouble) :: alpha, beta, rho, rho1, omega, B2
    real(kdouble), allocatable :: R(:), Z(:), Q(:), P(:)
    real(kdouble), pointer :: B(:), X(:)
    logical :: is_converge

#ifdef DEBUG
    call monolis_std_debug_log_header("monolis_solver_COCG")
#endif

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NNDOF = N*NDOF
    X => monoMAT%R%X
    B => monoMAT%R%B
    iter_RR = 200

    !if(monoPRM%is_init_x) X = 0.0d0

    call monolis_alloc_R_1d(R, NDOF*NP)
    call monolis_alloc_R_1d(Z, NDOF*NP)
    call monolis_alloc_R_1d(Q, NDOF*NP)
    call monolis_alloc_R_1d(P, NDOF*NP)

!    call monolis_residual(monoCOM, monoMAT, X, B, R, monoPRM%tspmv, monoPRM%tcomm_spmv)
!    call monolis_set_converge(monoPRM, monoCOM, monoMAT, R, B2, is_converge, monoPRM%tdotp, monoPRM%tcomm_dotp)
!    if(is_converge) return
!
!    do iter = 1, monoPRM%maxiter
!      call monolis_precond_apply(monoPRM, monoCOM, monoMAT, R, Z)
!      call monolis_inner_product_R(monoCOM, N, NDOF, R, Z, rho, monoPRM%tdotp, monoPRM%tcomm_dotp)
!
!      if(1 < iter)then
!        beta = rho/rho1
!        call monolis_vec_AXPY(N, NDOF, beta, P, Z, P)
!      else
!        call monolis_vec_copy_R(N, NDOF, Z, P)
!      endif
!
!      call monolis_matvec(monoCOM, monoMAT, P, Q, monoPRM%tspmv, monoPRM%tcomm_spmv)
!      call monolis_inner_product_R(monoCOM, N, NDOF, P, Q, omega, monoPRM%tdotp, monoPRM%tcomm_dotp)
!      alpha = rho/omega
!
!      call monolis_vec_AXPY(N, NDOF, alpha, P, X, X)
!
!      if(mod(iter, iter_RR) == 0)then
!        call monolis_residual(monoCOM, monoMAT, X, B, R, monoPRM%tspmv, monoPRM%tcomm_spmv)
!      else
!        call monolis_vec_AXPY(N, NDOF, -alpha, Q, R, R)
!      endif
!
!      call monolis_check_converge(monoPRM, monoCOM, monoMAT, R, B2, iter, is_converge, monoPRM%tdotp, monoPRM%tcomm_dotp)
!      if(is_converge) exit
!
!      rho1 = rho
!    enddo
!
!    call monolis_update_R(monoCOM, NDOF, X, monoPRM%tcomm_spmv)

    call monolis_dealloc_R_1d(R)
    call monolis_dealloc_R_1d(Z)
    call monolis_dealloc_R_1d(Q)
    call monolis_dealloc_R_1d(P)
  end subroutine monolis_solver_COCG

end module mod_monolis_solver_COCG
