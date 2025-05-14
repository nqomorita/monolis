!> IDR(s) 法モジュール
module mod_monolis_solver_IDRS
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_precond
  use mod_monolis_matvec
  use mod_monolis_inner_product
  use mod_monolis_vec_util
  use mod_monolis_converge
  use mod_monolis_lapack

  implicit none

contains

  !> @ingroup solver
  !> IDRS 法
  subroutine monolis_solver_IDRS(monoPRM, monoCOM, monoMAT, monoPREC)
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
    integer(kint) :: i, j, iter, iter_RR, S
    real(kdouble) :: B2, norm, beta, a1, a2
    real(kdouble) :: tspmv, tdotp, tcomm_spmv, tcomm_dotp
    logical :: is_converge
    real(kdouble), allocatable :: R(:), P(:,:), U(:,:), G(:,:)
    real(kdouble), allocatable :: E(:,:), F(:), Y(:), alpha(:), Z(:), V(:), T(:), C(:)
    real(kdouble), pointer :: B(:), X(:)

    call monolis_std_debug_log_header("monolis_solver_IDRS")

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NNDOF = N*NDOF
    X => monoMAT%R%X
    B => monoMAT%R%B
    iter_RR = 200
    S = monoPRM%Iarray(MONOLIS_PRM_I_IDRS_DIM)

    tspmv = monoPRM%Rarray(monolis_R_time_spmv)
    tcomm_spmv = monoPRM%Rarray(monolis_R_time_comm_spmv)
    tdotp = monoPRM%Rarray(monolis_R_time_dotp)
    tcomm_dotp = monoPRM%Rarray(monolis_R_time_comm_dotp)

    if(monoPRM%Iarray(monolis_prm_I_is_init_x) == monolis_I_true)then
      X = 0.0d0
    endif

    call monolis_alloc_R_1d(R, NDOF*NP)
    call monolis_alloc_R_1d(Z, NDOF*NP)
    call monolis_alloc_R_1d(V, NDOF*NP)
    call monolis_alloc_R_1d(T, NDOF*NP)
    call monolis_alloc_R_2d(P, NDOF*NP, S)
    call monolis_alloc_R_2d(U, NDOF*NP, S)
    call monolis_alloc_R_2d(G, NDOF*NP, S)
    call monolis_alloc_R_2d(E, S, S)
    call monolis_alloc_R_1d(F, S)
    call monolis_alloc_R_1d(alpha, S)
    call monolis_alloc_R_1d(C, S*S + S)

    call random_seed()
    call random_number(P)

    do i = 1, S
      call monolis_inner_product_main_R_no_comm(N*NDOF, P(:,i), P(:,i), C(i))
    end do
    call monolis_allreduce_R(S, C, monolis_mpi_sum, monoCOM%comm)
    do i = 1, S
      P(:,i) = P(:,i) / C(i)
    end do

    call monolis_precond_apply_R(monoPRM, monoCOM, monoMAT, monoPREC, R, U(:,1))
    call monolis_matvec_product_main_R(monoCOM, monoMAT, U(:,1), G(:,1), tspmv, tcomm_spmv)
    do i = 2, S
      U(:,i) = U(:,1)
      G(:,i) = G(:,1)
    end do

    call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
    call monolis_set_converge_R(monoCOM, monoMAT, R, B2, is_converge, tdotp, tcomm_dotp)
    if(is_converge) return

    do iter = 1, monoPRM%Iarray(monolis_prm_I_max_iter)
      do i = 1, S
      do j = 1, S
        call monolis_inner_product_main_R_no_comm(N*NDOF, P(:,i), G(:,j), C(S*(i-1) + j))
      enddo
      enddo
      do i = 1, S
        call monolis_inner_product_main_R_no_comm(N*NDOF, P(:,i), R, C(S*S + i))
      enddo
      call monolis_allreduce_R(S*S + S, C, monolis_mpi_sum, monoCOM%comm)
      do i = 1, S
      do j = 1, S
        E(i,j) = (S*(i-1) + j)
      enddo
      enddo
      do i = 1, S
        F(i) = C(S*S + i)
      enddo

      call monolis_lapack_dsysv(S, E, F, alpha)

      call monolis_vec_copy_R(N*NDOF, R, V)

      do i = 1, S
        call monolis_vec_AXPBY_R(N*NDOF, -alpha(i), G(:,i), 1.0d0, V, V)
      enddo

      call monolis_precond_apply_R(monoPRM, monoCOM, monoMAT, monoPREC, V, Z)

      call monolis_matvec_product_main_R(monoCOM, monoMAT, Z, T, tspmv, tcomm_spmv)

      call monolis_inner_product_main_R_no_comm(N*NDOF, T, V, C(1))
      call monolis_inner_product_main_R_no_comm(N*NDOF, T, T, C(2))
      call monolis_allreduce_R(2, C, monolis_mpi_sum, monoCOM%comm)
      a1 = C(1)
      a2 = C(2)
      beta = a1/a2

      do i = 1, S
        call monolis_vec_AXPBY_R(N*NDOF, alpha(i), U(:,i), 1.0d0, X, X)
      enddo

      call monolis_vec_AXPBY_R(N*NDOF, beta, Z, 1.0d0, X, X)

      if(mod(iter, iter_RR) == 0)then
        call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
      else
        call monolis_vec_AXPBY_R(N*NDOF, -beta, T, 1.0d0, V, R)
      endif

      call monolis_check_converge_R(monoPRM, monoCOM, monoMAT, R, B2, iter, is_converge, tdotp, tcomm_dotp)
      if(is_converge) exit

      do i = 1, S - 1
        G(:,i) = G(:,i + 1)
        U(:,i) = U(:,i + 1)
      enddo

      call monolis_precond_apply_R(monoPRM, monoCOM, monoMAT, monoPREC, R, U(:,S))
      call monolis_matvec_product_main_R(monoCOM, monoMAT, U(:,S), G(:,S), tspmv, tcomm_spmv)
    enddo

    call monolis_mpi_update_R(monoCOM, NDOF, X, tcomm_spmv)

    monoPRM%Rarray(monolis_R_time_spmv) = tspmv
    monoPRM%Rarray(monolis_R_time_comm_spmv) = tcomm_spmv
    monoPRM%Rarray(monolis_R_time_dotp) = tdotp
    monoPRM%Rarray(monolis_R_time_comm_dotp) = tcomm_dotp

    call monolis_dealloc_R_1d(R)
    call monolis_dealloc_R_1d(Z)
    call monolis_dealloc_R_1d(V)
    call monolis_dealloc_R_1d(T)
    call monolis_dealloc_R_2d(P)
    call monolis_dealloc_R_2d(U)
    call monolis_dealloc_R_2d(G)
    call monolis_dealloc_R_2d(E)
    call monolis_dealloc_R_1d(F)
    call monolis_dealloc_R_1d(alpha)
    call monolis_dealloc_R_1d(C)
  end subroutine monolis_solver_IDRS

end module mod_monolis_solver_IDRS
