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
    integer(kint) :: NNDOF, NPNDOF, NGDOF
    integer(kint) :: i, j, k, ipass, iter, iter_RR, S
    real(kdouble) :: B2, R2, omega, alpha, beta, Q(3), rho, kappa, tol
    real(kdouble) :: tspmv, tdotp, tcomm_spmv, tcomm_dotp
    logical :: is_converge, is_candidate, is_restart
    real(kdouble), allocatable :: R(:), P(:,:), U(:,:), G(:,:), M(:,:)
    real(kdouble), allocatable :: F(:), Z(:), V(:), T(:), C(:)
    real(kdouble), pointer, contiguous :: B(:), X(:)

    call monolis_std_debug_log_header("monolis_solver_IDRS")

    X => monoMAT%R%X
    B => monoMAT%R%B
    iter_RR = 50
    if(monoPRM%Iarray(monolis_prm_I_iter_RR) > 0) iter_RR = monoPRM%Iarray(monolis_prm_I_iter_RR)
    S = monoPRM%Iarray(MONOLIS_PRM_I_IDRS_DIM)
    tol = monoPRM%Rarray(monolis_prm_R_tol)

    tspmv = monoPRM%Rarray(monolis_R_time_spmv)
    tcomm_spmv = monoPRM%Rarray(monolis_R_time_comm_spmv)
    tdotp = monoPRM%Rarray(monolis_R_time_dotp)
    tcomm_dotp = monoPRM%Rarray(monolis_R_time_comm_dotp)

    if(monoPRM%Iarray(monolis_prm_I_is_init_x) == monolis_I_true)then
      X = 0.0d0
    endif

    call monolis_get_vec_size(monoMAT%N, monoMAT%NP, monoMAT%NDOF, &
      monoMAT%n_dof_index, NNDOF, NPNDOF)

    !# 直交基底数は袖を除いた全体自由度数以下とする
    NGDOF = NNDOF
    call monolis_allreduce_I1(NGDOF, monolis_mpi_sum, monoCOM%comm)
    S = min(S, NGDOF)
    if(S < 1)then
      call monolis_std_error_string("IDR(s) requires a positive number of basis vectors and unknowns")
      call monolis_std_error_stop()
    endif

    call monolis_alloc_R_1d(R, NPNDOF)
    call monolis_alloc_R_1d(Z, NPNDOF)
    call monolis_alloc_R_1d(V, NPNDOF)
    call monolis_alloc_R_1d(T, NPNDOF)
    call monolis_alloc_R_2d(P, NPNDOF, S)
    call monolis_alloc_R_2d(U, NPNDOF, S)
    call monolis_alloc_R_2d(G, NPNDOF, S)
    call monolis_alloc_R_2d(M, S, S)
    call monolis_alloc_R_1d(C, S)
    call monolis_alloc_R_1d(F, S)

    call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
    call monolis_set_converge_R(monoCOM, monoMAT, B, B2, is_converge, tdotp, tcomm_dotp)
    if(is_converge)then
      monoPRM%Iarray(monolis_prm_I_cur_iter) = 0
      monoPRM%Rarray(monolis_prm_R_cur_resid) = 0.0d0
    else
      !# 初期解が既に収束していれば基底生成や除算を行わない
      call monolis_check_converge_R(monoPRM, monoCOM, monoMAT, R, B2, 0, is_converge, tdotp, tcomm_dotp)
    endif

    if(.not. is_converge)then

    !# P は残差と独立な乱数基底とする。残差への直交化は P^T R を退化させる
    call random_number(P)
    P = P - 0.5d0
    do i = 1, S
      !# 修正 Gram-Schmidt を二回適用し、各列を直交化した後に正規化する
      do ipass = 1, 2
        do j = 1, i - 1
          call monolis_inner_product_main_R(monoCOM, NNDOF, P(:,j), P(:,i), alpha, tdotp, tcomm_dotp)
          call monolis_vec_AXPBY_R(NNDOF, -alpha, P(:,j), 1.0d0, P(:,i), P(:,i))
        enddo
      enddo
      call monolis_inner_product_main_R(monoCOM, NNDOF, P(:,i), P(:,i), alpha, tdotp, tcomm_dotp)
      if(alpha == 0.0d0)then
        call monolis_std_error_string("IDR(s) shadow space is rank deficient")
        call monolis_std_error_stop()
      endif
      P(:,i) = P(:,i)/dsqrt(alpha)
    enddo

    do i = 1, S
      M(i,i) = 1.0d0
    enddo

    omega = 1.0d0
    kappa = 0.7d0

    do iter = 1, monoPRM%Iarray(monolis_prm_I_max_iter)
      !# f の更新
      do i = 1, S
        call monolis_inner_product_main_R_no_comm(NNDOF, P(:,i), R, F(i))
      enddo
      call monolis_allreduce_R(S, F, monolis_mpi_sum, monoCOM%comm)

      do k = 1, S
        !# 既に消去した成分を除き M(k:S,k:S) C(k:S) = F(k:S) を解く
        C = 0.0d0
        C(k:S) = F(k:S)
        do i = k, S
          do j = k, i - 1
            C(i) = C(i) - M(i,j) * C(j)
          enddo
          if(M(i,i) == 0.0d0) stop "zero divide A"
          C(i) = C(i) / M(i,i)
        enddo

        !# V の更新
        V = R
        do i = k, S
          V = V - C(i)*G(:,i)
        enddo

        !# Z = M V
        call monolis_precond_apply_R(monoPRM, monoCOM, monoMAT, monoPREC, V, Z)

        !# U の更新
        U(:,k) = omega*Z + C(k)*U(:,k)
        do i = k + 1, S
          U(:,k) = U(:,k) + C(i)*U(:,i)
        enddo

        !# G_k = A U_k
        call monolis_matvec_product_main_R(monoCOM, monoMAT, U(:,k), G(:,k), tspmv, tcomm_spmv)

        do i = 1, k - 1
          call monolis_inner_product_main_R(monoCOM, NNDOF, P(:,i), G(:,k), alpha, tdotp, tcomm_dotp)
          if(M(i,i) == 0.0d0) stop "zero divide B"
          alpha = alpha / M(i,i)

          G(:,k) = G(:,k) - alpha*G(:,i)
          U(:,k) = U(:,k) - alpha*U(:,i)
        enddo

        do i = k, S
          call monolis_inner_product_main_R(monoCOM, NNDOF, P(:,i), G(:,k), M(i,k), tdotp, tcomm_dotp)
        enddo

        if(M(k,k) == 0.0d0) stop "IDR(s) breakdown: zero shadow-space pivot"
        beta = F(k) / M(k,k)

        call monolis_vec_AXPBY_R(NNDOF,  beta, U(:,k), 1.0d0, X, X)
        call monolis_vec_AXPBY_R(NNDOF, -beta, G(:,k), 1.0d0, R, R)

        call monolis_inner_product_main_R(monoCOM, NNDOF, R, R, R2, tdotp, tcomm_dotp)
        is_candidate = dsqrt(R2/B2) < tol
        if(is_candidate) exit

        if(k < S)then
          do i = k + 1, S
            F(i) = F(i) - beta*M(i,k)
          enddo
        endif
      enddo

      if(.not. is_candidate)then
        call monolis_precond_apply_R(monoPRM, monoCOM, monoMAT, monoPREC, R, Z)
        call monolis_matvec_product_main_R(monoCOM, monoMAT, Z, T, tspmv, tcomm_spmv)

        call monolis_inner_product_main_R_no_comm(NNDOF, T, R, Q(1))
        call monolis_inner_product_main_R_no_comm(NNDOF, T, T, Q(2))
        call monolis_inner_product_main_R_no_comm(NNDOF, R, R, Q(3))
        call monolis_allreduce_R(3, Q, monolis_mpi_sum, monoCOM%comm)
        if(Q(2) == 0.0d0) stop "zero divide E"
        rho = Q(1)/( dsqrt(Q(2))*dsqrt(Q(3)) )
        if(rho == 0.0d0) stop "zero divide D"

        omega = Q(1)/Q(2)
        if(dabs(rho) < kappa)then
          omega = omega*kappa/dabs(rho)
        endif

        call monolis_vec_AXPBY_R(NNDOF,  omega, Z, 1.0d0, X, X)
        call monolis_vec_AXPBY_R(NNDOF, -omega, T, 1.0d0, R, R)

        call monolis_inner_product_main_R(monoCOM, NNDOF, R, R, R2, tdotp, tcomm_dotp)
        is_candidate = dsqrt(R2/B2) < tol
      endif

      !# 収束候補と最終反復では必ず真の残差で判定する。最大反復判定は omega 更新後に行う
      is_restart = is_candidate .or. mod(iter, iter_RR) == 0 .or. iter == monoPRM%Iarray(monolis_prm_I_max_iter)
      if(is_restart)then
        call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
      endif
      call monolis_check_converge_R(monoPRM, monoCOM, monoMAT, R, B2, iter, is_converge, tdotp, tcomm_dotp)
      if(is_converge) exit

      if(is_restart)then
        !# 残差置換後は古い F と直交性を持ち越さず、新しい残差から再開する
        U = 0.0d0
        G = 0.0d0
        M = 0.0d0
        do i = 1, S
          M(i,i) = 1.0d0
        enddo
        omega = 1.0d0
      endif
    enddo

    endif

    call monolis_mpi_update_R_wrapper(monoCOM, monoMAT%NDOF, monoMAT%n_dof_index, X, tcomm_spmv)

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
    call monolis_dealloc_R_2d(M)
    call monolis_dealloc_R_1d(C)
    call monolis_dealloc_R_1d(F)
  end subroutine monolis_solver_IDRS

end module mod_monolis_solver_IDRS
