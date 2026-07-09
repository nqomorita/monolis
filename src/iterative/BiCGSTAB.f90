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
    integer(kint) :: NNDOF, NPNDOF
    integer(kint) :: i, iter, iter_RR
    real(kdouble) :: alpha, beta, rho, rho1, c2, omega
    real(kdouble) :: B2, CG(2)
    real(kdouble) :: tspmv, tdotp, tcomm_spmv, tcomm_dotp
    real(kdouble), allocatable :: R(:), RT(:), P(:), PT(:), S(:), ST(:), T(:), V(:)
    real(kdouble), pointer, contiguous :: B(:), X(:)
    logical :: is_converge

    call monolis_std_debug_log_header("monolis_solver_BiCGSTAB")

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

    call monolis_get_vec_size(monoMAT%N, monoMAT%NP, monoMAT%NDOF, &
      monoMAT%n_dof_index, NNDOF, NPNDOF)

    call monolis_alloc_R_1d(R , NPNDOF)
    call monolis_alloc_R_1d(RT, NPNDOF)
    call monolis_alloc_R_1d(P , NPNDOF)
    call monolis_alloc_R_1d(PT, NPNDOF)
    call monolis_alloc_R_1d(S , NPNDOF)
    call monolis_alloc_R_1d(ST, NPNDOF)
    call monolis_alloc_R_1d(T , NPNDOF)
    call monolis_alloc_R_1d(V , NPNDOF)

    !# OpenACC: ソルバ固有のワーク配列のみデバイスに確保（X/B/precD は外側で常駐）
    !$acc enter data create(R(1:NPNDOF), RT(1:NPNDOF), P(1:NPNDOF), PT(1:NPNDOF), &
    !$acc                   S(1:NPNDOF), ST(1:NPNDOF), T(1:NPNDOF), V(1:NPNDOF))

    call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
    call monolis_set_converge_R(monoCOM, monoMAT, R, B2, is_converge, tdotp, tcomm_dotp)
    if(is_converge)then
      !$acc update self(X(1:NPNDOF))
      !$acc exit data delete(R, RT, P, PT, S, ST, T, V)
      call monolis_dealloc_R_1d(R )
      call monolis_dealloc_R_1d(RT)
      call monolis_dealloc_R_1d(P )
      call monolis_dealloc_R_1d(PT)
      call monolis_dealloc_R_1d(S )
      call monolis_dealloc_R_1d(ST)
      call monolis_dealloc_R_1d(T )
      call monolis_dealloc_R_1d(V )
      return
    endif

    call monolis_vec_copy_R(NNDOF, R, RT)

    do iter = 1, monoPRM%Iarray(monolis_prm_I_max_iter)
      call monolis_inner_product_main_R(monoCOM, NNDOF, R, RT, rho, tdotp, tcomm_dotp)
      if(dabs(rho) == 0.0d0)then
        !$acc update self(R(1:NNDOF), RT(1:NNDOF))
        call monolis_global_sorted_inner_product_main_R(monoCOM, &
          NNDOF, R, RT, rho, tdotp, tcomm_dotp)
      endif

      if(1 < iter)then
        beta = (rho/rho1) * (alpha/omega)
!$omp parallel default(none) &
!$omp & shared(P, R, V) &
!$omp & firstprivate(NNDOF, beta, omega) &
!$omp & private(i)
!$omp do
!$acc parallel loop present(P, R, V)
        do i = 1, NNDOF
          P(i) = R(i) + beta * (P(i) - omega * V(i))
        enddo
!$acc end parallel loop
!$omp end do
!$omp end parallel
      else
        call monolis_vec_copy_R(NNDOF, R, P)
      endif

      call monolis_precond_apply_R(monoPRM, monoCOM, monoMAT, monoPREC, P, PT)
      call monolis_matvec_product_main_R(monoCOM, monoMAT, PT, V, tspmv, tcomm_spmv)
      call monolis_inner_product_main_R(monoCOM, NNDOF, RT, V, c2, tdotp, tcomm_dotp)
      if(dabs(c2) == 0.0d0)then
        !$acc update self(RT(1:NNDOF), V(1:NNDOF))
        call monolis_global_sorted_inner_product_main_R(monoCOM, &
          NNDOF, RT, V, c2, tdotp, tcomm_dotp)
      endif

      alpha = rho / c2
      call monolis_vec_AXPBY_R(NNDOF, -alpha, V, 1.0d0, R, S)

      call monolis_precond_apply_R(monoPRM, monoCOM, monoMAT, monoPREC, S, ST)
      call monolis_matvec_product_main_R(monoCOM, monoMAT, ST, T, tspmv, tcomm_spmv)

      call monolis_inner_product_main_R_no_comm(NNDOF, T, S, CG(1))
      call monolis_inner_product_main_R_no_comm(NNDOF, T, T, CG(2))
      !$acc data copy(CG)
      call monolis_allreduce_R(2, CG, monolis_mpi_sum, monoCOM%comm)
      !$acc end data
      if(dabs(CG(2)) == 0.0d0)then
        !$acc update self(T(1:NNDOF))
        call monolis_global_sorted_inner_product_main_R(monoCOM, &
          NNDOF, T, T, CG(2), tdotp, tcomm_dotp)
      endif

      if(dabs(CG(2)) == 0.0d0)then
        !# T = 0 の場合（完全 LU 前処理で S = 0 となる場合など）は
        !# 0/0 を避けて omega = 0 とする（R = S = 0 のため次の収束判定で終了する）
        omega = 0.0d0
      else
        omega = CG(1) / CG(2)
      endif

!$omp parallel default(none) &
!$omp & shared(X, PT, ST) &
!$omp & firstprivate(NNDOF, alpha, omega) &
!$omp & private(i)
!$omp do
!$acc parallel loop present(X, PT, ST)
      do i = 1, NNDOF
        X(i) = X(i) + alpha*PT(i) + omega*ST(i)
      enddo
!$acc end parallel loop
!$omp end do
!$omp end parallel

      if(mod(iter, iter_RR) == 0)then
        call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
      else
        call monolis_vec_AXPBY_R(NNDOF, -omega, T, 1.0d0, S, R)
      endif

      call monolis_check_converge_R(monoPRM, monoCOM, monoMAT, R, B2, iter, is_converge, tdotp, tcomm_dotp)
      if(is_converge) exit

      rho1 = rho
    enddo

    !$acc update self(X(1:NPNDOF))

    call monolis_mpi_update_R_wrapper(monoCOM, monoMAT%NDOF, monoMAT%n_dof_index, X, tcomm_spmv)

    monoPRM%Rarray(monolis_R_time_spmv) = tspmv
    monoPRM%Rarray(monolis_R_time_comm_spmv) = tcomm_spmv
    monoPRM%Rarray(monolis_R_time_dotp) = tdotp
    monoPRM%Rarray(monolis_R_time_comm_dotp) = tcomm_dotp

    !# OpenACC: ワーク配列のみ破棄（X/B/precD の破棄は外側で実施）
    !$acc exit data delete(R, RT, P, PT, S, ST, T, V)

    call monolis_dealloc_R_1d(R )
    call monolis_dealloc_R_1d(RT)
    call monolis_dealloc_R_1d(P )
    call monolis_dealloc_R_1d(PT)
    call monolis_dealloc_R_1d(S )
    call monolis_dealloc_R_1d(ST)
    call monolis_dealloc_R_1d(T )
    call monolis_dealloc_R_1d(V )
  end subroutine monolis_solver_BiCGSTAB

  !> @ingroup solver
  !> BiCGStab 法（GPU 向けカーネル融合版）
  !> @details
  !> 数理アルゴリズムは monolis_solver_BiCGSTAB と完全に等価。
  !> OpenACC カーネル起動回数を削減するため、以下の演算を融合する。
  !>   - Fusion A : 内積 (T,S) と (T,T) を 1 カーネルに統合（T を 1 回だけ読む）
  !>   - Fusion B : 解ベクトル更新 X、残差更新 R、残差ノルム (R,R) を 1 カーネルに統合
  !>                し、収束判定をインライン化する。
  subroutine monolis_solver_BiCGSTAB_GPU(monoPRM, monoCOM, monoMAT, monoPREC)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in,out] 行列構造体
    type(monolis_mat), target, intent(inout) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC
    integer(kint) :: NNDOF, NPNDOF
    integer(kint) :: i, iter, iter_RR
    real(kdouble) :: alpha, beta, rho, rho1, c2, omega
    real(kdouble) :: B2, CG(2)
    real(kdouble) :: dot_ts, dot_tt, R2, resid, rtmp
    real(kdouble) :: t1, t2, t3
    real(kdouble) :: tspmv, tdotp, tcomm_spmv, tcomm_dotp
    real(kdouble), allocatable :: R(:), RT(:), P(:), PT(:), S(:), ST(:), T(:), V(:)
    real(kdouble), pointer, contiguous :: B(:), X(:)
    logical :: is_converge

    call monolis_std_debug_log_header("monolis_solver_BiCGSTAB_GPU")

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

    call monolis_get_vec_size(monoMAT%N, monoMAT%NP, monoMAT%NDOF, &
      monoMAT%n_dof_index, NNDOF, NPNDOF)

    call monolis_alloc_R_1d(R , NPNDOF)
    call monolis_alloc_R_1d(RT, NPNDOF)
    call monolis_alloc_R_1d(P , NPNDOF)
    call monolis_alloc_R_1d(PT, NPNDOF)
    call monolis_alloc_R_1d(S , NPNDOF)
    call monolis_alloc_R_1d(ST, NPNDOF)
    call monolis_alloc_R_1d(T , NPNDOF)
    call monolis_alloc_R_1d(V , NPNDOF)

    !# OpenACC: ソルバ固有のワーク配列のみデバイスに確保（X/B/precD は外側で常駐）
    !$acc enter data create(R(1:NPNDOF), RT(1:NPNDOF), P(1:NPNDOF), PT(1:NPNDOF), &
    !$acc                   S(1:NPNDOF), ST(1:NPNDOF), T(1:NPNDOF), V(1:NPNDOF))

    call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
    call monolis_set_converge_R(monoCOM, monoMAT, R, B2, is_converge, tdotp, tcomm_dotp)
    if(is_converge)then
      !$acc update self(X(1:NPNDOF))
      !$acc exit data delete(R, RT, P, PT, S, ST, T, V)
      call monolis_dealloc_R_1d(R )
      call monolis_dealloc_R_1d(RT)
      call monolis_dealloc_R_1d(P )
      call monolis_dealloc_R_1d(PT)
      call monolis_dealloc_R_1d(S )
      call monolis_dealloc_R_1d(ST)
      call monolis_dealloc_R_1d(T )
      call monolis_dealloc_R_1d(V )
      return
    endif

    call monolis_vec_copy_R(NNDOF, R, RT)

    do iter = 1, monoPRM%Iarray(monolis_prm_I_max_iter)
      call monolis_inner_product_main_R(monoCOM, NNDOF, R, RT, rho, tdotp, tcomm_dotp)
      if(dabs(rho) == 0.0d0)then
        !$acc update self(R(1:NNDOF), RT(1:NNDOF))
        call monolis_global_sorted_inner_product_main_R(monoCOM, &
          NNDOF, R, RT, rho, tdotp, tcomm_dotp)
      endif

      if(1 < iter)then
        beta = (rho/rho1) * (alpha/omega)
!$omp parallel default(none) &
!$omp & shared(P, R, V) &
!$omp & firstprivate(NNDOF, beta, omega) &
!$omp & private(i)
!$omp do
!$acc parallel loop present(P, R, V)
        do i = 1, NNDOF
          P(i) = R(i) + beta * (P(i) - omega * V(i))
        enddo
!$acc end parallel loop
!$omp end do
!$omp end parallel
      else
        call monolis_vec_copy_R(NNDOF, R, P)
      endif

      call monolis_precond_apply_R(monoPRM, monoCOM, monoMAT, monoPREC, P, PT)
      call monolis_matvec_product_main_R(monoCOM, monoMAT, PT, V, tspmv, tcomm_spmv)
      call monolis_inner_product_main_R(monoCOM, NNDOF, RT, V, c2, tdotp, tcomm_dotp)
      if(dabs(c2) == 0.0d0)then
        !$acc update self(RT(1:NNDOF), V(1:NNDOF))
        call monolis_global_sorted_inner_product_main_R(monoCOM, &
          NNDOF, RT, V, c2, tdotp, tcomm_dotp)
      endif

      alpha = rho / c2
      call monolis_vec_AXPBY_R(NNDOF, -alpha, V, 1.0d0, R, S)

      call monolis_precond_apply_R(monoPRM, monoCOM, monoMAT, monoPREC, S, ST)
      call monolis_matvec_product_main_R(monoCOM, monoMAT, ST, T, tspmv, tcomm_spmv)

      !# Fusion A: 内積 (T,S) と (T,T) を 1 カーネルに融合（T を 1 回だけロード）
      t1 = monolis_get_time()
      dot_ts = 0.0d0
      dot_tt = 0.0d0
!$omp parallel default(none) &
!$omp & shared(T, S, dot_ts, dot_tt) &
!$omp & firstprivate(NNDOF) &
!$omp & private(i)
!$omp do reduction(+:dot_ts, dot_tt)
!$acc parallel loop present(T, S) reduction(+:dot_ts, dot_tt)
      do i = 1, NNDOF
        dot_ts = dot_ts + T(i)*S(i)
        dot_tt = dot_tt + T(i)*T(i)
      enddo
!$acc end parallel loop
!$omp end do
!$omp end parallel
      CG(1) = dot_ts
      CG(2) = dot_tt
      t2 = monolis_get_time()
      !$acc data copy(CG)
      call monolis_allreduce_R(2, CG, monolis_mpi_sum, monoCOM%comm)
      !$acc end data
      t3 = monolis_get_time()
      tdotp = tdotp + t3 - t1
      tcomm_dotp = tcomm_dotp + t3 - t2
      if(dabs(CG(2)) == 0.0d0)then
        !$acc update self(T(1:NNDOF))
        call monolis_global_sorted_inner_product_main_R(monoCOM, &
          NNDOF, T, T, CG(2), tdotp, tcomm_dotp)
      endif

      if(dabs(CG(2)) == 0.0d0)then
        !# T = 0 の場合（完全 LU 前処理で S = 0 となる場合など）は
        !# 0/0 を避けて omega = 0 とする（R = S = 0 のため次の収束判定で終了する）
        omega = 0.0d0
      else
        omega = CG(1) / CG(2)
      endif

      if(mod(iter, iter_RR) == 0)then
        !# 周期的な真の残差再計算（数値安定化）はオリジナルと同一経路
!$omp parallel default(none) &
!$omp & shared(X, PT, ST) &
!$omp & firstprivate(NNDOF, alpha, omega) &
!$omp & private(i)
!$omp do
!$acc parallel loop present(X, PT, ST)
        do i = 1, NNDOF
          X(i) = X(i) + alpha*PT(i) + omega*ST(i)
        enddo
!$acc end parallel loop
!$omp end do
!$omp end parallel

        call monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm_spmv)
        call monolis_check_converge_R(monoPRM, monoCOM, monoMAT, R, B2, iter, is_converge, tdotp, tcomm_dotp)
      else
        !# Fusion B: X 更新・R 更新・残差ノルム (R,R) を 1 カーネルに融合
        t1 = monolis_get_time()
        R2 = 0.0d0
!$omp parallel default(none) &
!$omp & shared(X, PT, ST, R, S, T, R2) &
!$omp & firstprivate(NNDOF, alpha, omega) &
!$omp & private(i, rtmp)
!$omp do reduction(+:R2)
!$acc parallel loop present(X, PT, ST, R, S, T) reduction(+:R2) private(rtmp)
        do i = 1, NNDOF
          X(i) = X(i) + alpha*PT(i) + omega*ST(i)
          rtmp = S(i) - omega*T(i)
          R(i) = rtmp
          R2 = R2 + rtmp*rtmp
        enddo
!$acc end parallel loop
!$omp end do
!$omp end parallel
        t2 = monolis_get_time()
        call monolis_allreduce_R1(R2, monolis_mpi_sum, monoCOM%comm)
        t3 = monolis_get_time()
        tdotp = tdotp + t3 - t1
        tcomm_dotp = tcomm_dotp + t3 - t2

        !# 収束判定のインライン化（monolis_check_converge_R と等価）
        resid = dsqrt(R2/B2)
        monoPRM%Iarray(monolis_prm_I_cur_iter) = iter
        monoPRM%Rarray(monolis_prm_R_cur_resid) = resid
        if(monoCOM%my_rank == 0 .and. monoPRM%Iarray(monolis_prm_I_show_iterlog) == monolis_I_true)then
          write (*,"(i7, 1pe16.6)") iter, resid
        endif
        if(resid < monoPRM%Rarray(monolis_prm_R_tol))then
          is_converge = .true.
        else
          is_converge = .false.
          if(iter == monoPRM%Iarray(monolis_prm_I_max_iter))then
            if(monoPRM%Iarray(monolis_prm_I_is_measurement) == monolis_I_false .and. &
               monoPRM%Iarray(monolis_prm_I_is_error_abort) == monolis_I_true)then
              call monolis_std_error_string("reached the maximum number of iterations")
              call monolis_std_error_stop()
            endif
          endif
        endif
      endif

      if(is_converge) exit

      rho1 = rho
    enddo

    !$acc update self(X(1:NPNDOF))

    call monolis_mpi_update_R_wrapper(monoCOM, monoMAT%NDOF, monoMAT%n_dof_index, X, tcomm_spmv)

    !# OpenACC: monolis_mpi_update_R_wrapper は device 常駐の X 上で ghost 成分を更新するため、
    !#          呼び出し元（Newton 更新等）が正しい解を読めるよう host に戻す
    !$acc update self(X(1:NPNDOF))

    monoPRM%Rarray(monolis_R_time_spmv) = tspmv
    monoPRM%Rarray(monolis_R_time_comm_spmv) = tcomm_spmv
    monoPRM%Rarray(monolis_R_time_dotp) = tdotp
    monoPRM%Rarray(monolis_R_time_comm_dotp) = tcomm_dotp

    !# OpenACC: ワーク配列のみ破棄（X/B/precD の破棄は外側で実施）
    !$acc exit data delete(R, RT, P, PT, S, ST, T, V)

    call monolis_dealloc_R_1d(R )
    call monolis_dealloc_R_1d(RT)
    call monolis_dealloc_R_1d(P )
    call monolis_dealloc_R_1d(PT)
    call monolis_dealloc_R_1d(S )
    call monolis_dealloc_R_1d(ST)
    call monolis_dealloc_R_1d(T )
    call monolis_dealloc_R_1d(V )
  end subroutine monolis_solver_BiCGSTAB_GPU

end module mod_monolis_solver_BiCGSTAB
