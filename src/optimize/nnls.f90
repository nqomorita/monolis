!> NNLS 関数
module mod_monolis_opt_nnls
  use mod_monolis_utils
  use mod_monolis_lapack
  use mod_monolis_scalapack
  implicit none

contains

  !> @ingroup linalg
  !> Non-Negative Least Squares（Lawson-Hanson 法）
  !> @details 行列 A は列ブロック分散とする（各ランクが m x n_loc の列ブロックを保持し、
  !> 全体の列順序はランク番号順に連結した順序とする）。右辺ベクトル b は全ランクに複製する。
  !> KKT 条件（非 active 列の双対変数が tol 以下）を満たすまで反復し、NNLS の最適解を求める。
  !> 逐次実行する場合はセルフコミュニケータを渡す。
  subroutine monolis_optimize_nnls_R(A, b, x, m, n_loc, max_iter, tol, residual, comm)
    implicit none
    !> [in] 行列（大きさ m x n_loc、列ブロック分散）
    real(kdouble), intent(in) :: A(:,:)
    !> [in] 右辺ベクトル（大きさ m、全ランク複製）
    real(kdouble), intent(in) :: b(:)
    !> [out] 解ベクトル（大きさ n_loc、自ランク保持列に対応）
    real(kdouble), intent(out) :: x(:)
    !> [in] 行列の大きさ（行数）
    integer(kint), intent(in) :: m
    !> [in] 行列の大きさ（自ランク保持列数）
    integer(kint), intent(in) :: n_loc
    !> [in] 最大反復回数
    integer(kint), intent(in) :: max_iter
    !> [in] 収束判定閾値（KKT 条件の判定と微小解成分の除外に使用）
    real(kdouble), intent(in) :: tol
    !> [out] 残差
    real(kdouble), intent(out) :: residual
    !> [in] MPI コミュニケータ
    integer(kint), intent(in) :: comm

    call monolis_std_debug_log_header("monolis_optimize_nnls_R")

    call monolis_optimize_nnls_LH_main(A, b, x, m, n_loc, max_iter, .false., 0.0d0, tol, residual, comm)
  end subroutine monolis_optimize_nnls_R

  !> @ingroup linalg
  !> Non-Negative Least Squares（疎な解ベクトル）
  !> @details 行列 A は列ブロック分散とする（各ランクが m x n_loc の列ブロックを保持し、
  !> 全体の列順序はランク番号順に連結した順序とする）。右辺ベクトル b は全ランクに複製する。
  !> 相対残差が tol_outer を下回った時点で反復を打ち切り、疎な解ベクトルを得る。
  !> 逐次実行する場合はセルフコミュニケータを渡す。
  subroutine monolis_optimize_nnls_R_with_sparse_solution(A, b, x, m, n_loc, max_iter, tol_outer, tol_inner, residual, comm)
    implicit none
    !> [in] 行列（大きさ m x n_loc、列ブロック分散）
    real(kdouble), intent(in) :: A(:,:)
    !> [in] 右辺ベクトル（大きさ m、全ランク複製）
    real(kdouble), intent(in) :: b(:)
    !> [out] 解ベクトル（大きさ n_loc、自ランク保持列に対応）
    real(kdouble), intent(out) :: x(:)
    !> [in] 行列の大きさ（行数）
    integer(kint), intent(in) :: m
    !> [in] 行列の大きさ（自ランク保持列数）
    integer(kint), intent(in) :: n_loc
    !> [in] 最大反復回数
    integer(kint), intent(in) :: max_iter
    !> [in] 外側反復の収束判定閾値（相対残差 r_norm/r0_norm に対する閾値）
    real(kdouble), intent(in) :: tol_outer
    !> [in] 内部 NNLS（Lawson-Hanson 反復）の閾値（この値以下の解成分を active 集合から除外）
    real(kdouble), intent(in) :: tol_inner
    !> [out] 残差
    real(kdouble), intent(out) :: residual
    !> [in] MPI コミュニケータ
    integer(kint), intent(in) :: comm

    call monolis_std_debug_log_header("monolis_optimize_nnls_R_with_sparse_solution")

    call monolis_optimize_nnls_LH_main(A, b, x, m, n_loc, max_iter, .true., tol_outer, tol_inner, residual, comm)
  end subroutine monolis_optimize_nnls_R_with_sparse_solution

  !> @ingroup dev_linalg
  !> Non-Negative Least Squares（Lawson-Hanson 法、メイン関数）
  !> @details active 列の候補探索は列分散で並列化し、active 集合に対する正規方程式は
  !> ScaLAPACK（PDGESV）で分散求解する。is_sparse が真の場合は相対残差 tol_outer による
  !> 早期終了を行い、偽の場合は KKT 条件（双対変数 <= tol_inner）まで反復する。
  subroutine monolis_optimize_nnls_LH_main(A, b, x, m, n_loc, max_iter, is_sparse, tol_outer, tol_inner, residual, comm)
    implicit none
    !> [in] 行列（大きさ m x n_loc、列ブロック分散）
    real(kdouble), intent(in) :: A(:,:)
    !> [in] 右辺ベクトル（大きさ m、全ランク複製）
    real(kdouble), intent(in) :: b(:)
    !> [out] 解ベクトル（大きさ n_loc、自ランク保持列に対応）
    real(kdouble), intent(out) :: x(:)
    !> [in] 行列の大きさ（行数）
    integer(kint), intent(in) :: m
    !> [in] 行列の大きさ（自ランク保持列数）
    integer(kint), intent(in) :: n_loc
    !> [in] 最大反復回数
    integer(kint), intent(in) :: max_iter
    !> [in] 相対残差による早期終了の有無（疎な解ベクトルの取得）
    logical, intent(in) :: is_sparse
    !> [in] 外側反復の収束判定閾値（is_sparse が真の場合に使用）
    real(kdouble), intent(in) :: tol_outer
    !> [in] KKT 条件の判定と微小解成分の除外に使用する閾値
    real(kdouble), intent(in) :: tol_inner
    !> [out] 残差
    real(kdouble), intent(out) :: residual
    !> [in] MPI コミュニケータ
    integer(kint), intent(in) :: comm
    integer(kint) :: scalapack_comm, my_rank, n_rank, n_global
    integer(kint) :: iter, inner_iter, i, k, in, idx, owner, p_loc, p_global, neg
    real(kdouble) :: r0_norm, r_norm, s_max, g_max, alpha, t, drop_tol
    integer(kint), allocatable :: map(:)
    real(kdouble), allocatable :: r(:)
    real(kdouble), allocatable :: s(:)
    real(kdouble), allocatable :: Ax(:)
    real(kdouble), allocatable :: A_z(:,:)
    real(kdouble), allocatable :: w_z(:)
    logical, allocatable :: is_nonzero(:)
    logical :: is_converge

    my_rank = monolis_mpi_get_local_my_rank(comm)
    n_rank = monolis_mpi_get_local_comm_size(comm)

    n_global = n_loc
    call monolis_allreduce_I1(n_global, monolis_mpi_sum, comm)

    call monolis_scalapack_comm_initialize(comm, scalapack_comm)

    !> メモリの確保
    call monolis_alloc_R_1d(r, m)
    call monolis_alloc_R_1d(s, n_loc)
    call monolis_alloc_R_1d(Ax, m)
    call monolis_alloc_L_1d(is_nonzero, n_loc)

    !> 収束判定のためのノルム計算（b は複製のため全ランクで同一値）
    call monolis_get_l2_norm_R(m, b, r0_norm)

    is_nonzero = .false.
    is_converge = .false.
    r = b
    x = 0.0d0
    p_global = 0
    drop_tol = tol_inner

    do iter = 1, max_iter
      !> 行列 A の転置と残差ベクトルをかける（自ランク保持列のみ）
      if(n_loc > 0)then
        call dgemv("T", m, n_loc, 1.0d0, A, m, r, 1, 0.0d0, s, 1)
      endif

      !> 非 active 列に対する最大値の全体探索
      idx = 0
      s_max = -huge(0.0d0)
      do i = 1, n_loc
        if(is_nonzero(i)) cycle
        if(s(i) > s_max)then
          s_max = s(i)
          idx = i
        endif
      enddo

      g_max = s_max
      call monolis_allreduce_R1(g_max, monolis_mpi_max, comm)

      !> KKT 条件の判定（完全 NNLS モード）: 非 active 列の双対変数が全て閾値以下なら最適
      if(.not. is_sparse)then
        if(g_max <= tol_inner)then
          is_converge = .true.
          exit
        endif
      endif

      !> 最大値を持つ最小ランクを選び、その列を active に指定する
      owner = n_rank
      if(idx > 0 .and. s_max == g_max) owner = my_rank
      call monolis_allreduce_I1(owner, monolis_mpi_min, comm)
      if(owner == my_rank) is_nonzero(idx) = .true.

      !> active 集合に対する最小二乗問題の求解（Lawson-Hanson 内部反復）
      inner_iter = 0
      do
        p_loc = 0
        do i = 1, n_loc
          if(is_nonzero(i)) p_loc = p_loc + 1
        enddo
        p_global = p_loc
        call monolis_allreduce_I1(p_global, monolis_mpi_sum, comm)
        if(p_global == 0) exit

        !> 行列 A から active 列を取得
        call monolis_alloc_I_1d(map, p_loc)
        call monolis_alloc_R_2d(A_z, m, p_loc)
        call monolis_alloc_R_1d(w_z, p_loc)

        in = 0
        do i = 1, n_loc
          if(.not. is_nonzero(i)) cycle
          in = in + 1
          map(in) = i
          do k = 1, m
            A_z(k,in) = A(k,i)
          enddo
        enddo

        !> 正規方程式 (A_z^T A_z) w_z = A_z^T b の分散求解
        call monolis_scalapack_gram_gesv_R(m, p_loc, A_z, b, w_z, comm, scalapack_comm)

        !> 負成分の有無の全体判定
        neg = 0
        do k = 1, p_loc
          if(w_z(k) < 0.0d0) neg = 1
        enddo
        call monolis_allreduce_I1(neg, monolis_mpi_max, comm)

        if(neg == 0)then
          !> x を更新し、零成分の active 指定を解除する
          do k = 1, p_loc
            x(map(k)) = w_z(k)
            if(w_z(k) <= 0.0d0) is_nonzero(map(k)) = .false.
          enddo
        else
          inner_iter = inner_iter + 1
          if(inner_iter > max_iter)then
            do k = 1, p_loc
              x(map(k)) = max(w_z(k), 0.0d0)
            enddo
            neg = 0
          else
            !> 実行可能領域内に留まる最大ステップ幅の計算
            alpha = huge(0.0d0)
            do k = 1, p_loc
              if(w_z(k) >= 0.0d0) cycle
              t = x(map(k))/(x(map(k)) - w_z(k))
              if(t < alpha) alpha = t
            enddo
            call monolis_allreduce_R1(alpha, monolis_mpi_min, comm)

            !> x の更新と微小成分の active 指定解除
            do k = 1, p_loc
              x(map(k)) = x(map(k)) + alpha*(w_z(k) - x(map(k)))
              if(x(map(k)) <= drop_tol)then
                is_nonzero(map(k)) = .false.
                x(map(k)) = 0.0d0
              endif
            enddo
          endif
        endif

        call monolis_dealloc_I_1d(map)
        call monolis_dealloc_R_2d(A_z)
        call monolis_dealloc_R_1d(w_z)

        if(neg == 0) exit
      enddo

      !> 残差を計算する（全ランクの寄与を集約）
      Ax = 0.0d0
      if(n_loc > 0)then
        call dgemv("N", m, n_loc, 1.0d0, A, m, x, 1, 0.0d0, Ax, 1)
      endif
      call monolis_allreduce_R(m, Ax, monolis_mpi_sum, comm)
      r = b - Ax
      call monolis_get_l2_norm_R(m, r, r_norm)

      if(p_global == n_global)then
        if(.not. is_sparse) is_converge = .true.
        exit
      endif
      if(is_sparse)then
        if(r_norm/r0_norm < tol_outer)then
          is_converge = .true.
          exit
        endif
      endif
    enddo

    !> 残差の計算
    Ax = 0.0d0
    if(n_loc > 0)then
      call dgemv("N", m, n_loc, 1.0d0, A, m, x, 1, 0.0d0, Ax, 1)
    endif
    call monolis_allreduce_R(m, Ax, monolis_mpi_sum, comm)
    call monolis_get_l2_norm_R(m, Ax - b, residual)

    if(.not. is_converge)then
      !call monolis_std_warning_string("monolis_optimize_nnls_LH_main")
      !call monolis_std_warning_string("Residual is not less than tolerance")
      !call monolis_std_error_stop()
    endif

    call monolis_scalapack_comm_finalize(scalapack_comm)

    call monolis_dealloc_R_1d(r)
    call monolis_dealloc_R_1d(s)
    call monolis_dealloc_R_1d(Ax)
    call monolis_dealloc_L_1d(is_nonzero)
  end subroutine monolis_optimize_nnls_LH_main

end module mod_monolis_opt_nnls