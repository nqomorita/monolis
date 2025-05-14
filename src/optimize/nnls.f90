!> NNLS 関数
module mod_monolis_opt_nnls
  use mod_monolis_utils
  use mod_monolis_lapack
  implicit none

contains

  !> @ingroup linalg
  !> Non-Negative Least Squares
  subroutine monolis_optimize_parallel_nnls_R_with_sparse_solution(A, b, x, m, n, max_iter, tol, residual, monoCOM)
    implicit none
    !> [in] 行列
    real(kdouble), intent(in) :: A(:,:)
    !> [in] 右辺ベクトル
    real(kdouble), intent(in) :: b(:)
    !> [out] 解ベクトル
    real(kdouble), intent(out) :: x(:)
    !> [in] 行列の大きさ
    integer(kint), intent(in) :: m
    !> [in] 行列の大きさ
    integer(kint), intent(in) :: n
    !> [in] 最大反復回数
    integer(kint), intent(in) :: max_iter
    !> [out] 収束判定閾値
    real(kdouble), intent(in) :: tol
    !> [out] 残差
    real(kdouble), intent(out) :: residual
    !> [in] COM 構造体
    type(monolis_COM), intent(in) :: monoCOM
    integer(kint) :: idx(1), iter, p, i, in
    real(kdouble) :: r0_norm, r_norm, res
    real(kdouble), allocatable :: r(:)
    real(kdouble), allocatable :: s(:)
    real(kdouble), allocatable :: A_z(:,:)
    real(kdouble), allocatable :: Q_z(:,:)
    real(kdouble), allocatable :: R_z(:,:)
    real(kdouble), allocatable :: c(:)
    real(kdouble), allocatable :: w_z(:)
    logical, allocatable :: is_nonzero(:)

    !> メモリの確保
    call monolis_alloc_R_1d(r, m)
    call monolis_alloc_R_1d(s, n)
    call monolis_alloc_L_1d(is_nonzero, n)

    !> 収束判定のためのノルム計算
    call monolis_get_l2_norm_R(m, b, r0_norm)

    is_nonzero = .false.
    r = b
    x = 0.0d0

    do iter =  1, max_iter
      !> 行列 A の転置と残差ベクトルをかける
      s = matmul(transpose(A), r)

      !> 最大値のインデックスを求め、その値を非零に指定する
      idx = maxloc(s, .not. is_nonzero)
      is_nonzero(idx(1)) = .true.

      !> 行列 A から is_nonzero 配列で非零に指定された列要素を取得
      p = 0
      do i = 1, n
        if(is_nonzero(i)) p = p + 1
      enddo

      call monolis_alloc_R_2d(A_z, m, p)
      call monolis_alloc_R_2d(Q_z, m, p)
      call monolis_alloc_R_2d(R_z, p, p)
      call monolis_alloc_R_1d(c, p)
      call monolis_alloc_R_1d(w_z, p)

      call get_column_matrix(m, n, A, A_z, is_nonzero)
      call monolis_lapack_dgeqrf(m, p, A_z, Q_z, R_z)

      c = matmul(transpose(Q_z), b)
      call monolis_optimize_nnls_R(R_z, c, w_z, p, p, max_iter, tol, res)

      !> x を更新
      in = 0
      do i = 1, n
        if(.not. is_nonzero(i)) cycle
        in = in + 1
        x(i) = w_z(in)
      enddo

      !> 解の中で値が負の要素の非零指定を解除する
      in = 0
      do i = 1, n
        if(is_nonzero(i))then
          in = in + 1
          if(w_z(in) <= 0.0d0)then
            is_nonzero(i) = .false.
          endif
        endif
      enddo

      !> 残差を計算する
      r = b - matmul(A_z, w_z)
      call monolis_get_l2_norm_R(m, r, r_norm)

      call monolis_dealloc_R_2d(A_z)
      call monolis_dealloc_R_2d(Q_z)
      call monolis_dealloc_R_2d(R_z)
      call monolis_dealloc_R_1d(w_z)
      call monolis_dealloc_R_1d(c)

!write(*,*)"r_norm/r0_norm", r_norm/r0_norm
      if(all(is_nonzero)) exit
      if(r_norm/r0_norm < tol) exit
    enddo

    call monolis_get_l2_norm_R(m, matmul(A, x) - b, residual)

    call monolis_dealloc_R_1d(r)
    call monolis_dealloc_R_1d(s)
    call monolis_dealloc_L_1d(is_nonzero)
  end subroutine monolis_optimize_parallel_nnls_R_with_sparse_solution

  !> @ingroup linalg
  !> Non-Negative Least Squares
  subroutine monolis_optimize_nnls_R_with_sparse_solution(A, b, x, m, n, max_iter, tol, residual)
    implicit none
    !> [in] 行列
    real(kdouble), intent(in) :: A(:,:)
    !> [in] 右辺ベクトル
    real(kdouble), intent(in) :: b(:)
    !> [out] 解ベクトル
    real(kdouble), intent(out) :: x(:)
    !> [in] 行列の大きさ
    integer(kint), intent(in) :: m
    !> [in] 行列の大きさ
    integer(kint), intent(in) :: n
    !> [in] 最大反復回数
    integer(kint), intent(in) :: max_iter
    !> [out] 収束判定閾値
    real(kdouble), intent(in) :: tol
    !> [out] 残差
    real(kdouble), intent(out) :: residual
    integer(kint) :: idx(1), iter, p, i, in
    real(kdouble) :: r0_norm, r_norm, res
    real(kdouble), allocatable :: r(:)
    real(kdouble), allocatable :: s(:)
    real(kdouble), allocatable :: A_z(:,:)
    real(kdouble), allocatable :: Q_z(:,:)
    real(kdouble), allocatable :: R_z(:,:)
    real(kdouble), allocatable :: c(:)
    real(kdouble), allocatable :: w_z(:)
    logical, allocatable :: is_nonzero(:)

    !> メモリの確保
    call monolis_alloc_R_1d(r, m)
    call monolis_alloc_R_1d(s, n)
    call monolis_alloc_L_1d(is_nonzero, n)

    !> 収束判定のためのノルム計算
    call monolis_get_l2_norm_R(m, b, r0_norm)

    is_nonzero = .false.
    r = b
    x = 0.0d0

    do iter =  1, max_iter
      !> 行列 A の転置と残差ベクトルをかける
      s = matmul(transpose(A), r)

      !> 最大値のインデックスを求め、その値を非零に指定する
      idx = maxloc(s, .not. is_nonzero)
      is_nonzero(idx(1)) = .true.

      !> 行列 A から is_nonzero 配列で非零に指定された列要素を取得
      p = 0
      do i = 1, n
        if(is_nonzero(i)) p = p + 1
      enddo

      call monolis_alloc_R_2d(A_z, m, p)
      call monolis_alloc_R_2d(Q_z, m, p)
      call monolis_alloc_R_2d(R_z, p, p)
      call monolis_alloc_R_1d(c, p)
      call monolis_alloc_R_1d(w_z, p)

      call get_column_matrix(m, n, A, A_z, is_nonzero)
      call monolis_lapack_dgeqrf(m, p, A_z, Q_z, R_z)

      c = matmul(transpose(Q_z), b)
      call monolis_optimize_nnls_R(R_z, c, w_z, p, p, max_iter, tol, res)

      !> x を更新
      in = 0
      do i = 1, n
        if(.not. is_nonzero(i)) cycle
        in = in + 1
        x(i) = w_z(in)
      enddo

      !> 解の中で値が負の要素の非零指定を解除する
      in = 0
      do i = 1, n
        if(is_nonzero(i))then
          in = in + 1
          if(w_z(in) <= 0.0d0)then
            is_nonzero(i) = .false.
          endif
        endif
      enddo

      !> 残差を計算する
      r = b - matmul(A_z, w_z)
      call monolis_get_l2_norm_R(m, r, r_norm)

      call monolis_dealloc_R_2d(A_z)
      call monolis_dealloc_R_2d(Q_z)
      call monolis_dealloc_R_2d(R_z)
      call monolis_dealloc_R_1d(w_z)
      call monolis_dealloc_R_1d(c)

!write(*,*)"r_norm/r0_norm", r_norm/r0_norm
      if(all(is_nonzero)) exit
      if(r_norm/r0_norm < tol) exit
    enddo

    call monolis_get_l2_norm_R(m, matmul(A, x) - b, residual)

    call monolis_dealloc_R_1d(r)
    call monolis_dealloc_R_1d(s)
    call monolis_dealloc_L_1d(is_nonzero)
  end subroutine monolis_optimize_nnls_R_with_sparse_solution

  subroutine get_column_matrix(m, n, A, A_z, is_nonzero)
    implicit none
    real(kdouble) :: A(:,:)
    real(kdouble) :: A_z(:,:)
    logical :: is_nonzero(:)
    integer(kint) :: m, n
    integer(kint) :: i, in, j

    in = 0
    do i = 1, n
      if(.not. is_nonzero(i)) cycle
      in = in + 1
      do j = 1, m
        A_z(j,in) = A(j,i)
      enddo
    enddo
  end subroutine get_column_matrix

  !> @ingroup linalg
  !> Non-Negative Least Squares
  subroutine monolis_optimize_nnls_R(A, b, x, m, n, max_iter, tol, residual)
    implicit none
    !> [in] 行列
    real(kdouble), intent(in) :: A(:,:)
    !> [in] 右辺ベクトル
    real(kdouble), intent(in) :: b(:)
    !> [out] 解ベクトル
    real(kdouble), intent(out) :: x(:)
    !> [in] 行列の大きさ
    integer(kint), intent(in) :: m
    !> [in] 行列の大きさ
    integer(kint), intent(in) :: n
    !> [in] 最大反復回数
    integer(kint), intent(in) :: max_iter
    !> [out] 収束判定閾値
    real(kdouble), intent(in) :: tol
    !> [out] 残差
    real(kdouble), intent(out) :: residual
    integer(kint) :: idx(1), iter
    real(kdouble) :: alpha
    logical :: is_converge_inner, is_converge, is_all_false
    real(kdouble), allocatable :: AtA(:,:)
    real(kdouble), allocatable :: Atb(:)
    real(kdouble), allocatable :: s(:)
    real(kdouble), allocatable :: w(:)
    logical, allocatable :: P(:)

    call monolis_alloc_R_2d(AtA, n, n)
    call monolis_alloc_R_1d(Atb, n)
    call monolis_alloc_R_1d(s, n)
    call monolis_alloc_R_1d(w, n)
    call monolis_alloc_L_1d(P, n)

    AtA = matmul(transpose(A), A)
    Atb = matmul(transpose(A), b)
    w = Atb
    iter = 0
    x = 0.0d0

    do
      !> Get the "most" active coeff index and move to inactive set
      idx = maxloc(w, .not. P)
      if(idx(1) == 0) exit
      P(idx(1)) = .true.

      !> Iteration solution
      s = 0.0d0
      call monolis_lapack_dsysv_with_select(n, AtA, Atb, s, P)

      !> Inner loop
      aa:do
        call check_tolerance_inner(n, s, P, is_converge_inner)
        if(is_converge_inner) exit aa

        iter = iter + 1
        if(iter > max_iter) exit

        call update_candidate(n, s, P)
        call get_minimum_alpha(n, x, s, P, alpha)

        x = (1.0d0 - alpha)*x
        x = x + alpha*s

        call check_is_all(n, x, tol, P, is_all_false)
        if(is_all_false)then
          s = 0.0d0
          exit
        endif

        call monolis_lapack_dsysv_with_select(n, AtA, Atb, s, P)
      enddo aa
      x = s
      w = Atb - matmul(AtA, x)

      call check_tolerance_outer(n, w, P, tol, is_converge)
      if(is_converge) exit
    enddo

    call monolis_get_l2_norm_R(m, matmul(A, x) - b, residual)

    call monolis_dealloc_R_2d(AtA)
    call monolis_dealloc_R_1d(Atb)
    call monolis_dealloc_R_1d(s)
    call monolis_dealloc_R_1d(w)
    call monolis_dealloc_L_1d(P)
  end subroutine monolis_optimize_nnls_R

  subroutine get_minimum_alpha(n, x, s, P, alpha)
    implicit none
    integer(kint) :: n
    real(kdouble) :: x(:)
    real(kdouble) :: s(:)
    logical :: P(:)
    real(kdouble) :: alpha
    real(kdouble) :: temp
    integer(kint) :: i

    alpha = huge(kdouble)
    do i = 1, n
      if(.not. P(i)) cycle
      temp = x(i)/(x(i) - s(i))
      if(temp < alpha) alpha = temp
    enddo
  end subroutine get_minimum_alpha

  subroutine check_tolerance_outer(n, w, P, tol, is_converge)
    implicit none
    integer(kint) :: n
    real(kdouble) :: w(:)
    logical :: P(:)
    real(kdouble) :: tol
    logical :: is_converge
    integer(kint) :: i

    is_converge = .true.
    do i = 1, n
      if(P(i)) cycle
      if(w(i) > tol) is_converge = .false.
    enddo
  end subroutine check_tolerance_outer

  subroutine check_tolerance_inner(n, s, P, is_converge_inner)
    implicit none
    integer(kint) :: n
    real(kdouble) :: s(:)
    logical :: P(:)
    logical :: is_converge_inner
    integer(kint) :: i

    is_converge_inner = .true.
    do i = 1, n
      if(.not. P(i)) cycle
      if(s(i) < 0.0d0) is_converge_inner = .false.
    enddo
  end subroutine check_tolerance_inner

  subroutine check_is_all(n, x, tol, P, is_all_false)
    implicit none
    integer(kint) :: n
    real(kdouble) :: x(:)
    real(kdouble) :: tol
    logical :: P(:)
    logical :: is_all_false
    integer(kint) :: i

    do i = 1, n
      if(x(i) < tol)then
        P(i) = .false.
      else
        P(i) = .true.
      endif
    enddo

    is_all_false = .true.
    do i = 1, n
      if(P(i))then
        is_all_false = .false.
        return
      endif
    enddo
  end subroutine check_is_all

  subroutine update_candidate(n, s, P)
    implicit none
    integer(kint) :: n
    real(kdouble) :: s(:)
    logical :: P(:)
    integer(kint) :: i

    do i = 1, n
      if(.not. P(i)) cycle
      if(s(i) > 0.0d0) P(i) = .false.
    enddo
  end subroutine update_candidate

  subroutine check_all_positive(n, s, P, is_all_positve)
    implicit none
    integer(kint) :: n
    real(kdouble) :: s(:)
    logical :: P(:)
    logical :: is_all_positve
    integer(kint) :: i

    is_all_positve = .true.
    do i = 1, n
      if(.not. P(i)) cycle
      if(s(i) > 0.0d0) is_all_positve = .false.
    enddo
  end subroutine check_all_positive

end module mod_monolis_opt_nnls