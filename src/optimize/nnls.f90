!> NNLS 関数
module mod_monolis_opt_nnls
  use mod_monolis_utils
  use mod_monolis_lapack
  implicit none

contains

  !> @ingroup linalg
  !> Non-Negative Least Squares
  subroutine monolis_optimize_nnls(A, b, x, m, n, max_iter, tol, residual)
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
    logical :: is_all_positve, is_converge_inner, is_converge, is_all_false
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
        call check_tolerance_inner(n, s, tol, P, is_converge_inner)
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
  end subroutine monolis_optimize_nnls

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

  subroutine check_tolerance_inner(n, s, tol, P, is_converge_inner)
    implicit none
    integer(kint) :: n
    real(kdouble) :: s(:)
    real(kdouble) :: tol
    logical :: P(:)
    logical :: is_converge_inner
    integer(kint) :: i

    is_converge_inner = .true.
    do i = 1, n
      if(.not. P(i)) cycle
      if(s(i) < 0) is_converge_inner = .false.
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
      if(s(i) > 0) P(i) = .false.
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
      if(s(i) > 0) is_all_positve = .false.
    enddo
  end subroutine check_all_positive

end module mod_monolis_opt_nnls