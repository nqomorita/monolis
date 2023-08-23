module mod_monolis_matmat
  use mod_monolis_prm
  use mod_monolis_mat
  implicit none

contains

  !> @ingroup linalg
  !> ベクトル内積
  subroutine monolis_matmat(monoCOM, monoMAT, X, Y, tspmv, tcomm)
    implicit none
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: monoMAT
    !> [in,out] 右辺ベクトル
    real(kdouble), intent(inout) :: X(:)
    !> [out] 結果ベクトル
    real(kdouble), intent(out) :: Y(:)
    !> [in,out] 計算時間
    real(kdouble), intent(inout) :: tspmv
    !> [in,out] 通信時間
    real(kdouble), intent(inout) :: tcomm
    real(kdouble) :: t1, t2

#ifdef DEBUG
    call monolis_std_debug_log_header("monolis_matmat")
#endif
    t1 = monolis_get_time()

    !call monolis_update_pre_R(monoCOM, monoMAT%NDOF, X, tcomm)

    call monolis_matmat_nn(monoCOM, monoMAT, X, Y, monoMAT%NDOF)

    !call monolis_update_post_R(monoCOM, monoMAT%NDOF, X, tcomm)

    t2 = monolis_get_time()
    tspmv = tspmv + t2 - t1
  end subroutine monolis_matmat

  !> @ingroup linalg
  !> ベクトル内積
  subroutine monolis_matmat_nn(monoCOM, monoMAT, X, Y, NDOF)
    implicit none
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: monoMAT
    !> [in,out] 右辺ベクトル
    real(kdouble), intent(inout) :: X(:)
    !> [out] 結果ベクトル
    real(kdouble), intent(out) :: Y(:)
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: NDOF
    integer(kint) :: i, j, k, l, in, N, NDOF2, jS, jE
    integer(kint), pointer :: index(:), item(:)
    real(kdouble) :: XT(NDOF), YT(NDOF)
    real(kdouble), pointer :: A(:)

    N = monoMAT%N
    NDOF2 = NDOF*NDOF
    A => monoMAT%A
    index => monoMAT%index
    item  => monoMAT%item

    do i = 1, N
      YT = 0.0d0
      jS = index(i-1) + 1
      jE = index(i  )
      do j = jS, jE
        in = item(j)
        do k = 1, NDOF
          XT(k) = X(NDOF*(in-1)+k)
        enddo
        do k = 1, NDOF
          do l = 1, NDOF
            YT(k) = YT(k) + A(NDOF2*(j-1)+NDOF*(k-1)+l) * XT(l)
          enddo
        enddo
      enddo
      do k = 1, NDOF
        Y(NDOF*(i-1)+k) = YT(k)
      enddo
    enddo
  end subroutine monolis_matmat_nn

end module mod_monolis_matmat