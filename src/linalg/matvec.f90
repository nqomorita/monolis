!> 疎行列ベクトル積関数群
module mod_monolis_matvec
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  implicit none

contains

  !> @ingroup linalg
  !> 疎行列ベクトル積（実数型）
  subroutine monolis_matvec_product_R(monolis, X, Y)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> 右辺ベクトル
    real(kdouble) :: X(:)
    !> 結果ベクトル
    real(kdouble) :: Y(:)
    real(kdouble) :: tspmv, tcomm

    call monolis_matvec_product_main_R(monolis%COM, monolis%MAT, X, Y, tspmv, tcomm)
  end subroutine monolis_matvec_product_R

  !> @ingroup linalg
  !> 疎行列ベクトル積（複素数型）
  subroutine monolis_matvec_product_C(monolis, X, Y)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> 右辺ベクトル
    complex(kdouble) :: X(:)
    !> 結果ベクトル
    complex(kdouble) :: Y(:)
    real(kdouble) :: tspmv, tcomm

    call monolis_matvec_product_main_C(monolis%COM, monolis%MAT, X, Y, tspmv, tcomm)
  end subroutine monolis_matvec_product_C

  !> @ingroup dev_linalg
  !> 疎行列ベクトル積（実数型、メイン関数）
  subroutine monolis_matvec_product_main_R(monoCOM, monoMAT, X, Y, tspmv, tcomm)
    implicit none
    !> monolis 構造体
    type(monolis_com) :: monoCOM
    !> monolis 構造体
    type(monolis_mat) :: monoMAT
    !> 右辺ベクトル
    real(kdouble) :: X(:)
    !> 結果ベクトル
    real(kdouble) :: Y(:)
    real(kdouble) :: t1, t2
    real(kdouble) :: tspmv, tcomm

#ifdef DEBUG
    call monolis_std_debug_log_header("monolis_matvec_product_main_R")
#endif
    t1 = monolis_get_time()

    call monolis_update_R(monoCOM, monoMAT%NDOF, X, tcomm)

    if(monoMAT%NDOF == 3)then
      call monolis_matvec_33_R(monoCOM, monoMAT, X, Y)
    elseif(monoMAT%NDOF == 1)then
      call monolis_matvec_11_R(monoCOM, monoMAT, X, Y)
    else
      call monolis_matvec_nn_R(monoCOM, monoMAT, X, Y, monoMAT%NDOF)
    endif

    t2 = monolis_get_time()
    tspmv = tspmv + t2 - t1
  end subroutine monolis_matvec_product_main_R

  !> @ingroup dev_linalg
  !> 疎行列ベクトル積（複素数型、メイン関数）
  subroutine monolis_matvec_product_main_C(monoCOM, monoMAT, X, Y, tspmv, tcomm)
    implicit none
    !> monolis 構造体
    type(monolis_com) :: monoCOM
    !> monolis 構造体
    type(monolis_mat) :: monoMAT
    !> 右辺ベクトル
    complex(kdouble) :: X(:)
    !> 結果ベクトル
    complex(kdouble) :: Y(:)
    real(kdouble) :: t1, t2
    real(kdouble) :: tspmv, tcomm

#ifdef DEBUG
    call monolis_std_debug_log_header("monolis_matvec_product_main_C")
#endif
    t1 = monolis_get_time()

    call monolis_update_C(monoCOM, monoMAT%NDOF, X, tcomm)

    if(monoMAT%NDOF == 3)then
      call monolis_matvec_33_C(monoCOM, monoMAT, X, Y)
    elseif(monoMAT%NDOF == 1)then
      call monolis_matvec_11_C(monoCOM, monoMAT, X, Y)
    else
      call monolis_matvec_nn_C(monoCOM, monoMAT, X, Y, monoMAT%NDOF)
    endif

    t2 = monolis_get_time()
    tspmv = tspmv + t2 - t1
  end subroutine monolis_matvec_product_main_C

  !> @ingroup dev_linalg
  !> 疎行列ベクトル積（実数型、nxn ブロック）
  subroutine monolis_matvec_nn_R(monoCOM, monoMAT, X, Y, NDOF)
    implicit none
    !> monolis 構造体
    type(monolis_com) :: monoCOM
    !> monolis 構造体
    type(monolis_mat), target :: monoMAT
    !> 右辺ベクトル
    real(kdouble) :: X(:)
    !> 結果ベクトル
    real(kdouble) :: Y(:)
    !> ブロックサイズ
    integer(kint) :: NDOF
    integer(kint) :: i, j, k, l, in, N, NDOF2, jS, jE
    real(kdouble) :: XT(NDOF), YT(NDOF)
    integer(kint), pointer :: index(:), item(:)
    real(kdouble), pointer :: A(:)

    N = monoMAT%N
    NDOF2 = NDOF*NDOF
    A => monoMAT%R%A
    index => monoMAT%CSR%index
    item  => monoMAT%CSR%item

!$omp parallel default(none) &
!$omp & shared(A, Y, X, index, item) &
!$omp & firstprivate(N, NDOF, NDOF2) &
!$omp & private(YT, XT, i, j, k, l, jS, jE, in)
!$omp do
    do i = 1, N
      YT = 0.0d0
      jS = index(i) + 1
      jE = index(i + 1)
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
!$omp end do
!$omp end parallel
  end subroutine monolis_matvec_nn_R

  !> @ingroup dev_linalg
  !> 疎行列ベクトル積（複素数、nxn ブロック）
  subroutine monolis_matvec_nn_C(monoCOM, monoMAT, X, Y, NDOF)
    implicit none
    !> monolis 構造体
    type(monolis_com) :: monoCOM
    !> monolis 構造体
    type(monolis_mat), target :: monoMAT
    !> 右辺ベクトル
    complex(kdouble) :: X(:)
    !> 結果ベクトル
    complex(kdouble) :: Y(:)
    !> ブロックサイズ
    integer(kint) :: NDOF
    integer(kint) :: i, j, k, l, in, N, NDOF2, jS, jE
    complex(kdouble) :: XT(NDOF), YT(NDOF)
    integer(kint), pointer :: index(:), item(:)
    complex(kdouble), pointer :: A(:)

    N = monoMAT%N
    NDOF2 = NDOF*NDOF
    A => monoMAT%C%A
    index => monoMAT%CSR%index
    item  => monoMAT%CSR%item

!$omp parallel default(none) &
!$omp & shared(A, Y, X, index, item) &
!$omp & firstprivate(N, NDOF, NDOF2) &
!$omp & private(YT, XT, i, j, k, l, jS, jE, in)
!$omp do
    do i = 1, N
      YT = 0.0d0
      jS = index(i) + 1
      jE = index(i + 1)
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
!$omp end do
!$omp end parallel
  end subroutine monolis_matvec_nn_C

  !> @ingroup dev_linalg
  !> 疎行列ベクトル積（実数型、1x1 ブロック）
  subroutine monolis_matvec_11_R(monoCOM, monoMAT, X, Y)
    implicit none
    !> monolis 構造体
    type(monolis_com) :: monoCOM
    !> monolis 構造体
    type(monolis_mat), target :: monoMAT
    !> 右辺ベクトル
    real(kdouble) :: X(:)
    !> 結果ベクトル
    real(kdouble) :: Y(:)
    integer(kint) :: i, j, in, N, jS, jE
    integer(kint), pointer :: index(:), item(:)
    real(kdouble) :: Y1
    real(kdouble), pointer :: A(:)

    N = monoMAT%N
    A => monoMAT%R%A
    index => monoMAT%CSR%index
    item  => monoMAT%CSR%item

!$omp parallel default(none) &
!$omp & shared(A, Y, X, index, item) &
!$omp & firstprivate(N) &
!$omp & private(Y1, i, j, jS, jE, in)
!$omp do
    do i = 1, N
      Y1 = 0.0d0
      jS = index(i) + 1
      jE = index(i + 1)
      do j = jS, jE
        in = item(j)
        Y1 = Y1 + A(j)*X(in)
      enddo
      Y(i) = Y1
    enddo
!$omp end do
!$omp end parallel
  end subroutine monolis_matvec_11_R

  !> @ingroup dev_linalg
  !> 疎行列ベクトル積（複素数型、1x1 ブロック）
  subroutine monolis_matvec_11_C(monoCOM, monoMAT, X, Y)
    implicit none
    !> monolis 構造体
    type(monolis_com) :: monoCOM
    !> monolis 構造体
    type(monolis_mat), target :: monoMAT
    !> 右辺ベクトル
    complex(kdouble) :: X(:)
    !> 結果ベクトル
    complex(kdouble) :: Y(:)
    integer(kint) :: i, j, in, N, jS, jE
    integer(kint), pointer :: index(:), item(:)
    complex(kdouble) :: Y1
    complex(kdouble), pointer :: A(:)

    N = monoMAT%N
    A => monoMAT%C%A
    index => monoMAT%CSR%index
    item  => monoMAT%CSR%item

!$omp parallel default(none) &
!$omp & shared(A, Y, X, index, item) &
!$omp & firstprivate(N) &
!$omp & private(Y1, i, j, jS, jE, in)
!$omp do
    do i = 1, N
      Y1 = 0.0d0
      jS = index(i) + 1
      jE = index(i + 1)
      do j = jS, jE
        in = item(j)
        Y1 = Y1 + A(j)*X(in)
      enddo
      Y(i) = Y1
    enddo
!$omp end do
!$omp end parallel
  end subroutine monolis_matvec_11_C

  !> @ingroup dev_linalg
  !> 疎行列ベクトル積（実数型、3x3 ブロック）
  subroutine monolis_matvec_33_R(monoCOM, monoMAT, X, Y)
    implicit none
    !> monolis 構造体
    type(monolis_com) :: monoCOM
    !> monolis 構造体
    type(monolis_mat), target :: monoMAT
    !> 右辺ベクトル
    real(kdouble) :: X(:)
    !> 結果ベクトル
    real(kdouble) :: Y(:)
    integer(kint) :: i, j, in, N, jS, jE
    integer(kint), pointer :: index(:), item(:)
    real(kdouble) :: X1, X2, X3, Y1, Y2, Y3
    real(kdouble), pointer :: A(:)

    N = monoMAT%N
    A => monoMAT%R%A
    index => monoMAT%CSR%index
    item  => monoMAT%CSR%item

!$omp parallel default(none) &
!$omp & shared(A, Y, X, index, item) &
!$omp & firstprivate(N) &
!$omp & private(Y1, Y2, Y3, X1, X2, X3, i, j, jS, jE, in)
!$omp do
    do i = 1, N
      Y1 = 0.0d0
      Y2 = 0.0d0
      Y3 = 0.0d0
      jS = index(i) + 1
      jE = index(i + 1)
      do j = jS, jE
        in = item(j)
        X1 = X(3*in-2)
        X2 = X(3*in-1)
        X3 = X(3*in  )
        Y1 = Y1 + A(9*j-8)*X1 + A(9*j-7)*X2 + A(9*j-6)*X3
        Y2 = Y2 + A(9*j-5)*X1 + A(9*j-4)*X2 + A(9*j-3)*X3
        Y3 = Y3 + A(9*j-2)*X1 + A(9*j-1)*X2 + A(9*j  )*X3
      enddo
      Y(3*i-2) = Y1
      Y(3*i-1) = Y2
      Y(3*i  ) = Y3
    enddo
!$omp end do
!$omp end parallel
  end subroutine monolis_matvec_33_R

  !> @ingroup dev_linalg
  !> 疎行列ベクトル積（複素数型、3x3 ブロック）
  subroutine monolis_matvec_33_C(monoCOM, monoMAT, X, Y)
    implicit none
    !> monolis 構造体
    type(monolis_com) :: monoCOM
    !> monolis 構造体
    type(monolis_mat), target :: monoMAT
    !> 右辺ベクトル
    complex(kdouble) :: X(:)
    !> 結果ベクトル
    complex(kdouble) :: Y(:)
    integer(kint) :: i, j, in, N, jS, jE
    integer(kint), pointer :: index(:), item(:)
    complex(kdouble) :: X1, X2, X3, Y1, Y2, Y3
    complex(kdouble), pointer :: A(:)

    N = monoMAT%N
    A => monoMAT%C%A
    index => monoMAT%CSR%index
    item  => monoMAT%CSR%item

!$omp parallel default(none) &
!$omp & shared(A, Y, X, index, item) &
!$omp & firstprivate(N) &
!$omp & private(Y1, Y2, Y3, X1, X2, X3, i, j, jS, jE, in)
!$omp do
    do i = 1, N
      Y1 = 0.0d0
      Y2 = 0.0d0
      Y3 = 0.0d0
      jS = index(i) + 1
      jE = index(i + 1)
      do j = jS, jE
        in = item(j)
        X1 = X(3*in-2)
        X2 = X(3*in-1)
        X3 = X(3*in  )
        Y1 = Y1 + A(9*j-8)*X1 + A(9*j-7)*X2 + A(9*j-6)*X3
        Y2 = Y2 + A(9*j-5)*X1 + A(9*j-4)*X2 + A(9*j-3)*X3
        Y3 = Y3 + A(9*j-2)*X1 + A(9*j-1)*X2 + A(9*j  )*X3
      enddo
      Y(3*i-2) = Y1
      Y(3*i-1) = Y2
      Y(3*i  ) = Y3
    enddo
!$omp end do
!$omp end parallel
  end subroutine monolis_matvec_33_C

  !> @ingroup dev_linalg
  !> 残差ベクトルの取得（実数型、メイン関数）
  subroutine monolis_residual_main_R(monoCOM, monoMAT, X, B, R, tspmv, tcomm)
    implicit none
    !> monolis 構造体
    type(monolis_com) :: monoCOM
    !> monolis 構造体
    type(monolis_mat) :: monoMAT
    !> 解ベクトル
    real(kdouble) :: X(:)
    !> 右辺ベクトル
    real(kdouble) :: B(:)
    !> 残差ベクトル
    real(kdouble) :: R(:)
    !> 計算時間
    real(kdouble) :: tspmv
    !> 通信時間
    real(kdouble) :: tcomm
    integer(kint) :: i

    call monolis_matvec_product_main_R(monoCOM, monoMAT, X, R, tspmv, tcomm)

!$omp parallel default(none) &
!$omp & shared(monoMAT, B, R) &
!$omp & private(i)
!$omp do
    do i = 1, monoMAT%N*monoMAT%NDOF
      R(i) = B(i) - R(i)
    enddo
!$omp end do
!$omp end parallel
  end subroutine monolis_residual_main_R

  !> @ingroup dev_linalg
  !> 残差ベクトルの取得（複素数型、メイン関数）
  subroutine monolis_residual_main_C(monoCOM, monoMAT, X, B, R, tspmv, tcomm)
    implicit none
    !> monolis 構造体
    type(monolis_com) :: monoCOM
    !> monolis 構造体
    type(monolis_mat) :: monoMAT
    !> 解ベクトル
    complex(kdouble) :: X(:)
    !> 右辺ベクトル
    complex(kdouble) :: B(:)
    !> 残差ベクトル
    complex(kdouble) :: R(:)
    !> 計算時間
    real(kdouble) :: tspmv
    !> 通信時間
    real(kdouble) :: tcomm
    integer(kint) :: i

    call monolis_matvec_product_main_C(monoCOM, monoMAT, X, R, tspmv, tcomm)

!$omp parallel default(none) &
!$omp & shared(monoMAT, B, R) &
!$omp & private(i)
!$omp do
    do i = 1, monoMAT%N*monoMAT%NDOF
      R(i) = B(i) - R(i)
    enddo
!$omp end do
!$omp end parallel
  end subroutine monolis_residual_main_C
end module mod_monolis_matvec