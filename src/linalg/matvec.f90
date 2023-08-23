!> 疎行列ベクトル積関数群
module mod_monolis_matvec
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  implicit none

contains

  !> @ingroup linalg
  !> 疎行列ベクトル積（実数型）
  subroutine monolis_matvec_product_R(monolis, monoCOM, X, Y)
    implicit none
    !> [in] monolis 構造体
    type(monolis_structure), intent(in) :: monolis
    !> [in] COM 構造体
    type(monolis_COM), intent(in) :: monoCOM
    !> [in,out] 右辺ベクトル
    real(kdouble), intent(inout) :: X(:)
    !> [out] 結果ベクトル
    real(kdouble), intent(out) :: Y(:)
    real(kdouble) :: tspmv, tcomm

    call monolis_std_debug_log_header("monolis_matvec_product_R")

    call monolis_matvec_product_main_R(monoCOM, monolis%MAT, X, Y, tspmv, tcomm)

    call monolis_mpi_update_R(monoCOM, monolis%MAT%NDOF, Y, tcomm)
  end subroutine monolis_matvec_product_R

  !> @ingroup linalg
  !> 疎行列ベクトル積（複素数型）
  subroutine monolis_matvec_product_C(monolis, monoCOM, X, Y)
    implicit none
    !> [in] monolis 構造体
    type(monolis_structure), intent(in) :: monolis
    !> [in] COM 構造体
    type(monolis_COM), intent(in) :: monoCOM
    !> [in,out] 右辺ベクトル
    complex(kdouble), intent(inout) :: X(:)
    !> [out] 結果ベクトル
    complex(kdouble), intent(out) :: Y(:)
    real(kdouble) :: tspmv, tcomm

    call monolis_std_debug_log_header("monolis_matvec_product_C")

    call monolis_matvec_product_main_C(monoCOM, monolis%MAT, X, Y, tspmv, tcomm)

    call monolis_mpi_update_C(monoCOM, monolis%MAT%NDOF, Y, tcomm)
  end subroutine monolis_matvec_product_C

  !> @ingroup dev_linalg
  !> 疎行列ベクトル積（実数型、メイン関数）
  subroutine monolis_matvec_product_main_R(monoCOM, monoMAT, X, Y, tspmv, tcomm)
    implicit none
    !> [in] monolis 構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] monolis 構造体
    type(monolis_mat), intent(in) :: monoMAT
    !> [in,out] 右辺ベクトル
    real(kdouble), intent(inout) :: X(:)
    !> [out] 結果ベクトル
    real(kdouble), intent(out) :: Y(:)
    !> [in,out] 計算時間
    real(kdouble) :: tspmv
    !> [in,out] 通信時間
    real(kdouble) :: tcomm
    real(kdouble) :: t1, t2

    call monolis_std_debug_log_header("monolis_matvec_product_main_R")

    t1 = monolis_get_time()

    call monolis_mpi_update_R(monoCOM, monoMAT%NDOF, X, tcomm)

    if(monoMAT%NDOF == 3)then
      call monolis_matvec_33_R(monoMAT%N, monoMAT%CSR%index, monoMAT%CSR%item, monoMAT%R%A, X, Y)
    elseif(monoMAT%NDOF == 1)then
      call monolis_matvec_11_R(monoMAT%N, monoMAT%CSR%index, monoMAT%CSR%item, monoMAT%R%A, X, Y)
    else
      call monolis_matvec_nn_R(monoMAT%N, monoMAT%CSR%index, monoMAT%CSR%item, monoMAT%R%A, X, Y, monoMAT%NDOF)
    endif

    t2 = monolis_get_time()
    tspmv = tspmv + t2 - t1
  end subroutine monolis_matvec_product_main_R

  !> @ingroup dev_linalg
  !> 疎行列ベクトル積（複素数型、メイン関数）
  subroutine monolis_matvec_product_main_C(monoCOM, monoMAT, X, Y, tspmv, tcomm)
    implicit none
    !> [in] monolis 構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] monolis 構造体
    type(monolis_mat), intent(in) :: monoMAT
    !> [in,out] 右辺ベクトル
    complex(kdouble), intent(inout) :: X(:)
    !> [out] 結果ベクトル
    complex(kdouble), intent(out) :: Y(:)
    !> [in,out] 計算時間
    real(kdouble) :: tspmv
    !> [in,out] 通信時間
    real(kdouble) :: tcomm
    real(kdouble) :: t1, t2

    call monolis_std_debug_log_header("monolis_matvec_product_main_C")

    t1 = monolis_get_time()

    call monolis_mpi_update_C(monoCOM, monoMAT%NDOF, X, tcomm)

    if(monoMAT%NDOF == 3)then
      call monolis_matvec_33_C(monoMAT%N, monoMAT%CSR%index, monoMAT%CSR%item, monoMAT%C%A, X, Y)
    elseif(monoMAT%NDOF == 1)then
      call monolis_matvec_11_C(monoMAT%N, monoMAT%CSR%index, monoMAT%CSR%item, monoMAT%C%A, X, Y)
    else
      call monolis_matvec_nn_C(monoMAT%N, monoMAT%CSR%index, monoMAT%CSR%item, monoMAT%C%A, X, Y, monoMAT%NDOF)
    endif

    t2 = monolis_get_time()
    tspmv = tspmv + t2 - t1
  end subroutine monolis_matvec_product_main_C

  !> @ingroup dev_linalg
  !> 疎行列ベクトル積（実数型、nxn ブロック）
  subroutine monolis_matvec_nn_R(N, index, item, A, X, Y, NDOF)
    implicit none
    !> [in] 行列サイズ
    integer(kint), intent(in) :: N
    !> [in] index 配列
    integer(kint), intent(in) :: index(:)
    !> [in] item 配列
    integer(kint), intent(in) :: item(:)
    !> [in] 行列成分
    real(kdouble), intent(in) :: A(:)
    !> [in] 右辺ベクトル
    real(kdouble), intent(in) :: X(:)
    !> [out] 結果ベクトル
    real(kdouble), intent(out) :: Y(:)
    !> [in] ブロックサイズ
    integer(kint), intent(in) :: NDOF
    integer(kint) :: i, j, k, l, in, NDOF2, jS, jE
    real(kdouble) :: XT(NDOF), YT(NDOF)

    NDOF2 = NDOF*NDOF

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
  subroutine monolis_matvec_nn_C(N, index, item, A, X, Y, NDOF)
    implicit none
    !> [in] 行列サイズ
    integer(kint), intent(in) :: N
    !> [in] index 配列
    integer(kint), intent(in) :: index(:)
    !> [in] item 配列
    integer(kint), intent(in) :: item(:)
    !> [in] 行列成分
    complex(kdouble), intent(in) :: A(:)
    !> [in] 右辺ベクトル
    complex(kdouble), intent(in) :: X(:)
    !> [out] 結果ベクトル
    complex(kdouble), intent(out) :: Y(:)
    !> [in] ブロックサイズ
    integer(kint), intent(in) :: NDOF
    integer(kint) :: i, j, k, l, in, NDOF2, jS, jE
    complex(kdouble) :: XT(NDOF), YT(NDOF)

    NDOF2 = NDOF*NDOF

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
  subroutine monolis_matvec_11_R(N, index, item, A, X, Y)
    implicit none
    !> [in] 行列サイズ
    integer(kint), intent(in) :: N
    !> [in] index 配列
    integer(kint), intent(in) :: index(:)
    !> [in] item 配列
    integer(kint), intent(in) :: item(:)
    !> [in] 行列成分
    real(kdouble), intent(in) :: A(:)
    !> [in] 右辺ベクトル
    real(kdouble), intent(in) :: X(:)
    !> [out] 結果ベクトル
    real(kdouble), intent(out) :: Y(:)
    integer(kint) :: i, j, in, jS, jE
    real(kdouble) :: Y1

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
  subroutine monolis_matvec_11_C(N, index, item, A, X, Y)
    implicit none
    !> [in] 行列サイズ
    integer(kint), intent(in) :: N
    !> [in] index 配列
    integer(kint), intent(in) :: index(:)
    !> [in] item 配列
    integer(kint), intent(in) :: item(:)
    !> [in] 行列成分
    complex(kdouble), intent(in) :: A(:)
    !> [in] 右辺ベクトル
    complex(kdouble), intent(in) :: X(:)
    !> [out] 結果ベクトル
    complex(kdouble), intent(out) :: Y(:)
    integer(kint) :: i, j, in, jS, jE
    complex(kdouble) :: Y1

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
  subroutine monolis_matvec_33_R(N, index, item, A, X, Y)
    implicit none
    !> [in] 行列サイズ
    integer(kint), intent(in) :: N
    !> [in] index 配列
    integer(kint), intent(in) :: index(:)
    !> [in] item 配列
    integer(kint), intent(in) :: item(:)
    !> [in] 行列成分
    real(kdouble), intent(in) :: A(:)
    !> [in] 右辺ベクトル
    real(kdouble), intent(in) :: X(:)
    !> [out] 結果ベクトル
    real(kdouble), intent(out) :: Y(:)
    integer(kint) :: i, j, in, jS, jE
    real(kdouble) :: X1, X2, X3, Y1, Y2, Y3

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
  subroutine monolis_matvec_33_C(N, index, item, A, X, Y)
    implicit none
    !> [in] 行列サイズ
    integer(kint), intent(in) :: N
    !> [in] index 配列
    integer(kint), intent(in) :: index(:)
    !> [in] item 配列
    integer(kint), intent(in) :: item(:)
    !> [in] 行列成分
    complex(kdouble), intent(in) :: A(:)
    !> [in] 右辺ベクトル
    complex(kdouble), intent(in) :: X(:)
    !> [out] 結果ベクトル
    complex(kdouble), intent(out) :: Y(:)
    integer(kint) :: i, j, in, jS, jE
    complex(kdouble) :: X1, X2, X3, Y1, Y2, Y3

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
    !> [in] monolis 構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] monolis 構造体
    type(monolis_mat), intent(in) :: monoMAT
    !> [in,out] 解ベクトル
    real(kdouble), intent(inout) :: X(:)
    !> [in] 右辺ベクトル
    real(kdouble), intent(in) :: B(:)
    !> [out] 残差ベクトル
    real(kdouble), intent(out) :: R(:)
    !> [in,out] 計算時間
    real(kdouble), intent(inout) :: tspmv
    !> [in,out] 通信時間
    real(kdouble), intent(inout) :: tcomm
    integer(kint) :: i

    call monolis_std_debug_log_header("monolis_residual_main_R")

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
    !> [in] monolis 構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] monolis 構造体
    type(monolis_mat), intent(in) :: monoMAT
    !> [in,out] 解ベクトル
    complex(kdouble), intent(inout) :: X(:)
    !> [in] 右辺ベクトル
    complex(kdouble), intent(in) :: B(:)
    !> [out] 残差ベクトル
    complex(kdouble), intent(out) :: R(:)
    !> [in,out] 計算時間
    real(kdouble), intent(inout) :: tspmv
    !> [in,out] 通信時間
    real(kdouble), intent(inout) :: tcomm
    integer(kint) :: i

    call monolis_std_debug_log_header("monolis_residual_main_C")

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